import os
import cv2
import numpy as np
from torch.utils.data import Dataset
import matplotlib.pyplot as plt
import torch
import multiprocessing

from PIL import Image
from torch.utils.data import DataLoader
from multiprocessing import freeze_support

class YoloDataset_original(Dataset):
    """YOLO Dataset. Read images, apply augmentation and preprocessing transformations.
    
    Args:
        data_dir (str): path to the directory containing images and YOLO annotations
        classes (list): list of class names
        augmentation (albumentations.Compose): data transformation pipeline 
            (e.g. flip, scale, etc.)
        preprocessing (albumentations.Compose): data preprocessing 
            (e.g. normalization, shape manipulation, etc.)
    """
    
    def __init__(self, data_dir, classes=None, mode="train", transform=None):
        self.data_dir = data_dir
        self.images_dir = os.path.join(data_dir, mode, "images")
        self.labels_dir = os.path.join(data_dir, mode, "labels")
        self.ids = [os.path.splitext(file)[0] for file in os.listdir(self.images_dir) if file.endswith('.jpg')]
        self.images_fps = [os.path.join(self.images_dir, image_id + '.jpg') for image_id in self.ids]
        self.annotations_fps = [os.path.join(self.labels_dir, image_id + '.txt') for image_id in self.ids]
        self.classes = classes
        self.class_values = {cls: idx + 1 for idx, cls in enumerate(classes)}  # class IDs start from 1
        
        self.mode = mode
        self.transform = transform
    
    def __len__(self):
        return len(self.ids)

    def __getitem__(self, idx):
        filename = self.ids[idx]
        image_path = os.path.join(self.images_dir, filename + ".jpg")
        annotation_path = os.path.join(self.labels_dir, filename + ".txt")
        image = np.array(Image.open(image_path).convert("RGB"))
        h, w, _ = image.shape
        mask = self._read_annotation(annotation_path, w, h)
        
        result = dict(image=image, mask=mask)
        if self.transform is not None:
            result = self.transform(**result)
        
        # resize images
        image = np.array(
            Image.fromarray(result["image"]).resize((512, 512), Image.BILINEAR)
        )
        if result["mask"] is not None:
            mask = np.array(
                Image.fromarray(result["mask"]).resize((512, 512), Image.NEAREST)
            )
            result["mask"] = np.expand_dims(mask, 0)
        else:
            result["mask"] = np.zeros((1, 512, 512), dtype=np.uint8)

        # convert to other format HWC -> CHW
        result["image"] = np.moveaxis(image, -1, 0)
            
        return result

    def _read_annotation(self, annotation_path, img_width, img_height):
        if not os.path.exists(annotation_path):
            return None
        
        mask = np.zeros((img_height, img_width), dtype=np.uint8)
        with open(annotation_path, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) < 3:
                    print(f"Skipping invalid annotation line: {line}")
                    continue
                
                class_id = int(parts[0]) + 1
                coords = list(map(float, parts[1:]))
                if len(coords) % 2 != 0:
                    print(f"Skipping invalid coordinate line: {line}")
                    continue
                
                coords = np.array(coords).reshape(-1, 2)
                
                # Denormalize the coordinates
                coords[:, 0] *= img_width
                coords[:, 1] *= img_height
                
                # Create a mask from the polygon
                coords = coords.astype(np.int32)
                cv2.fillPoly(mask, [coords], class_id)
        

        return mask

import os
import torch
import torch.optim as optim
import pytorch_lightning as pl
import segmentation_models_pytorch as smp
from torch.utils.data import DataLoader
from torch.cuda.amp import GradScaler, autocast
from torchmetrics.classification import BinaryPrecision, BinaryRecall

class SegModel(pl.LightningModule):
    def __init__(self, arch, encoder_name, encoder_weights, in_channels, out_classes, **kwargs):
        super().__init__()
        self.save_hyperparameters
        self.model = smp.create_model(
            arch, encoder_name=encoder_name, encoder_weights=encoder_weights, in_channels=in_channels, classes=out_classes, **kwargs
        )

        # preprocessing parameters for image
        params = smp.encoders.get_preprocessing_params(encoder_name, encoder_weights)
        self.register_buffer("std", torch.tensor(params["std"]).view(1, 3, 1, 1))
        self.register_buffer("mean", torch.tensor(params["mean"]).view(1, 3, 1, 1))
        self.training_step_outputs = []
        self.validing_step_outputs = []
        self.testing_step_outputs = []

        # for image segmentation dice loss could be the best first choice
        self.loss_fn = smp.losses.DiceLoss(smp.losses.BINARY_MODE, from_logits=True)
        
        # Metrics
        self.precision = BinaryPrecision(threshold=0.5)
        self.recall = BinaryRecall(threshold=0.5)

    def forward(self, image):
        # normalize image here
        image = (image - self.mean) / self.std
        mask = self.model(image)
        return mask

    def training_step(self, batch, batch_idx):
        # training_step defines the train loop.
        # it is independent of forward
        
        image = batch["image"]
        # Shape of the image should be (batch_size, num_channels, height, width)
        # if you work with grayscale images, expand channels dim to have [batch_size, 1, height, width]
        assert image.ndim == 4
        # Check that image dimensions are divisible by 32, 
        h, w = image.shape[2:]
        assert h % 32 == 0 and w % 32 == 0

        mask = batch["mask"]
        # Shape of the mask should be [batch_size, num_classes, height, width]
        # for binary segmentation num_classes = 1
        assert mask.ndim == 4

        # Check that mask values in between 0 and 1, NOT 0 and 255 for binary segmentation
        assert mask.max() <= 1.0 and mask.min() >= 0
        
        logits_mask = self.forward(image)
        
        # Predicted mask contains logits, and loss_fn param `from_logits` is set to True
        loss = self.loss_fn(logits_mask, mask)
        # Lets compute metrics for some threshold
        # first convert mask values to probabilities, then 
        # apply thresholding
        prob_mask = logits_mask.sigmoid()
        pred_mask = (prob_mask > 0.5).float()
        # We will compute IoU metric by two ways
        #   1. dataset-wise
        #   2. image-wise
        # but for now we just compute true positive, false positive, false negative and
        # true negative 'pixels' for each image and class
        # these values will be aggregated in the end of an epoch
        tp, fp, fn, tn = smp.metrics.get_stats(pred_mask.long(), mask.long(), mode="binary")
        self.training_step_outputs.append({
            "loss": loss.detach(),  # Ensure loss is detached
            "tp": tp,
            "fp": fp,
            "fn": fn,
            "tn": tn,
            "pred_mask": pred_mask,
            "mask": mask
        })
        return loss

    def on_train_epoch_end(self):
        tp = torch.cat([x["tp"] for x in self.training_step_outputs])
        fp = torch.cat([x["fp"] for x in self.training_step_outputs])
        fn = torch.cat([x["fn"] for x in self.training_step_outputs])
        tn = torch.cat([x["tn"] for x in self.training_step_outputs])
        pred_masks = torch.cat([x["pred_mask"] for x in self.training_step_outputs])
        masks = torch.cat([x["mask"] for x in self.training_step_outputs])

        # per image IoU means that we first calculate IoU score for each image 
        # and then compute mean over these scores
        per_image_iou = smp.metrics.iou_score(tp, fp, fn, tn, reduction="micro-imagewise")
        
        # dataset IoU means that we aggregate intersection and union over whole dataset
        # and then compute IoU score. The difference between dataset_iou and per_image_iou scores
        # in this particular case will not be much, however for dataset 
        # with "empty" images (images without target class) a large gap could be observed. 
        # Empty images influence a lot on per_image_iou and much less on dataset_iou.
        dataset_iou = smp.metrics.iou_score(tp, fp, fn, tn, reduction="micro")

        # Calculate precision and recall
        precision = self.precision(pred_masks, masks)
        recall = self.recall(pred_masks, masks)

        # Calculate mAP@50
        map50 = (precision * recall).mean().item()

        # Log the metrics
        avg_loss = torch.stack([x["loss"] for x in self.training_step_outputs]).mean()
        self.log('tr_loss', avg_loss, prog_bar=True,on_epoch=True)
        self.log('tr_image_iou', per_image_iou, prog_bar=True,on_epoch=True)
        self.log('tr_dataset_iou', dataset_iou, prog_bar=True,on_epoch=True)
        self.log('tr_map50', map50, prog_bar=True,on_epoch=True)
        self.training_step_outputs.clear()  # free memory

    def validation_step(self, batch, batch_idx):
        # training_step defines the train loop.
        # it is independent of forward
        image = batch["image"]

        # Shape of the image should be (batch_size, num_channels, height, width)
        # if you work with grayscale images, expand channels dim to have [batch_size, 1, height, width]
        assert image.ndim == 4
        # Check that image dimensions are divisible by 32, 
        h, w = image.shape[2:]
        assert h % 32 == 0 and w % 32 == 0

        mask = batch["mask"]
        # Shape of the mask should be [batch_size, num_classes, height, width]
        # for binary segmentation num_classes = 1
        assert mask.ndim == 4
        
        # Check that mask values in between 0 and 1, NOT 0 and 255 for binary segmentation
        assert mask.max() <= 1.0 and mask.min() >= 0
              

        logits_mask = self.forward(image)
           

        # Predicted mask contains logits, and loss_fn param `from_logits` is set to True
        loss = self.loss_fn(logits_mask, mask)


        # Lets compute metrics for some threshold
        # first convert mask values to probabilities, then 
        # apply thresholding
        prob_mask = logits_mask.sigmoid()
        pred_mask = (prob_mask > 0.5).float()

        # We will compute IoU metric by two ways
        #   1. dataset-wise
        #   2. image-wise
        # but for now we just compute true positive, false positive, false negative and
        # true negative 'pixels' for each image and class
        # these values will be aggregated in the end of an epoch
        tp, fp, fn, tn = smp.metrics.get_stats(pred_mask.long(), mask.long(), mode="binary")

        # Calculate precision and recall
        precision = self.precision(pred_mask, mask)
        recall = self.recall(pred_mask, mask)

        # Calculate mAP@50
        map50 = (precision * recall).mean().item()        
        self.log('val_loss', loss, prog_bar=True,on_epoch=True)
        self.log('val_map50', map50, prog_bar=True,on_epoch=True)
        
    def test_step(self, batch, batch_idx):
        # training_step defines the train loop.
        # it is independent of forward
        image = batch["image"]

        # Shape of the image should be (batch_size, num_channels, height, width)
        # if you work with grayscale images, expand channels dim to have [batch_size, 1, height, width]
        assert image.ndim == 4
        # Check that image dimensions are divisible by 32, 
        h, w = image.shape[2:]
        assert h % 32 == 0 and w % 32 == 0

        mask = batch["mask"]
        # Shape of the mask should be [batch_size, num_classes, height, width]
        # for binary segmentation num_classes = 1
        assert mask.ndim == 4
        
        # Check that mask values in between 0 and 1, NOT 0 and 255 for binary segmentation
        assert mask.max() <= 1.0 and mask.min() >= 0
              

        logits_mask = self.forward(image)
           

        # Predicted mask contains logits, and loss_fn param `from_logits` is set to True
        loss = self.loss_fn(logits_mask, mask)

        # Lets compute metrics for some threshold
        # first convert mask values to probabilities, then 
        # apply thresholding
        prob_mask = logits_mask.sigmoid()
        pred_mask = (prob_mask > 0.5).float()

        # We will compute IoU metric by two ways
        #   1. dataset-wise
        #   2. image-wise
        # but for now we just compute true positive, false positive, false negative and
        # true negative 'pixels' for each image and class
        # these values will be aggregated in the end of an epoch
        tp, fp, fn, tn = smp.metrics.get_stats(pred_mask.long(), mask.long(), mode="binary")

        # Calculate precision and recall
        precision = self.precision(pred_mask, mask)
        recall = self.recall(pred_mask, mask)

        # Calculate mAP@50
        map50 = (precision * recall).mean().item()        
        self.log('test_loss', loss, prog_bar=True,on_epoch=True)
        self.log('test_map50', map50, prog_bar=True,on_epoch=True)


    def configure_optimizers(self):
        optimizer = optim.Adam(self.parameters(), lr=1e-4)
        return optimizer


    def train_dataloader(self):
        return DataLoader(
            self.train_dataset,
            batch_size=16,
            shuffle=True,
            num_workers=os.cpu_count(),
            persistent_workers=True,
            drop_last=True  # Ensures each batch is the same size
        )

    def val_dataloader(self):
        return DataLoader(
            self.valid_dataset,
            batch_size=16,
            shuffle=False,
            num_workers=os.cpu_count(),
            persistent_workers=True,
            drop_last=False
        )


    def test_dataloader(self):
        return DataLoader(
            self.test_dataset,
            batch_size=16,
            shuffle=False,
            num_workers=os.cpu_count(),
            persistent_workers=True,
            drop_last=False
        )

def main():
    DATA_DIR = './DRG123_batch2/'
    image_height, image_width = 512, 512
    CLASSES = ['cell']

    # Read combinations from the text file
    with open('combinations.txt', 'r') as f:
        combinations = [line.strip() for line in f if line.strip()]

    for combination in combinations:
        arch, encoder_name, encoder_weights = combination.split(',')
        arch = arch.strip()
        encoder_name = encoder_name.strip()
        encoder_weights = encoder_weights.strip()
        model_name = f"{arch}_result_{encoder_name}"
        

        print(f"Training model with architecture {arch} and encoder {encoder_name}")

        # Initialize train, val, test sets
        train_dataset = YoloDataset_original(
            DATA_DIR,
            mode="train",
            classes=CLASSES,
        )

        valid_dataset = YoloDataset_original(
            DATA_DIR,
            mode="val",
            classes=CLASSES,
        )

        test_dataset = YoloDataset_original(
            DATA_DIR,
            mode="test",
            classes=CLASSES,
        )

        print(f"Train size: {len(train_dataset)}")
        print(f"Valid size: {len(valid_dataset)}")
        print(f"Test size: {len(test_dataset)}")
        n_cpu = os.cpu_count()

        train_loader = DataLoader(train_dataset, batch_size=16, shuffle=True)
        valid_loader = DataLoader(valid_dataset, batch_size=16, shuffle=False)
        test_loader = DataLoader(test_dataset, batch_size=16, shuffle=False)

        model = SegModel(arch, encoder_name, encoder_weights, in_channels=3, out_classes=1, )

        # Create unique logger and CSV logger for each combination
        logger = pl.loggers.TensorBoardLogger(save_dir=model_name, name='tb_logs')
        csv_logger = pl.loggers.CSVLogger(save_dir=model_name, name='csv_logs')

        trainer = pl.Trainer(
            default_root_dir=model_name,
            max_epochs=300,
            logger=[logger, csv_logger]
        )
        trainer.fit(model, train_loader, valid_loader)

        valid_metrics = trainer.validate(model, dataloaders=valid_loader, verbose=False)
        print(valid_metrics)

        test_metrics = trainer.test(model, dataloaders=test_loader, verbose=False)
        print(test_metrics)
        
if __name__ == '__main__':
    import multiprocessing
    multiprocessing.freeze_support()
    main()