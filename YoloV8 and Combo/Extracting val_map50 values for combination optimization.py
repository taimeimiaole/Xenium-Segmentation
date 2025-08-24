import os
import pandas as pd

# List of CSV file paths (use raw string literals for Windows paths)
csv_files = [
    r'path to metrics.csv for each combination'

# List to store the extracted dataframes
extracted_dfs = []

# List to store error information
error_files = []

# Loop through each file and extract the val_map50 value
for file in csv_files:
    try:
        print(f"Processing file: {file}")
        df = pd.read_csv(file)
        print(f"Columns in the file: {df.columns.tolist()}")
        if 'val_map50' in df.columns:
            # Use the full file path as the column name
            column_name = file
            print(f"Extracting 'val_map50' values from {file}")
            extracted_df = df[['val_map50']].copy()
            extracted_df.columns = [column_name]
            extracted_dfs.append(extracted_df)
        else:
            print(f"'val_map50' column not found in {file}")
            error_files.append((file, 'val_map50 column not found'))
    except pd.errors.EmptyDataError:
        print(f"EmptyDataError: {file} is empty or has no data.")
        error_files.append((file, 'EmptyDataError'))
    except pd.errors.ParserError:
        print(f"ParserError: {file} is corrupted or has parsing issues.")
        error_files.append((file, 'ParserError'))
    except FileNotFoundError:
        print(f"FileNotFoundError: {file} not found.")
        error_files.append((file, 'FileNotFoundError'))
    except Exception as e:
        print(f"Unexpected error processing {file}: {e}")
        error_files.append((file, str(e)))

# Concatenate all the extracted dataframes
if extracted_dfs:
    output_df = pd.concat(extracted_dfs, axis=1)
else:
    output_df = pd.DataFrame()

print("Output DataFrame:")
print(output_df)

# Write the DataFrame to a new CSV file
output_df.to_csv('name of output csv file', index=False)

# Save error log to a separate CSV file
error_df = pd.DataFrame(error_files, columns=['File', 'Error'])
error_df.to_csv('error_log.csv', index=False)

print("Extraction complete. Data saved to 'extracted_val_map50_values_secondcombos.csv'")
print("Error log saved to 'error_log.csv'")
