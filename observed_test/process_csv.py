import csv
import os

# --- Configuration ---
INPUT_FILENAME = 'spectra.csv'
OUTPUT_FILENAME = 'spectra_processed.csv'
# ---------------------

def process_csv(input_file, output_file):
    """
    Reads a CSV with 'id' and 'snr' columns, adds a 'file' column
    (id + '.fits'), and writes to a new CSV.
    """
    print(f"Starting to process '{input_file}'...")
    
    try:
        # Check if input file exists
        if not os.path.exists(input_file):
            print(f"Error: Input file '{input_file}' not found.")
            print("Please make sure the file is in the same directory as the script.")
            return

        with open(input_file, mode='r', encoding='utf-8') as infile:
            # Use DictReader to read based on column names
            reader = csv.DictReader(infile)
            
            # Get the fieldnames from the reader
            # We want to preserve the original order and add 'file'
            original_fieldnames = reader.fieldnames
            if 'id' not in original_fieldnames or 'snr' not in original_fieldnames:
                print(f"Error: Input file must have 'id' and 'snr' columns.")
                return

            # Define the new fieldnames
            # This inserts 'file' right after 'id'
            new_fieldnames = []
            for name in original_fieldnames:
                new_fieldnames.append(name)
                if name == 'id':
                    new_fieldnames.append('file')
            
            # If 'id' was the last column, just append 'file'
            if 'file' not in new_fieldnames:
                 new_fieldnames.append('file')


            with open(output_file, mode='w', encoding='utf-8', newline='') as outfile:
                # Use DictWriter for easy writing
                writer = csv.DictWriter(outfile, fieldnames=new_fieldnames)
                
                # Write the new header
                writer.writeheader()
                
                # Process each row
                for row in reader:
                    # Create the new 'file' value
                    file_id = row['id']
                    row['file'] = f"{file_id}.fits"
                    
                    # Write the modified row to the new file
                    writer.writerow(row)

        print(f"Successfully processed data and saved to '{output_file}'.")

    except Exception as e:
        print(f"An unexpected error occurred: {e}")

# --- Main execution ---
if __name__ == "__main__":
    process_csv(INPUT_FILENAME, OUTPUT_FILENAME)