import csv

input_file_path = r"C:\Users\phch\OneDrive - Danmarks Tekniske Universitet\Mapper\Python_data\MolFormatic\worklist\procul\missing_trans.txt"
output_file_path = r"C:\Users\phch\OneDrive - Danmarks Tekniske Universitet\Mapper\Python_data\MolFormatic\worklist\procul\missing_trans_done_2.txt"

# Default values
default_procul = 'procul_001'
default_source = 'source_2'

# Open the input and output files
with open(input_file_path, 'r') as infile, open(output_file_path, 'w', newline='') as outfile:
    reader = infile.readlines()
    csv_writer = csv.writer(outfile, delimiter=';')

    # Write header to the CSV file
    csv_writer.writerow(["destination_plates", "destination_well", "volume", "source_well", "source_plates"])

    # Process each line in the input file
    for row, line in enumerate(reader):
        if row < 3 :
            continue

        # Split the line based on spaces
        parts = line.split()
        print(parts)
        # Extract relevant information
        destination = parts[4]
        print(f"destination - {destination}")
        volume = parts[6]
        print(f"volume - {volume}")
        source = "O6"
        print(f"source - {source}")
        # Write the extracted information to the CSV file
        csv_writer.writerow([default_procul, destination, volume, source, default_source])

print(f"Conversion complete. Output saved to {output_file_path}")
