import pandas as pd

def process_tsv(input_file, output_file):
    # Load the TSV into a DataFrame
    df = pd.read_csv(input_file, sep='\t')

    # Filter rows where the third column (z) equals 5
    df_filtered = df[df['z'] == -10].copy()

    # Set the third column (z) to 0
    df_filtered['z'] = 0

    # Save the processed DataFrame to a new TSV file
    df_filtered.to_csv(output_file, sep='\t', index=False)

# Example usage
input_file = 'setup_ecm.csv'  # Replace with your actual input file path
output_file = 'ecm_2_walls_cilindro_2d.tsv'  # Replace with the desired output file path

process_tsv(input_file, output_file)

