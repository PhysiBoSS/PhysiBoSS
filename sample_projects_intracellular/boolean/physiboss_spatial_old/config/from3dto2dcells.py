import pandas as pd

def process_csv(input_file, output_file):
    # Load the CSV into a DataFrame
    df = pd.read_csv(input_file, header=None)

    # Filter the rows where the third column is 10
    df_filtered = df[df[2] == -10].copy()  # Use .copy() to avoid the SettingWithCopyWarning

    # Set the 3rd column to 0 for the remaining rows
    df_filtered.loc[:, 2] = 0  # Use .loc[] to safely modify the DataFrame

    # Save the filtered DataFrame to a new CSV file
    df_filtered.to_csv(output_file, header=False, index=False)


# Example usage
input_file = 'setup_cells_norm.csv'  # Replace with your actual input file path
output_file = 'cilindro_cells2d.csv'  # Replace with the desired output file path

process_csv(input_file, output_file)

