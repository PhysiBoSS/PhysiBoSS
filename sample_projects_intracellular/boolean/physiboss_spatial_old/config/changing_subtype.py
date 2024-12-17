import pandas as pd
import random

# Load the CSV file into a DataFrame
file_path = 'cells_prisma.csv'
df = pd.read_csv(file_path)

# Calculate the number of rows where 'id' is 0
num_zeros = (df['id'] == 0).sum()

# Determine how many zeros to change to ones (50% of the total zeros)
num_to_change = int(num_zeros * 0.5)

# Get indices of rows where 'id' is 0
zero_indices = df.index[df['id'] == 0].tolist()

# Randomly select indices to change to 1
indices_to_change = random.sample(zero_indices, num_to_change)

# Change the selected rows to 1
df.loc[indices_to_change, 'id'] = 2

# Save the modified DataFrame back to a new CSV file
output_file = 'cells_prisma_modified.csv'
df.to_csv(output_file, index=False)

print(f"Successfully modified {num_to_change} '0's to '1's and saved to {output_file}.")

