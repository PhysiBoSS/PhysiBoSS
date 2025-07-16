import pandas as pd

# Load the CSV
df = pd.read_csv("cilindro_cells.csv", header=None)  # Replace with your actual file name

# Get the last column (by index)
last_col = df.columns[-1]

# Count where last column is 0 or 2
count_0 = (df[last_col] == 0).sum()
count_2 = (df[last_col] == 2).sum()

print(f"Rows where last column is 0: {count_0}")
print(f"Rows where last column is 2: {count_2}")

