import pandas as pd

# Load the CSV file
df = pd.read_csv("cilindro_cells.csv")

# Add a new column filled with 0
df["id"] = 0

# Save the modified CSV
df.to_csv("your_file_modified.csv", index=False)
