import pandas as pd

# Cargar el archivo CSV
df = pd.read_csv('cilindro_cells_original.csv')

# Reemplazar la última columna por 0 (como entero)
df.iloc[:, -1] = 0

# Guardar el nuevo archivo CSV
df.to_csv('cilindro_cells.csv', index=False)

