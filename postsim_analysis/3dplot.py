import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from ipywidgets import widgets


df = pd.read_csv("./total_info_sim.csv")


# Simple 3D scatter

# fig = px.scatter_3d(
#     df,
#     x='position_x',
#     y='position_y',
#     z='time',
#     color='drug_X_internal_density',
#     # size='drug_Y_internal_density',  # Use 'A' to determine the size of the points
#     symbol='pi3k_node',  # Use 'Category' to determine the symbol of the points
#     hover_name='ID',  # Show 'Category' name when hovering
#     title='Interactive 3D Scatter Plot with Customization'
#     # labels={'A': 'X Axis', 'B': 'Y Axis', 'C': 'Z Axis'}
# )

# # Show plot in browser
# fig.show()


# Function to update the scatter plot based on slider value
def update_scatter(z_multiplier):
    fig = go.Figure(data=[go.Scatter3d(
        x=df['position_x'],
        y=df['position_y'],
        z=df['time'] * z_multiplier,
        mode='markers',
        marker=dict(size=10, color=df['pi3k_node'], colorscale='Viridis'),
        text=df['ID']
    )])
    
    fig.update_layout(
        title=f'Interactive 3D Scatter Plot (Z scaled by {z_multiplier})',
        scene=dict(
            xaxis_title='pos X',
            yaxis_title='pos Y',
            zaxis_title='time (scaled)'
        )
    )
    fig.show()

# Create a slider widget
slider = widgets.FloatSlider(
    value=1.0,
    min=4200,
    max=4200.0,
    step=40,
    description='Z scale:'
)

# Display the interactive plot with slider
widgets.interact(update_scatter, z_multiplier=slider)