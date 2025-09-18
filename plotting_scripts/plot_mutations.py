# this script will read the .mut files from the output directory and plot the mutations

import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

# -----------------------
# Matplotlib global settings for publication quality
# -----------------------
mpl.rcParams.update({
    "figure.dpi": 300,        # High-resolution figures
    "savefig.dpi": 300,       # High-resolution saved images
    "pdf.fonttype": 42,       # Embed fonts as TrueType in PDF/PS
    "ps.fonttype": 42,
    "font.size": 10           # Base font size
})

output_dir = "output"

# Helper to extract time index from filename
def extract_time(filename):
    match = re.search(r'output(\d+)_mut', filename)
    return int(match.group(1)) if match else -1

# Get and sort all .mut files by time
mut_files = [f for f in os.listdir(output_dir) if f.endswith('mut.csv')]
mut_files = sorted(mut_files, key=extract_time)

# Build a list of DataFrames with time info
records = []
for mut_file in mut_files:
    time_idx = extract_time(mut_file)
    df = pd.read_csv(os.path.join(output_dir, mut_file))
    df['mutations'] = df['mutations'].fillna('')
    print(f"File: {mut_file}")
    print(df.head())
    for lineage, count in df['mutations'].value_counts().items():
        records.append({'time': time_idx, 'lineage': lineage, 'count': count})

print("Records:", records)

# Combine into a DataFrame
lineage_df = pd.DataFrame(records)

print(lineage_df.head())

if not lineage_df.empty:
    # Pivot: rows=time, columns=lineage, values=count (fill missing with 0)
    pivot_df = lineage_df.pivot_table(index='time', columns='lineage', values='count', fill_value=0)

    # Sort by time
    pivot_df = pivot_df.sort_index()

    # Ensure empty lineage label shows up in legend
    if '' in pivot_df.columns:
        pivot_df = pivot_df.rename(columns={'': 'No_mutations'})

    # Compute total cell count per lineage to rank abundance
    lineage_totals = pivot_df.sum().sort_values(ascending=False)
    top_n = 8  # number of dominant lineages to show in the main legend
    majority_lineages = list(lineage_totals.head(top_n).index)

    # Plot as stacked area chart
    fig, ax = plt.subplots(figsize=(10, 5))
    pivot_df.plot.area(ax=ax, colormap="tab20")

    # Axis labels and title with balanced font sizes
    ax.set_xlabel("Time Step", fontsize=12)
    ax.set_ylabel("Number of Cells", fontsize=12)
    ax.set_title("Lineage Abundance Over Time", fontsize=14)

    # Gather handles & labels after plotting
    handles, labels = ax.get_legend_handles_labels()

    # Map labels to handles for easy lookup
    label_to_handle = dict(zip(labels, handles))

    # Prepare legend entries for majority lineages
    handles_major = [label_to_handle[l] for l in majority_lineages if l in label_to_handle]
    labels_major = [l for l in majority_lineages if l in label_to_handle]

    # Add compact legend for majority lineages to the main figure
    ax.legend(handles_major, labels_major, title="Lineage (top)", bbox_to_anchor=(1.02, 1),
              loc="upper left", borderaxespad=0, fontsize=8, title_fontsize=9)

    # Prepare base path for saving outputs once, use for both main plot and legend
    base_path = os.path.join(os.path.dirname(__file__), "mutations")

    # -----------------------------
    # Save full legend as a separate figure (use fresh patch handles to avoid artist duplication)
    # Build new handles with the same facecolors
    patch_handles = []
    for h, lab in zip(handles, labels):
        # Extract facecolor from original collection/patch
        if hasattr(h, 'get_facecolor'):
            fc = h.get_facecolor()
            # facecolor returns an array-like; take the first entry if needed
            color = fc[0] if len(fc) else fc
        elif hasattr(h, 'get_color'):
            color = h.get_color()
        else:
            color = 'gray'
        patch_handles.append(mpl.patches.Patch(facecolor=color, label=lab))

    # Estimate height (0.3 inch per entry)
    fig_leg = plt.figure(figsize=(3, max(1, 0.3 * len(patch_handles))))
    fig_leg.legend(handles=patch_handles, labels=[ph.get_label() for ph in patch_handles],
                   ncol=1, frameon=False, loc='center', title='Lineage')
    fig_leg.tight_layout()
    legend_png = f"{base_path}_legend.png"
    legend_pdf = f"{base_path}_legend.pdf"
    fig_leg.savefig(legend_png, bbox_inches='tight')
    fig_leg.savefig(legend_pdf, bbox_inches='tight')
    plt.close(fig_leg)

    # Tight layout to ensure everything fits nicely
    fig.tight_layout()

    # Save the figure (both raster and vector formats) next to this script
    png_path = f"{base_path}.png"
    pdf_path = f"{base_path}.pdf"

    fig.savefig(png_path, bbox_inches="tight")  # 300 DPI from rcParams
    fig.savefig(pdf_path, bbox_inches="tight")  # Vector format for publications

    # --------------------------------------------------
    # Normalized (100%) stacked AREA plot of lineage proportions (time on X, proportions on Y)
    # --------------------------------------------------
    # Compute percentage per time step
    pivot_pct = pivot_df.div(pivot_df.sum(axis=1), axis=0) * 100

    # Build a color mapping from the original handles
    lineage_colors = {}
    for lbl, hndl in label_to_handle.items():
        if hasattr(hndl, "get_facecolor"):
            fc = hndl.get_facecolor()
            lineage_colors[lbl] = fc[0] if len(fc) else fc
        elif hasattr(hndl, "get_color"):
            lineage_colors[lbl] = hndl.get_color()

    # Color list matching column order
    color_list = [lineage_colors.get(col, None) for col in pivot_pct.columns]

    fig_prop, ax_prop = plt.subplots(figsize=(10, 5))
    pivot_pct.plot.area(ax=ax_prop, color=color_list)

    # Axis formatting – journal style
    ax_prop.set_xlabel("Time Step", fontsize=12)
    ax_prop.set_ylabel("Lineage Proportion (%)", fontsize=12)
    ax_prop.set_title("Normalized Lineage Proportions Over Time", fontsize=14, pad=12)
    ax_prop.set_ylim(0, 100)

    # Clean aesthetic: remove extra spines and add subtle grid
    for spine in ["top", "right", "left", "bottom"]:
        ax_prop.spines[spine].set_visible(False)
    ax_prop.grid(axis="y", color="0.9", linewidth=0.7)

    # Legend for majority lineages (fresh Patch handles)
    prop_patch_handles = [
        mpl.patches.Patch(facecolor=lineage_colors.get(l, "gray"), label=l)
        for l in majority_lineages if l in lineage_colors
    ]
    ax_prop.legend(handles=prop_patch_handles,
                   labels=[p.get_label() for p in prop_patch_handles],
                   title="Lineage (top)", bbox_to_anchor=(1.02, 1), loc="upper left",
                   borderaxespad=0, fontsize=8, title_fontsize=9)

    fig_prop.tight_layout()
    fig_prop.savefig(f"{base_path}_proportion.png", bbox_inches="tight")
    fig_prop.savefig(f"{base_path}_proportion.pdf", bbox_inches="tight")
    plt.close(fig_prop)

    print(f"Figures saved to: {png_path} and {pdf_path}")
else:
    print("No lineage data found. Check your .mut.csv files and their contents.")

# Suppose you have your pivot_df as before
# You need to assign a unique integer ID to each lineage
lineage_ids = {lin: i for i, lin in enumerate(pivot_df.columns)}
populations = []
for time, row in pivot_df.iterrows():
    for lin, count in row.items():
        populations.append({'Id': lineage_ids[lin], 'Step': time, 'Pop': count})
pop_df = pd.DataFrame(populations)

# If you have parent-child relationships, build a parent_df as well
# For now, if you don't, you can set all parents to 0 or None
parent_df = pd.DataFrame({'ParentId': [0]*len(lineage_ids), 'ChildId': list(lineage_ids.values())})

# Plot
# plot_fish(pop_df, parent_df)


