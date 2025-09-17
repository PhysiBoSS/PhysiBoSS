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
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "font.size": 10
})

output_dir = "output"

# Extract time index from filename
def extract_time(filename):
    match = re.search(r'output(\d+)_mut', filename)
    return int(match.group(1)) if match else -1

# Transform full mutation string into final genotype
def final_genotype(mutation_str):
    if not mutation_str or pd.isna(mutation_str):
        return 'No_mutations'
    mutations = mutation_str.split('.')
    gene_state = {}
    for mut in mutations:
        if '_' in mut:
            gene, state = mut.split('_')
            gene_state[gene] = state  # keep latest mutation per gene
    # Sort genes alphabetically for consistency
    final = [f"{gene}_{gene_state[gene]}" for gene in sorted(gene_state.keys())]
    return '.'.join(final)

# Read all .mut.csv files
mut_files = [f for f in os.listdir(output_dir) if f.endswith('mut.csv')]
mut_files = sorted(mut_files, key=extract_time)

records = []
for mut_file in mut_files:
    time_idx = extract_time(mut_file)
    df = pd.read_csv(os.path.join(output_dir, mut_file))
    df['mutations'] = df['mutations'].fillna('')
    df['final_lineage'] = df['mutations'].apply(final_genotype)

    print(f"File: {mut_file}")
    print(df[['mutations', 'final_lineage']].head())

    for lineage, count in df['final_lineage'].value_counts().items():
        records.append({'time': time_idx, 'lineage': lineage, 'count': count})

# Combine into DataFrame
lineage_df = pd.DataFrame(records)
print(lineage_df.head())

if not lineage_df.empty:
    # Pivot: time x lineage
    pivot_df = lineage_df.pivot_table(index='time', columns='lineage', values='count', fill_value=0)
    pivot_df = pivot_df.sort_index()

    # Rename empty lineage
    if 'No_mutations' in pivot_df.columns:
        pivot_df = pivot_df.rename(columns={'No_mutations': 'No_mutations'})

    # Sort by total abundance
    lineage_totals = pivot_df.sum().sort_values(ascending=False)
    top_n = 8
    majority_lineages = list(lineage_totals.head(top_n).index)

    # Plot
    fig, ax = plt.subplots(figsize=(10, 5))
    pivot_df.plot.area(ax=ax, colormap="tab20")
    ax.set_xlabel("Time Step", fontsize=12)
    ax.set_ylabel("Number of Cells", fontsize=12)
    ax.set_title("Lineage Abundance Over Time", fontsize=14)

    # Legend
    handles, labels = ax.get_legend_handles_labels()
    label_to_handle = dict(zip(labels, handles))
    handles_major = [label_to_handle[l] for l in majority_lineages if l in label_to_handle]
    labels_major = [l for l in majority_lineages if l in label_to_handle]

    ax.legend(handles_major, labels_major, title="Lineage (top)", bbox_to_anchor=(1.02, 1),
              loc="upper left", borderaxespad=0, fontsize=8, title_fontsize=9)

    # Save paths
    base_path = os.path.join(os.path.dirname(__file__), "mutations")

    # Full legend
    patch_handles = []
    for h, lab in zip(handles, labels):
        if hasattr(h, 'get_facecolor'):
            fc = h.get_facecolor()
            color = fc[0] if len(fc) else fc
        elif hasattr(h, 'get_color'):
            color = h.get_color()
        else:
            color = 'gray'
        patch_handles.append(mpl.patches.Patch(facecolor=color, label=lab))

    fig_leg = plt.figure(figsize=(3, max(1, 0.3 * len(patch_handles))))
    fig_leg.legend(handles=patch_handles, labels=[ph.get_label() for ph in patch_handles],
                   ncol=1, frameon=False, loc='center', title='Lineage')
    fig_leg.tight_layout()
    fig_leg.savefig(f"{base_path}_legend.png", bbox_inches='tight')
    fig_leg.savefig(f"{base_path}_legend.pdf", bbox_inches='tight')
    plt.close(fig_leg)

    # Save main figure
    fig.tight_layout()
    fig.savefig(f"{base_path}_aggregated.png", bbox_inches="tight")
    fig.savefig(f"{base_path}_aggregated.pdf", bbox_inches="tight")

    print(f"Figures saved to: {base_path}_aggregated.png and {base_path}_aggregated.pdf")

    # --------------------------------------------------
    # Normalized (100%) stacked AREA plot of lineage proportions (time on X, proportions on Y)
    # --------------------------------------------------
    # Compute percentage per time step
    pivot_pct = pivot_df.div(pivot_df.sum(axis=1), axis=0) * 100

    # Consistent colors: map each lineage to its area color (reuse fishplot colors)
    lineage_colors = {}
    for lbl, hndl in label_to_handle.items():
        if hasattr(hndl, "get_facecolor"):
            fc = hndl.get_facecolor()
            lineage_colors[lbl] = fc[0] if len(fc) else fc
        elif hasattr(hndl, "get_color"):
            lineage_colors[lbl] = hndl.get_color()

    # Define color list in column order for area plot
    color_list = [lineage_colors.get(col, None) for col in pivot_pct.columns]

    fig_prop, ax_prop = plt.subplots(figsize=(10, 5))
    pivot_pct.plot.area(ax=ax_prop, color=color_list)

    # Formatting – journal style
    ax_prop.set_xlabel("Time Step", fontsize=12)
    ax_prop.set_ylabel("Lineage Proportion (%)", fontsize=12)
    ax_prop.set_title("Normalized Lineage Proportions Over Time", fontsize=14, pad=12)
    ax_prop.set_ylim(0, 100)

    # Clean up spines for a minimalist look
    for spine in ["top", "right", "left", "bottom"]:
        ax_prop.spines[spine].set_visible(False)
    ax_prop.grid(axis="y", color="0.9", linewidth=0.7)

    # Legend outside – top lineages only (fresh Patch handles)
    prop_patch_handles = [
        mpl.patches.Patch(facecolor=lineage_colors.get(l, "gray"), label=l)
        for l in majority_lineages if l in lineage_colors
    ]
    ax_prop.legend(handles=prop_patch_handles,
                   labels=[p.get_label() for p in prop_patch_handles],
                   title="Lineage (top)", bbox_to_anchor=(1.02, 1), loc="upper left",
                   borderaxespad=0, fontsize=8, title_fontsize=9)

    fig_prop.tight_layout()
    fig_prop.savefig(f"{base_path}_aggregated_proportion.png", bbox_inches="tight")
    fig_prop.savefig(f"{base_path}_aggregated_proportion.pdf", bbox_inches="tight")
    plt.close(fig_prop)

else:
    print("No lineage data found. Check your .mut.csv files and their contents.")

# Build fishplot-compatible DataFrames
lineage_ids = {lin: i for i, lin in enumerate(pivot_df.columns)}
populations = []
for time, row in pivot_df.iterrows():
    for lin, count in row.items():
        populations.append({'Id': lineage_ids[lin], 'Step': time, 'Pop': count})
pop_df = pd.DataFrame(populations)

# Set all parent IDs to 0 for now
parent_df = pd.DataFrame({'ParentId': [0]*len(lineage_ids), 'ChildId': list(lineage_ids.values())})

# You can use pop_df and parent_df for fishplot if needed
# plot_fish(pop_df, parent_df)

