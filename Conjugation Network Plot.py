import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

# =========================
# Load data
# =========================
file_path = r""
df = pd.read_excel(file_path)

required_cols = ["Metal_form", "Concentration_mg_per_L", "Duration", "Value_Log10",
                 "Donor", "Recipient", "Plasmid_ID"]
missing = [c for c in required_cols if c not in df.columns]
if missing:
    raise ValueError(f"Missing columns: {missing}\nFound: {list(df.columns)}")

def clean_str(s: pd.Series) -> pd.Series:
    return (s.astype(str)
              .str.replace("\n", " ", regex=True)
              .str.replace("\r", " ", regex=True)
              .str.replace(r"\s+", " ", regex=True)
              .str.strip())

# =========================
# Clean + standardize
# =========================
# Forward fill ONLY non-metal fields 
df[["Duration", "Donor", "Recipient", "Plasmid_ID"]] = df[["Duration", "Donor", "Recipient", "Plasmid_ID"]].ffill()

for col in ["Metal_form", "Duration", "Donor", "Recipient", "Plasmid_ID"]:
    df[col] = clean_str(df[col])

# Treat blanks / nan-like strings in Metal_form as NA
df["Metal_form"] = df["Metal_form"].replace({"nan": "NA", "none": "NA", "None": "NA", "": "NA"}).fillna("NA")
df["Metal_form"] = df["Metal_form"].str.replace(r"\s+", " ", regex=True).str.strip()

# Numeric conversions
df["Value_Log10"] = clean_str(df["Value_Log10"]).str.replace("−", "-", regex=False)
df["Value_Log10"] = pd.to_numeric(df["Value_Log10"], errors="coerce")

df["Concentration_mg_per_L"] = clean_str(df["Concentration_mg_per_L"]).str.replace("−", "-", regex=False)
df["Concentration_mg_per_L"] = pd.to_numeric(df["Concentration_mg_per_L"], errors="coerce")

df = df.dropna(subset=["Value_Log10", "Concentration_mg_per_L"])

# =========================
# NA CONTROL definition:
# NA + 0 mg/L is the control, matched by Donor+Recipient+Plasmid+Duration
# =========================
control_keys = ["Donor", "Recipient", "Plasmid_ID", "Duration"]

control_df = df[(df["Metal_form"].str.upper() == "NA") & (df["Concentration_mg_per_L"] == 0)].copy()

# If duplicate NA controls exist per block, average them
control_df = (
    control_df.groupby(control_keys, as_index=False)["Value_Log10"]
              .mean()
              .rename(columns={"Value_Log10": "Control_Log10"})
)

# Attach control to all rows in that block (inner keeps only rows with valid NA control)
merged = pd.merge(df, control_df, on=control_keys, how="inner")

# Treated = everything except the NA/0 row itself
treated_df = merged[~((merged["Metal_form"].str.upper() == "NA") & (merged["Concentration_mg_per_L"] == 0))].copy()

# ΔLog10 for ALL metals vs NA control (within each block)
treated_df["Delta_Log10"] = treated_df["Value_Log10"] - treated_df["Control_Log10"]
treated_df["Effect"] = np.where(treated_df["Delta_Log10"] > 0, "Induced", "Inhibited")

# =========================
# KEEP ALL concentrations (ALL log changes)
# =========================
significant_df = treated_df.copy().reset_index(drop=True)

# Optional: remove exact duplicates
significant_df = significant_df.drop_duplicates(
    subset=["Metal_form", "Donor", "Recipient", "Plasmid_ID", "Duration", "Concentration_mg_per_L", "Value_Log10"]
).reset_index(drop=True)

# =========================
# Pair node label includes donor/recipient + plasmid + time + concentration
# =========================
def fmt_conc(x):
    return f"{x:g} mg/L"

significant_df["PairNode"] = (
    significant_df["Donor"] + " → " + significant_df["Recipient"]
    + "\n" + significant_df["Plasmid_ID"]
    + "\n" + significant_df["Duration"]
    + "\n" + significant_df["Concentration_mg_per_L"].apply(fmt_conc)
)

# =========================
# Network graph
# =========================
G = nx.DiGraph()
color_map = {"Induced": "#0077BB", "Inhibited": "#EE7733"}

for _, r in significant_df.iterrows():
    metal = r["Metal_form"]
    pair = r["PairNode"]
    effect = r["Effect"]

    G.add_node(metal, type="metal")
    G.add_node(pair, type="pair")
    G.add_edge(metal, pair, color=color_map[effect], weight=float(abs(r["Delta_Log10"])))

# =========================
# Visual style
# =========================
pale_yellow = "#FFFFB3"
pair_color = "#D3D3D3"

node_colors = [pair_color if G.nodes[n]["type"] == "pair" else pale_yellow for n in G.nodes()]
edge_colors = [G[u][v]["color"] for u, v in G.edges()]
edge_widths = [2.5 for _ in G.edges()]

# =========================
# Positioning
# =========================
metals = [n for n in G.nodes if G.nodes[n]["type"] == "metal"]
pairs = [n for n in G.nodes if G.nodes[n]["type"] == "pair"]

pos = {}

# Metals on inner circle
metal_angles = {m: a for m, a in zip(metals, np.linspace(0, 2*np.pi, max(len(metals), 1), endpoint=False))}
radius_metals = 3.5
for m in metals:
    a = metal_angles[m]
    pos[m] = (radius_metals * np.cos(a), radius_metals * np.sin(a))

# Pairs on outer circle
radius_pairs = 8.5  
for i, p in enumerate(pairs):
    a = i * (2*np.pi / max(len(pairs), 1))
    pos[p] = (radius_pairs * np.cos(a), radius_pairs * np.sin(a))

# =========================
# Draw
# =========================
plt.figure(figsize=(40, 28))
node_sizes = [40000 if G.nodes[n]["type"] == "pair" else 40000 for n in G.nodes()]

nx.draw(
    G, pos,
    with_labels=False,
    node_color=node_colors,
    edge_color=edge_colors,
    node_size=node_sizes,
    width=edge_widths,
    arrows=True,
    arrowstyle="-|>",
    arrowsize=25,
    connectionstyle="arc3,rad=0.15"
)

# Labels
for node, (x, y) in pos.items():
    plt.text(
        x, y, node,
        fontsize=20,  
        fontweight="bold",
        fontfamily="Arial",
        ha="center",
        va="center"
    )

# Legend
legend_elements = [
    Line2D([0], [0], color="#0077BB", lw=2, label="Induced Conjugation"),
    Line2D([0], [0], color="#EE7733", lw=2, label="Inhibited Conjugation"),
    Line2D([0], [0], marker="o", color="w", label="Metal", markerfacecolor=pale_yellow, markersize=20),
    Line2D([0], [0], marker="o", color="w",
           label="Donor→Recipient + Plasmid + Time + Conc.",
           markerfacecolor=pair_color, markersize=20),
]
plt.legend(handles=legend_elements, loc="upper left", fontsize=22)

plt.axis("off")
plt.savefig("", format="jpg", dpi=300, bbox_inches="tight")
plt.show()

# =========================
# Export all computed deltas (so you can check)
# =========================
significant_df.to_excel("", index=False)
