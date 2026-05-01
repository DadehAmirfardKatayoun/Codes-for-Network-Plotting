import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import re
import textwrap
from collections import defaultdict

# =========================
# USER SETTINGS
# =========================
FILE_PATH = r""
SHEET_NAME = ""

FIGSIZE = (42, 26)

# Fonts
FONT_FAMILY = "Arial"
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = [FONT_FAMILY, "Liberation Sans", "DejaVu Sans", "sans-serif"]
plt.rcParams["axes.unicode_minus"] = False

# Colors (colorblind-safe)
BLUE = "#0077BB"    # ↑
ORANGE = "#EE7733"  # ↓
METAL_COLOR = "#FFFFB3"
BACT_COLOR = "#D3D3D3"

# Sizes / fonts
METAL_NODE_SIZE = 26000
BACT_NODE_SIZE = 78000

METAL_FONT_SIZE = 22
BACT_NAME_FONT_SIZE = 19
BACT_LINE_FONT_SIZE = 15

EDGE_WIDTH = 2.8
ARROW_SIZE = 42

TITLE_FONT_SIZE = 20
LEGEND_FONT_SIZE = 18

# Curvature
RAD_SINGLE = 0.12
RAD_PAIR = (-0.22, 0.22)

# Wrapping inside bacteria nodes
LINE_WRAP_WIDTH = 34


# =========================
# Helpers
# =========================
def norm_text(x):
    if pd.isna(x):
        return ""
    s = str(x).strip()
    s = re.sub(r"\s+", " ", s)
    return s

def parse_numeric(x):
    """Extract first numeric value from strings like '< 13', '>30', '13–38', '26-131'."""
    if pd.isna(x):
        return float("nan")
    s = str(x).strip().replace("–", "-").replace("—", "-")
    m = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", s)
    return float(m.group()) if m else float("nan")

def format_range(min_v, max_v):
    if pd.isna(min_v) and pd.isna(max_v):
        return "NA"
    if pd.isna(min_v):
        return f"{max_v:g}"
    if pd.isna(max_v):
        return f"{min_v:g}"
    if float(min_v) == float(max_v):
        return f"{min_v:g}"
    return f"{min_v:g}–{max_v:g}"

def clean_conc_text(text, min_v, max_v):
    """Prefer original text (keeps <, >, ranges); otherwise fallback to numeric."""
    if pd.notna(text):
        s = str(text).strip()
        if s and s.lower() != "nan":
            # keep en-dash style
            s = s.replace("—", "–").replace("-", "–")
            return s
    return format_range(min_v, max_v)

def normalize_duration(d):
    """Standardize duration text a bit; keep the original words."""
    d = norm_text(d)
    if not d:
        return "Not_reported"
    # small clean-ups
    d = d.replace("Not reported", "Not_reported")
    d = d.replace("Not_reported", "Not_reported")
    return d

def unique_preserve_order(seq):
    """Works for strings or tuples; preserves first occurrence order."""
    seen = set()
    out = []
    for item in seq:
        if item in seen:
            continue
        seen.add(item)
        out.append(item)
    return out

def wrap_lines(s, width=34):
    return textwrap.wrap(s, width=width)

# Fix “cubic” ion rendering by mapping to roman-valence style
ION_DISPLAY_MAP = {
    "Cu²⁺": "Cu (II)", "Cu2+": "Cu (II)", "Cu²": "Cu (II)",
    "Zn²⁺": "Zn (II)", "Zn2+": "Zn (II)", "Zn²": "Zn (II)",
    "Ni²⁺": "Ni (II)", "Ni2+": "Ni (II)", "Ni²": "Ni (II)",
    "Co²⁺": "Co (II)", "Co2+": "Co (II)", "Co²": "Co (II)",
    "Mn²⁺": "Mn (II)", "Mn2+": "Mn (II)", "Mn²": "Mn (II)",
    "Fe³⁺": "Fe (III)", "Fe3+": "Fe (III)", "Fe³": "Fe (III)",
}

def metal_display_name(metal_raw: str) -> str:
    m = norm_text(metal_raw)
    return ION_DISPLAY_MAP.get(m, m)


# =========================
# Load Excel
# =========================
df = pd.read_excel(FILE_PATH, sheet_name=SHEET_NAME)

required = [
    "Heavy_metal(loid)",
    "Bacterial_identity",
    "Concentration_mg_per_L_text",
    "Concentration_min_mg_per_L",
    "Concentration_max_mg_per_L",
    "Effect_direction_symbol",
    "Duration",
]
missing = [c for c in required if c not in df.columns]
if missing:
    raise ValueError(f"Missing columns in Excel: {missing}")

df_clean = df.copy()

# Normalize
df_clean["Heavy_metal(loid)"] = df_clean["Heavy_metal(loid)"].apply(norm_text)
df_clean["Bacterial_identity"] = df_clean["Bacterial_identity"].apply(norm_text)
df_clean["Effect_direction_symbol"] = df_clean["Effect_direction_symbol"].apply(norm_text)
df_clean["Concentration_mg_per_L_text"] = df_clean["Concentration_mg_per_L_text"].apply(norm_text)
df_clean["Duration"] = df_clean["Duration"].apply(normalize_duration)

# Parse bounds
df_clean["min_mgL"] = df_clean["Concentration_min_mg_per_L"].apply(parse_numeric)
df_clean["max_mgL"] = df_clean["Concentration_max_mg_per_L"].apply(parse_numeric)

# If one bound missing, copy the other
df_clean.loc[df_clean["min_mgL"].isna(), "min_mgL"] = df_clean.loc[df_clean["min_mgL"].isna(), "max_mgL"]
df_clean.loc[df_clean["max_mgL"].isna(), "max_mgL"] = df_clean.loc[df_clean["max_mgL"].isna(), "min_mgL"]

# Keep only valid rows
df_clean = df_clean[
    (df_clean["Heavy_metal(loid)"] != "") &
    (df_clean["Bacterial_identity"] != "") &
    (df_clean["Effect_direction_symbol"].isin(["↑", "↓"]))
].copy()

# Conc text (keep <, >, ranges)
df_clean["Conc_text"] = df_clean.apply(
    lambda r: clean_conc_text(r["Concentration_mg_per_L_text"], r["min_mgL"], r["max_mgL"]),
    axis=1
)

# =========================
# Build bacteria inside-node text:
# per bacterium, per metal, split by effect (blue/orange)
# Store entries as (conc_text, duration)
# =========================
bact_effects = defaultdict(lambda: defaultdict(lambda: {"↑": [], "↓": []}))

for _, r in df_clean.iterrows():
    b = r["Bacterial_identity"]
    m = metal_display_name(r["Heavy_metal(loid)"])
    eff = r["Effect_direction_symbol"]
    conc = r["Conc_text"]
    dur = r["Duration"]
    bact_effects[b][m][eff].append((conc, dur))

# Deduplicate while preserving order
for b in bact_effects:
    for m in bact_effects[b]:
        for eff in ("↑", "↓"):
            bact_effects[b][m][eff] = unique_preserve_order(bact_effects[b][m][eff])

# =========================
# Build graph:
# - nodes: metals + bacteria
# - edges: for each metal->bacterium, create ↑ edge and/or ↓ edge
# =========================
G = nx.MultiDiGraph()

metals_raw = list(dict.fromkeys(df_clean["Heavy_metal(loid)"].tolist()))
metals = [metal_display_name(m) for m in metals_raw]
bacteria = list(dict.fromkeys(df_clean["Bacterial_identity"].tolist()))

for m in metals:
    G.add_node(m, type="metal")
for b in bacteria:
    G.add_node(b, type="bacterium")

# Collect which effects exist per pair
pair_effects = defaultdict(set)
for _, r in df_clean.iterrows():
    m = metal_display_name(r["Heavy_metal(loid)"])
    b = r["Bacterial_identity"]
    eff = r["Effect_direction_symbol"]
    pair_effects[(m, b)].add(eff)

for (m, b), effs in pair_effects.items():
    if "↑" in effs:
        G.add_edge(m, b, eff="↑", color=BLUE)
    if "↓" in effs:
        G.add_edge(m, b, eff="↓", color=ORANGE)

# =========================
# Layout
# =========================
metal_nodes = [n for n in G.nodes() if G.nodes[n].get("type") == "metal"]
bact_nodes = [n for n in G.nodes() if G.nodes[n].get("type") == "bacterium"]

pos = nx.shell_layout(G, nlist=[metal_nodes, bact_nodes])

# =========================
# Draw nodes
# =========================
plt.figure(figsize=FIGSIZE)

node_colors = [METAL_COLOR if G.nodes[n]["type"] == "metal" else BACT_COLOR for n in G.nodes()]
node_sizes = [METAL_NODE_SIZE if G.nodes[n]["type"] == "metal" else BACT_NODE_SIZE for n in G.nodes()]

nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes)

# =========================
# Draw edges with curvature:
# - one edge between pair -> RAD_SINGLE
# - two edges (↑ and ↓) -> RAD_PAIR
# =========================
pair_to_keys = defaultdict(list)
for u, v, k in G.edges(keys=True):
    pair_to_keys[(u, v)].append(k)

for (u, v), keys in pair_to_keys.items():
    if len(keys) == 1:
        k = keys[0]
        color = G.edges[u, v, k].get("color", "gray")
        nx.draw_networkx_edges(
            G, pos,
            edgelist=[(u, v)],
            edge_color=[color],
            width=EDGE_WIDTH,
            arrows=True,
            arrowstyle="-|>",
            arrowsize=ARROW_SIZE,
            connectionstyle=f"arc3,rad={RAD_SINGLE}",
            min_source_margin=12,
            min_target_margin=40
        )
    else:
        # two edges: assign consistent curvature by eff
        k_up = None
        k_dn = None
        for k in keys:
            eff = G.edges[u, v, k].get("eff")
            if eff == "↑":
                k_up = k
            elif eff == "↓":
                k_dn = k

        ordered = []
        if k_up is not None:
            ordered.append((k_up, RAD_PAIR[0]))
        if k_dn is not None:
            ordered.append((k_dn, RAD_PAIR[1]))

        for k, rad in ordered:
            color = G.edges[u, v, k].get("color", "gray")
            nx.draw_networkx_edges(
                G, pos,
                edgelist=[(u, v)],
                edge_color=[color],
                width=EDGE_WIDTH,
                arrows=True,
                arrowstyle="-|>",
                arrowsize=ARROW_SIZE,
                connectionstyle=f"arc3,rad={rad}",
                min_source_margin=12,
                min_target_margin=40
            )

# =========================
# Labels
# - Metals: label on nodes
# - Bacteria: name + colored (conc, duration) lines INSIDE nodes
# =========================
nx.draw_networkx_labels(
    G, pos,
    labels={m: m for m in metal_nodes},
    font_size=METAL_FONT_SIZE,
    font_weight="bold",
    font_family=FONT_FAMILY
)

for b in bact_nodes:
    x, y = pos[b]

    # Bacteria name
    plt.text(
        x, y + 0.070,
        b,
        ha="center", va="center",
        fontsize=BACT_NAME_FONT_SIZE,
        fontweight="bold",
        fontfamily=FONT_FAMILY,
        color="black"
    )

    line_y = y + 0.040
    line_step = 0.028

    # Print lines grouped by metal: first ↑ (blue), then ↓ (orange)
    for m in metal_nodes:
        if m not in bact_effects[b]:
            continue

        # Increase entries
        if bact_effects[b][m]["↑"]:
            parts = []
            for (conc, dur) in bact_effects[b][m]["↑"]:
                parts.append(f"{conc} ({dur})")
            conc_str = "; ".join(parts)
            text_line = f"{m}: {conc_str}"
            for w in wrap_lines(text_line, width=LINE_WRAP_WIDTH):
                plt.text(
                    x, line_y, w,
                    ha="center", va="center",
                    fontsize=BACT_LINE_FONT_SIZE,
                    fontfamily=FONT_FAMILY,
                    color=BLUE
                )
                line_y -= line_step

        # Decrease entries
        if bact_effects[b][m]["↓"]:
            parts = []
            for (conc, dur) in bact_effects[b][m]["↓"]:
                parts.append(f"{conc} ({dur})")
            conc_str = "; ".join(parts)
            text_line = f"{m}: {conc_str}"
            for w in wrap_lines(text_line, width=LINE_WRAP_WIDTH):
                plt.text(
                    x, line_y, w,
                    ha="center", va="center",
                    fontsize=BACT_LINE_FONT_SIZE,
                    fontfamily=FONT_FAMILY,
                    color=ORANGE
                )
                line_y -= line_step

# =========================
# Legend & Title
# =========================
legend_elements = [
    Line2D([0], [0], color=BLUE, lw=3, label="Enhances biofilm (↑)"),
    Line2D([0], [0], color=ORANGE, lw=3, label="Inhibits biofilm (↓)"),
    Line2D([0], [0], marker="o", color="w", label="Metal",
           markerfacecolor=METAL_COLOR, markersize=14),
    Line2D([0], [0], marker="o", color="w", label="Bacteria (colored conc + duration)",
           markerfacecolor=BACT_COLOR, markersize=14),
    Line2D([0], [0], color="black", lw=2, label="Arrow: Metal → Bacterium"),
]
plt.legend(handles=legend_elements, loc="upper left", fontsize=LEGEND_FONT_SIZE, frameon=True)

plt.title(
    "Network of Metal Effects on Bacterial Biofilm Formation\n(Concentrations + duration are shown only inside bacteria nodes; edges have no labels)",
    fontsize=TITLE_FONT_SIZE, fontfamily=FONT_FAMILY
)
plt.axis("off")
plt.subplots_adjust(left=0.03, right=0.97, top=0.92, bottom=0.03)

# =========================
# Save
# =========================
plt.savefig("biofilm_network_duration_inside_nodes.jpg", format="jpg", bbox_inches="tight")
plt.savefig(r"",
            format="jpg", dpi=600, bbox_inches="tight")
plt.show()
