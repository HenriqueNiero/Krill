#!/usr/bin/env python3
"""
krill_all_plots.py

Generate the consolidated set of plots/statistical summaries previously
produced from rede_node_table_completa_organizada.csv, using the current
Krill table DBs_BGCs_with_Hits_BiGSCAPE.tsv directly.

Primary input columns
---------------------
Database
product
product_bigscape
BiGSCAPE_Class
BiGSCAPE_Category
KnownResistanceHit_product
KnownResistanceHit_Resfam
regulatory_genes
genes
BGCs_Hits_Mean_Similarities(%)
Size
completeness

The script deliberately does NOT recreate obsolete intermediate CSV tables.
All transformations are performed in memory and the final plots/tables are
written directly to the output directory.

Usage
-----
python krill_all_plots.py DBs_BGCs_with_Hits_BiGSCAPE.tsv -o krill_plots

The script also reads Krill's DBs_normalized_info.tsv from
<input>/DBsReportOutput/ to generate database size/BGC-density plots.

Dependencies
------------
pandas, numpy, matplotlib, seaborn, scipy
"""

#!/usr/bin/env python3

import argparse, ast, itertools, math, re
from pathlib import Path
import ast
from upsetplot import UpSet, from_contents

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import pingouin as pg


from scipy import stats
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist


SIM_COL = "BGCs_Hits_Mean_Similarities(%)"


# ============================================================
# General functions
# ============================================================

def savefig(fig, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def create_dataset_colors(df):
    datasets = df["Database"].dropna().astype(str).str.strip().drop_duplicates().tolist()
    palette = sns.color_palette("husl", n_colors=len(datasets))
    return dict(zip(datasets, palette))


def dataset_colors_for(values, colors):
    return [colors.get(str(v).strip(), (0.5, 0.5, 0.5)) for v in values]


def split_cell(value):
    if pd.isna(value) or str(value).strip().lower() in {"", "nan", "none"}:
        return []
    return [x.strip() for x in str(value).split(",") if x.strip()]


# ============================================================
# Read main Krill table
# ============================================================

def read_input(path):
    df = pd.read_csv(path, sep="\t", low_memory=False)

    if "Database" not in df.columns:
        raise ValueError("The input table must contain a 'Database' column.")

    optional = [
        "product", "product_bigscape", "BiGSCAPE_Class",
        "BiGSCAPE_Category", "KnownResistanceHit_product",
        "KnownResistanceHit_Resfam", "regulatory_genes",
        "genes", SIM_COL, "Size", "completeness",
        "BGCs_Hits", "KnownResistenceHit"
    ]

    for c in optional:
        if c not in df.columns:
            df[c] = np.nan

    df["Database"] = df["Database"].astype(str).str.strip()

    for c in [SIM_COL, "Size"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    df["regulatory_genes"] = (
        df["regulatory_genes"].astype(str).str.lower().str.strip()
        .map({"true": True, "false": False, "1": True, "0": False})
        .fillna(False)
    )

    return df


# ============================================================
# Locate DBs_normalized_info.tsv
# ============================================================

def find_input_path(table_path, input_path=None):
    if input_path:
        return Path(input_path).resolve()

    table_path = Path(table_path).resolve()

    if table_path.parent.name.lower() == "dbsreportoutput":
        return table_path.parent.parent

    return table_path.parent


def read_normalized_info(table_path, input_path=None):
    input_path = find_input_path(table_path, input_path)
    path = input_path / "DBsReportOutput" / "DBs_normalized_info.tsv"

    if not path.exists():
        path = Path(table_path).resolve().parent / "DBs_normalized_info.tsv"

    if not path.exists():
        raise FileNotFoundError(
            f"Could not find DBs_normalized_info.tsv. Expected: {path}"
        )

    info = pd.read_csv(path, sep="\t", low_memory=False)

    db_col = next(
        (c for c in info.columns if str(c).strip().lower() in {"database", "db"}),
        None
    )

    nt_col = next(
        (
            c for c in info.columns
            if "nt" in str(c).lower().replace(" ", "").replace("_", "")
            and "kb" in str(c).lower().replace(" ", "").replace("_", "")
        ),
        None
    )

    if db_col is None or nt_col is None:
        raise ValueError(
            f"Could not identify Database or NT (KB) in {path}. "
            f"Columns: {list(info.columns)}"
        )

    info = info[[db_col, nt_col]].copy()
    info.columns = ["Database", "NT_KB"]
    info["Database"] = info["Database"].astype(str).str.strip()
    info["NT_KB"] = pd.to_numeric(info["NT_KB"], errors="coerce")
    info["NT_MB"] = info["NT_KB"] / 1000

    return (
        info.dropna(subset=["Database", "NT_MB"])
        .groupby("Database", as_index=False)[["NT_KB", "NT_MB"]]
        .sum()
    )


# ============================================================
# Database size / BGC / BGC density
# ============================================================

def plot_database_summary(table_path, input_path, outdir, df, colors):
    info = read_normalized_info(table_path, input_path)

    counts = (
        df.groupby("Database")
        .size()
        .reset_index(name="BGC_count")
    )

    summary = info.merge(counts, on="Database", how="left")
    summary["BGC_count"] = summary["BGC_count"].fillna(0).astype(int)
    summary["BGCs_per_MB"] = summary["BGC_count"] / summary["NT_MB"]

    summary.to_csv(
        outdir / "database_size_bgc_summary.csv",
        index=False
    )

    plots = [
        ("database_megabases.png", "NT (MB)", "Database size", "NT_MB"),
        ("database_bgc_count.png", "Number of BGCs", "BGC count", "BGC_count"),
        ("database_bgcs_per_megabase.png", "BGCs / MB", "BGC density", "BGCs_per_MB")
    ]

    n = summary["Database"].nunique()
    figsize = (10, max(2, n))

    for filename, xlabel, title, value_col in plots:
        fig, ax = plt.subplots(figsize=figsize)
        ax.barh(
            summary["Database"],
            summary[value_col],
            color=dataset_colors_for(summary["Database"], colors)
        )
        ax.set_ylabel("Database")
        ax.set_xlabel(xlabel)
        ax.set_title(title)
        ax.invert_yaxis()
        savefig(fig, outdir / filename)

    pd.DataFrame({
        "Database": list(colors.keys()),
        "Color": [matplotlib.colors.to_hex(c) for c in colors.values()]
    }).to_csv(outdir / "dataset_colors.csv", index=False)

    return summary



# ============================================================
# Database BGC completeness / similarity to MIBIG / Resistance hits
# ============================================================

def plot_cumulative_percentage_bars(df, outdir, colors):
    """
    Create 100%-stacked horizontal percentage bar plots for:
      1. Complete vs Fragmented BGCs
      2. BGCs with vs without MIBiG hits
      3. BGCs with vs without resistance hits
    """

    plots = [
        {
            "filename": "bgcs_completeness_percentage.png",
            "column": "completeness",
            "conditions": ["Complete", "Fragmented"],
            "labels": ["Complete", "Fragmented"],
            "title": "BGC completeness",
        },
        {
            "filename": "bgcs_mibig_similarity_percentage.png",
            "column": "BGCs_Hits",
            "conditions": ["With hits", "Without hits"],
            "title": "BGCs with similarity to MIBiG",
        },
        {
            "filename": "bgcs_resistance_percentage.png",
            "column": "KnownResistenceHit",
            "conditions": ["With hits", "Without hits"],
            "title": "BGCs with resistance hits",
        },
    ]

    datasets = (df["Database"].dropna().astype(str).str.strip().drop_duplicates().tolist())
    n = len(datasets)
    figsize = (10, max(3, 0.7 * n + 1))

    for plot in plots:
        column = plot["column"]

        if column not in df.columns:
            print(
                f"WARNING: Column '{column}' was not found. "
                f"Skipping {plot['filename']}."
            )
            continue

        data = df.copy()
        data[column] = data[column].fillna("").astype(str).str.strip()

        # ------------------------------------------------------------
        # Build the two categories for each plot
        # ------------------------------------------------------------

        if column == "completeness":
            data["Condition"] = data[column].str.lower().map(
                {
                    "complete": "Complete",
                    "fragmented": "Fragmented",
                }
            )

        elif column == "BGCs_Hits":
            data["Condition"] = np.where(
                data[column].str.strip() != "",
                "With hits",
                "Without hits",
            )

        elif column == "KnownResistenceHit":
            data["Condition"] = np.where(
                data[column].str.strip() != "",
                "With hits",
                "Without hits",
            )

        # Ignore unexpected completeness values
        data = data[data["Condition"].notna()].copy()

        # ------------------------------------------------------------
        # Count BGCs
        # ------------------------------------------------------------

        counts = (data.groupby(["Database", "Condition"]).size().unstack(fill_value=0))

        for condition in plot["conditions"]:
            if condition not in counts.columns:
                counts[condition] = 0

        counts = counts[plot["conditions"]]

        # Make sure every dataset is represented
        counts = counts.reindex(datasets, fill_value=0)

        # ------------------------------------------------------------
        # Convert counts to percentages
        # ------------------------------------------------------------

        totals = counts.sum(axis=1)

        percentages = counts.div(totals.replace(0, np.nan),axis=0) * 100
        percentages = percentages.fillna(0)

        # ------------------------------------------------------------
        # Plot
        # ------------------------------------------------------------

        fig, ax = plt.subplots(figsize=figsize)
        left = np.zeros(len(datasets))

        for i, condition in enumerate(plot["conditions"]):
            values = percentages[condition].values

            # --------------------------------------------------------
            # Same dataset color, but different opacity:
            # first condition = stronger
            # second condition = lighter
            # --------------------------------------------------------

            alpha = 0.95 if i == 0 else 0.50

            bar_colors = [
                colors.get(
                    dataset,
                    (0.5, 0.5, 0.5)
                )
                for dataset in datasets
            ]

            bars = ax.barh(
                datasets,
                values,
                left=left,
                color=bar_colors,
                alpha=alpha,
                edgecolor="white",
                linewidth=0.8,
                label=condition,
            )

            # --------------------------------------------------------
            # Number of BGCs inside each segment
            # --------------------------------------------------------

            for j, bar in enumerate(bars):
                count = int(counts.iloc[j][condition])
                percentage = percentages.iloc[j][condition]

                if count == 0:
                    continue

                # Only write the number if the segment is wide enough
                # to contain readable text.
                if percentage >= 3:

                    x_position = (left[j] + percentage / 2)

                    ax.text(
                        x_position,
                        bar.get_y() + bar.get_height() / 2,
                        str(count),
                        ha="center",
                        va="center",
                        color="black",
                        fontsize=9,
                        fontweight="bold",
                    )

            left += values

        # ------------------------------------------------------------
        # Axis formatting
        # ------------------------------------------------------------

        ax.set_xlim(0, 100)

        # ax.set_xlabel("Percentage of BGCs (%)")
        ax.set_ylabel("Database")
        ax.set_title(plot["title"])

        ax.set_xticks(np.arange(0, 101, 10))
        ax.set_xticklabels(
            [f"{x}%" for x in range(0, 101, 10)]
        )

        ax.invert_yaxis()

        # Legend describing the conditions
        legend_handles = [
            Patch(
                facecolor="gray",
                alpha=0.95,
                label=plot["conditions"][0]
            ),
            Patch(
                facecolor="gray",
                alpha=0.50,
                label=plot["conditions"][1]
            )
        ]

        ax.legend(
            handles=legend_handles,
            title="",
            loc="upper center",
            bbox_to_anchor=(0.5, -0.12),
            ncol=2,
            frameon=False,
        )

        savefig(fig, outdir / plot["filename"])


# ============================================================
# Database BGC size distribution
# ============================================================


def plot_size_violin_plots(df, outdir, colors, alpha=0.05):
    """
    Create three violin plots with boxplots inside:

    1. Size distribution by dataset.
    2. Size distribution by completeness (Complete vs Fragmented).
    3. Size distribution by completeness, separated by dataset.

    Also performs:
      - Kruskal-Wallis
      - Mann-Whitney U for Complete vs Fragmented
      - Games-Howell for dataset pairwise comparisons when
        Kruskal-Wallis is significant.

    Statistical results are written to TSV files.
    """

    from scipy.stats import kruskal, mannwhitneyu
    from statsmodels.stats.multitest import multipletests

    try:
        import pingouin as pg
        has_pingouin = True
    except ImportError:
        has_pingouin = False

    data = df.copy()

    # ------------------------------------------------------------
    # Prepare columns
    # ------------------------------------------------------------

    data["Database"] = (
        data["Database"]
        .astype(str)
        .str.strip()
    )

    data["Size"] = pd.to_numeric(
        data["Size"],
        errors="coerce"
    )

    data["completeness"] = (
        data["completeness"]
        .astype(str)
        .str.strip()
        .str.lower()
    )

    data["completeness"] = data["completeness"].map(
        {
            "complete": "Complete",
            "fragmented": "Fragmented"
        }
    )

    data = data.dropna(
        subset=["Database", "Size"]
    )

    # Remove zero/negative sizes
    data = data[data["Size"] > 0].copy()

    # ------------------------------------------------------------
    # Dataset order
    # ------------------------------------------------------------

    datasets = (
        data["Database"]
        .drop_duplicates()
        .tolist()
    )

    # ============================================================
    # PLOT 1
    # Size distribution by dataset
    # ============================================================

    fig, ax = plt.subplots(figsize=(max(7, 1.0 * len(datasets)), 7))

    sns.violinplot(
        data=data,
        x="Database",
        y="Size",
        hue="Database",
        palette=colors,
        inner=None,
        fill=False,
        cut=0,
        gap=0.1,
        linewidth=1,
        legend=False,
        ax=ax,
        zorder=2
    )

    sns.boxplot(
        data=data,
        x="Database",
        y="Size",
        hue="Database",
        fill=False,
        linewidth=1, 
        linecolor='black',
        gap=0.1,
        zorder=3,
        palette={
            d: "gray"
            for d in datasets
        },
        width=0.18,
        whis=(0, 100),
        legend=False,
        showcaps=True,
        boxprops={
            "zorder": 3,
        },
        whiskerprops={
            "color": "black"
        },
        capprops={
            "color": "black"
        },
        medianprops={
            "color": "black",
            "linewidth": 2
        },
        ax=ax,
    )

    sns.stripplot(
        data=df, 
        x="Database", 
        y="Size", 
        palette=colors, 
        hue="Database", 
        s=7, 
        linewidth=0.1, 
        edgecolor="black", 
        alpha=.3, 
        jitter=0.15, 
        ax=ax, 
        zorder=1
    )


    ax.set_xlabel("")
    ax.set_ylabel("BGC size (b)")
    ax.set_title("BGC size distribution by dataset")

    plt.xticks(rotation=0, ha="center")

    savefig(
        fig,
        outdir / "violin_size_by_dataset.png"
    )

    # ============================================================
    # PLOT 2
    # Size distribution: Complete vs Fragmented
    # ============================================================

    completeness_data = data.dropna(
        subset=["completeness"]
    ).copy()

    completeness_order = [
        "Complete",
        "Fragmented"
    ]

    fig, ax = plt.subplots(
        figsize=(7, 7)
    )

    sns.violinplot(
        data=completeness_data,
        x="completeness",
        y="Size",
        order=completeness_order,
        color="gray",
        inner=None,
        fill=False,
        cut=0,
        gap=0.1,
        linewidth=1,
        ax=ax,
        zorder=2
    )

    sns.boxplot(
        data=completeness_data,
        x="completeness",
        y="Size",
        fill=False,
        gap=0.1,
        linewidth=1,
        color="gray",
        zorder=3,
        order=completeness_order,
        width=0.18,
        whis=(0, 100),
        showcaps=True,
        boxprops={
            "color": "black",
            "zorder": 3,
        },
        whiskerprops={
            "color": "black"
        },
        capprops={
            "color": "black"
        },
        medianprops={
            "color": "black",
            "linewidth": 2
        },
        ax=ax,
    )

    sns.stripplot(
        data=completeness_data, 
        x="completeness", 
        y="Size",
        color="gray", 
        s=7, 
        linewidth=0.1, 
        edgecolor="black", 
        alpha=.3, 
        jitter=0.15, 
        ax=ax, 
        zorder=1
    )

    ax.set_xlabel("")
    ax.set_ylabel("BGC size (b)")
    ax.set_title("BGC size distribution by completeness")

    savefig(fig, outdir / "violin_size_by_completeness.png")

    # ------------------------------------------------------------
    # Mann-Whitney U: Complete vs Fragmented
    # ------------------------------------------------------------

    complete_values = completeness_data.loc[
        completeness_data["completeness"] == "Complete",
        "Size"
    ]

    fragmented_values = completeness_data.loc[
        completeness_data["completeness"] == "Fragmented",
        "Size"
    ]

    mann_whitney_results = []

    if len(complete_values) > 0 and len(fragmented_values) > 0:

        U, p = mannwhitneyu(
            complete_values,
            fragmented_values,
            alternative="two-sided"
        )

        mann_whitney_results.append(
            {
                "Comparison": "Complete vs Fragmented",
                "N_Complete": len(complete_values),
                "N_Fragmented": len(fragmented_values),
                "U": U,
                "p_value": p,
                "alpha": alpha,
                "Significant": p < alpha,
            }
        )

    pd.DataFrame(
        mann_whitney_results
    ).to_csv(
        outdir / "size_mann_whitney_complete_fragmented.tsv",
        sep="\t",
        index=False
    )

    # ============================================================
    # PLOT 3
    # Completeness + dataset
    # ============================================================

    fig, ax = plt.subplots(
        figsize=(max(9, 1.2 * len(datasets)), 7)
    )

    
    sns.violinplot(
        data=completeness_data,
        x="completeness",
        y="Size",
        hue="Database",
        order=completeness_order,
        hue_order=datasets,
        palette=colors,
        fill=False,
        cut=0,
        linewidth=1,
        ax=ax,
        zorder=2,
        inner=None,
        width=0.8,
        gap=0,
        dodge=True
    )

    sns.boxplot(
        data=completeness_data,
        x="completeness",
        y="Size",
        hue="Database",
        fill=False,
        order=completeness_order,
        hue_order=datasets,
        color="gray",
        palette={
            d: "gray"
            for d in datasets
        },
        width=0.8,
        gap=0.7,
        linewidth=1,
        whis=(0, 100),
        dodge=True,
        legend=False,
        showcaps=True,
        boxprops={
            "zorder": 3,
        },
        whiskerprops={
            "color": "black"
        },
        capprops={
            "color": "black"
        },
        medianprops={
            "color": "black",
            "linewidth": 1.5
        },
        ax=ax,
        zorder=3
    )
    
    sns.stripplot(
        data=completeness_data,
        x="completeness",
        y="Size",
        hue="Database",
        order=completeness_order,
        hue_order=datasets,
        palette=colors, 
        s=7,
        linewidth=0.1, 
        edgecolor="black", 
        alpha=0.3, 
        jitter=0.15,
        legend=False,
        ax=ax, 
        zorder=1,
        dodge=True,
    )

    ax.set_xlabel("")
    ax.set_ylabel("BGC size (bp)")
    ax.set_title(
        "BGC size distribution by completeness and dataset"
    )

    ax.legend(
        title="Dataset",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        frameon=False,
    )

    savefig(fig, outdir / "violin_size_completeness_by_dataset.png")

    # ============================================================
    # KRUSKAL-WALLIS
    # ============================================================

    dataset_groups = []

    for dataset in datasets:

        values = data.loc[
            data["Database"] == dataset,
            "Size"
        ].dropna()

        if len(values) > 0:
            dataset_groups.append(values)

    kruskal_results = []

    if len(dataset_groups) >= 2:

        H, p = kruskal(
            *dataset_groups
        )

        kruskal_results.append(
            {
                "Comparison": "All datasets",
                "K": len(dataset_groups),
                "H": H,
                "p_value": p,
                "alpha": alpha,
                "Significant": p < alpha,
            }
        )

    pd.DataFrame(
        kruskal_results
    ).to_csv(
        outdir / "size_kruskal_wallis_datasets.tsv",
        sep="\t",
        index=False
    )

    # ============================================================
    # GAMES-HOWELL
    # Only if Kruskal-Wallis is significant
    # ============================================================

    games_howell_results = []

    if (
        len(kruskal_results) > 0
        and kruskal_results[0]["p_value"] < alpha
    ):

        if not has_pingouin:

            print(
                "\nWARNING: pingouin is not installed."
            )

            print(
                "Games-Howell was not calculated."
            )

            print(
                "Install it with:"
            )

            print(
                "pip install pingouin"
            )

        else:

            gh_data = data[
                ["Database", "Size"]
            ].dropna().copy()

            gh = pg.pairwise_gameshowell(
                data=gh_data,
                dv="Size",
                between="Database"
            )

            gh = gh.rename(
                columns={
                    "A": "Dataset_1",
                    "B": "Dataset_2",
                    "mean(A)": "Mean_1",
                    "mean(B)": "Mean_2",
                    "diff": "Mean_Difference",
                    "se": "SE",
                    "T": "T",
                    "df": "df",
                    "pval": "p_value",
                    "hedges": "Hedges_g",
                }
            )

            gh["alpha"] = alpha
            gh["Significant"] = (
                gh["p_value"] < alpha
            )

            gh.to_csv(
                outdir / "size_games_howell_datasets.tsv",
                sep="\t",
                index=False
            )

            games_howell_results = gh.to_dict(
                orient="records"
            )

    # ============================================================
    # Summary of statistics
    # ============================================================

    stats_summary = []

    for result in kruskal_results:
        stats_summary.append(
            {
                "Test": "Kruskal-Wallis",
                "Comparison": result["Comparison"],
                "Statistic": result["H"],
                "p_value": result["p_value"],
                "Significant": result["Significant"],
            }
        )

    for result in mann_whitney_results:
        stats_summary.append(
            {
                "Test": "Mann-Whitney U",
                "Comparison": result["Comparison"],
                "Statistic": result["U"],
                "p_value": result["p_value"],
                "Significant": result["Significant"],
            }
        )

    for result in games_howell_results:

        stats_summary.append(
            {
                "Test": "Games-Howell",
                "Comparison": (
                    f'{result["Dataset_1"]} vs '
                    f'{result["Dataset_2"]}'
                ),
                "Statistic": result.get("T", np.nan),
                "p_value": result["p_value"],
                "Significant": result["Significant"],
            }
        )

    pd.DataFrame(
        stats_summary
    ).to_csv(
        outdir / "size_statistics_summary.tsv",
        sep="\t",
        index=False
    )


# ============================================================
# Database BGC similarity distribution
# ============================================================


def plot_similarity_violin_plots(df, outdir, colors, alpha=0.05):
    """
    Create three violin plots with boxplots inside:

    1. similarity distribution by dataset.
    2. similarity distribution by completeness (Complete vs Fragmented).
    3. similarity distribution by completeness, separated by dataset.

    Also performs:
      - Kruskal-Wallis
      - Mann-Whitney U for Complete vs Fragmented
      - Games-Howell for dataset pairwise comparisons when
        Kruskal-Wallis is significant.

    Statistical results are written to TSV files.
    """

    from scipy.stats import kruskal, mannwhitneyu
    from statsmodels.stats.multitest import multipletests

    try:
        import pingouin as pg
        has_pingouin = True
    except ImportError:
        has_pingouin = False

    data = df.copy()

    # ------------------------------------------------------------
    # Prepare columns
    # ------------------------------------------------------------

    data["Database"] = (
        data["Database"]
        .astype(str)
        .str.strip()
    )

    data["BGCs_Hits_BestSimilarity"] = pd.to_numeric(
        data["BGCs_Hits_BestSimilarity"],
        errors="coerce"
    )

    data["completeness"] = (
        data["completeness"]
        .astype(str)
        .str.strip()
        .str.lower()
    )

    data["completeness"] = data["completeness"].map(
        {
            "complete": "Complete",
            "fragmented": "Fragmented"
        }
    )

    data = data.dropna(
        subset=["Database", "BGCs_Hits_BestSimilarity"]
    )

    # Remove zero/negative similarities
    data = data[data["BGCs_Hits_BestSimilarity"] > 0].copy()

    # ------------------------------------------------------------
    # Dataset order
    # ------------------------------------------------------------

    datasets = (
        data["Database"]
        .drop_duplicates()
        .tolist()
    )

    # ============================================================
    # PLOT 1
    # Similarity distribution by dataset
    # ============================================================

    fig, ax = plt.subplots(figsize=(max(7, 1.0 * len(datasets)), 7))

    sns.violinplot(
        data=data,
        x="Database",
        y="BGCs_Hits_BestSimilarity",
        hue="Database",
        palette=colors,
        inner=None,
        fill=False,
        cut=0,
        gap=0.1,
        linewidth=1,
        legend=False,
        ax=ax,
        zorder=2
    )

    sns.boxplot(
        data=data,
        x="Database",
        y="BGCs_Hits_BestSimilarity",
        hue="Database",
        fill=False,
        linewidth=1, 
        linecolor='black',
        gap=0.1,
        zorder=3,
        palette={
            d: "gray"
            for d in datasets
        },
        width=0.18,
        whis=(0, 100),
        legend=False,
        showcaps=True,
        boxprops={
            "zorder": 3,
        },
        whiskerprops={
            "color": "black"
        },
        capprops={
            "color": "black"
        },
        medianprops={
            "color": "black",
            "linewidth": 2
        },
        ax=ax,
    )

    sns.stripplot(
        data=df, 
        x="Database", 
        y="BGCs_Hits_BestSimilarity", 
        palette=colors, 
        hue="Database", 
        s=7, 
        linewidth=0.1, 
        edgecolor="black", 
        alpha=.3, 
        jitter=0.15, 
        ax=ax, 
        zorder=1
    )


    ax.set_xlabel("")
    ax.set_ylabel("BGC similarity (%)")
    ax.set_title("BGC MIBIG similarity distribution by dataset")

    plt.xticks(rotation=0, ha="center")

    savefig(fig, outdir / "violin_similarity_by_dataset.png")

    # ============================================================
    # PLOT 2
    # Similarity distribution: Complete vs Fragmented
    # ============================================================

    completeness_data = data.dropna(
        subset=["completeness"]
    ).copy()

    completeness_order = [
        "Complete",
        "Fragmented"
    ]

    fig, ax = plt.subplots(
        figsize=(7, 7)
    )

    sns.violinplot(
        data=completeness_data,
        x="completeness",
        y="BGCs_Hits_BestSimilarity",
        order=completeness_order,
        color="gray",
        inner=None,
        fill=False,
        cut=0,
        gap=0.1,
        linewidth=1,
        ax=ax,
        zorder=2
    )

    sns.boxplot(
        data=completeness_data,
        x="completeness",
        y="BGCs_Hits_BestSimilarity",
        fill=False,
        gap=0.1,
        linewidth=1,
        color="gray",
        zorder=3,
        order=completeness_order,
        width=0.18,
        whis=(0, 100),
        showcaps=True,
        boxprops={
            "color": "black",
            "zorder": 3,
        },
        whiskerprops={
            "color": "black"
        },
        capprops={
            "color": "black"
        },
        medianprops={
            "color": "black",
            "linewidth": 2
        },
        ax=ax,
    )

    sns.stripplot(
        data=completeness_data, 
        x="completeness", 
        y="BGCs_Hits_BestSimilarity",
        color="gray", 
        s=7, 
        linewidth=0.1, 
        edgecolor="black", 
        alpha=.3, 
        jitter=0.15, 
        ax=ax, 
        zorder=1
    )

    ax.set_xlabel("")
    ax.set_ylabel("BGC similarity (%)")
    ax.set_title("BGC MIBIG similarity distribution by completeness")

    savefig(fig, outdir / "violin_similarity_by_completeness.png")

    # ------------------------------------------------------------
    # Mann-Whitney U: Complete vs Fragmented
    # ------------------------------------------------------------

    complete_values = completeness_data.loc[
        completeness_data["completeness"] == "Complete",
        "BGCs_Hits_BestSimilarity"
    ]

    fragmented_values = completeness_data.loc[
        completeness_data["completeness"] == "Fragmented",
        "BGCs_Hits_BestSimilarity"
    ]

    mann_whitney_results = []

    if len(complete_values) > 0 and len(fragmented_values) > 0:

        U, p = mannwhitneyu(
            complete_values,
            fragmented_values,
            alternative="two-sided"
        )

        mann_whitney_results.append(
            {
                "Comparison": "Complete vs Fragmented",
                "N_Complete": len(complete_values),
                "N_Fragmented": len(fragmented_values),
                "U": U,
                "p_value": p,
                "alpha": alpha,
                "Significant": p < alpha,
            }
        )

    pd.DataFrame(
        mann_whitney_results
    ).to_csv(
        outdir / "similarity_mann_whitney_complete_fragmented.tsv",
        sep="\t",
        index=False
    )

    # ============================================================
    # PLOT 3
    # Completeness + dataset
    # ============================================================

    fig, ax = plt.subplots(
        figsize=(max(9, 1.2 * len(datasets)), 7)
    )

    
    sns.violinplot(
        data=completeness_data,
        x="completeness",
        y="BGCs_Hits_BestSimilarity",
        hue="Database",
        order=completeness_order,
        hue_order=datasets,
        palette=colors,
        fill=False,
        cut=0,
        linewidth=1,
        ax=ax,
        zorder=2,
        inner=None,
        width=0.8,
        gap=0,
        dodge=True
    )

    sns.boxplot(
        data=completeness_data,
        x="completeness",
        y="BGCs_Hits_BestSimilarity",
        hue="Database",
        fill=False,
        order=completeness_order,
        hue_order=datasets,
        color="gray",
        palette={
            d: "gray"
            for d in datasets
        },
        width=0.8,
        gap=0.7,
        linewidth=1,
        whis=(0, 100),
        dodge=True,
        legend=False,
        showcaps=True,
        boxprops={
            "zorder": 3,
        },
        whiskerprops={
            "color": "black"
        },
        capprops={
            "color": "black"
        },
        medianprops={
            "color": "black",
            "linewidth": 1.5
        },
        ax=ax,
        zorder=3
    )
    
    sns.stripplot(
        data=completeness_data,
        x="completeness",
        y="BGCs_Hits_BestSimilarity",
        hue="Database",
        order=completeness_order,
        hue_order=datasets,
        palette=colors, 
        s=7,
        linewidth=0.1, 
        edgecolor="black", 
        alpha=0.3, 
        jitter=0.15,
        legend=False,
        ax=ax, 
        zorder=1,
        dodge=True,
    )

    ax.set_xlabel("")
    ax.set_ylabel("BGC similarity (%)")
    ax.set_title(
        "BGC MIBIG similarity distribution by completeness and dataset"
    )

    ax.legend(
        title="Dataset",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        frameon=False,
    )

    savefig(fig, outdir / "violin_similarity_completeness_by_dataset.png")

    # ============================================================
    # KRUSKAL-WALLIS
    # ============================================================

    dataset_groups = []

    for dataset in datasets:

        values = data.loc[
            data["Database"] == dataset,
            "BGCs_Hits_BestSimilarity"
        ].dropna()

        if len(values) > 0:
            dataset_groups.append(values)

    kruskal_results = []

    if len(dataset_groups) >= 2:

        H, p = kruskal(
            *dataset_groups
        )

        kruskal_results.append(
            {
                "Comparison": "All datasets",
                "K": len(dataset_groups),
                "H": H,
                "p_value": p,
                "alpha": alpha,
                "Significant": p < alpha,
            }
        )

    pd.DataFrame(
        kruskal_results
    ).to_csv(
        outdir / "similarity_kruskal_wallis_datasets.tsv",
        sep="\t",
        index=False
    )

    # ============================================================
    # GAMES-HOWELL
    # Only if Kruskal-Wallis is significant
    # ============================================================

    games_howell_results = []

    if (
        len(kruskal_results) > 0
        and kruskal_results[0]["p_value"] < alpha
    ):

        if not has_pingouin:

            print(
                "\nWARNING: pingouin is not installed."
            )

            print(
                "Games-Howell was not calculated."
            )

            print(
                "Install it with:"
            )

            print(
                "pip install pingouin"
            )

        else:

            gh_data = data[
                ["Database", "BGCs_Hits_BestSimilarity"]
            ].dropna().copy()

            gh = pg.pairwise_gameshowell(
                data=gh_data,
                dv="BGCs_Hits_BestSimilarity",
                between="Database"
            )

            gh = gh.rename(
                columns={
                    "A": "Dataset_1",
                    "B": "Dataset_2",
                    "mean(A)": "Mean_1",
                    "mean(B)": "Mean_2",
                    "diff": "Mean_Difference",
                    "se": "SE",
                    "T": "T",
                    "df": "df",
                    "pval": "p_value",
                    "hedges": "Hedges_g",
                }
            )

            gh["alpha"] = alpha
            gh["Significant"] = (
                gh["p_value"] < alpha
            )

            gh.to_csv(
                outdir / "similarity_games_howell_datasets.tsv",
                sep="\t",
                index=False
            )

            games_howell_results = gh.to_dict(
                orient="records"
            )

    # ============================================================
    # Summary of statistics
    # ============================================================

    stats_summary = []

    for result in kruskal_results:
        stats_summary.append(
            {
                "Test": "Kruskal-Wallis",
                "Comparison": result["Comparison"],
                "Statistic": result["H"],
                "p_value": result["p_value"],
                "Significant": result["Significant"],
            }
        )

    for result in mann_whitney_results:
        stats_summary.append(
            {
                "Test": "Mann-Whitney U",
                "Comparison": result["Comparison"],
                "Statistic": result["U"],
                "p_value": result["p_value"],
                "Significant": result["Significant"],
            }
        )

    for result in games_howell_results:

        stats_summary.append(
            {
                "Test": "Games-Howell",
                "Comparison": (
                    f'{result["Dataset_1"]} vs '
                    f'{result["Dataset_2"]}'
                ),
                "Statistic": result.get("T", np.nan),
                "p_value": result["p_value"],
                "Significant": result["Significant"],
            }
        )

    pd.DataFrame(
        stats_summary
    ).to_csv(
        outdir / "similarity_statistics_summary.tsv",
        sep="\t",
        index=False
    )


# ============================================================
# Database BGC product distribution
# ============================================================


def plot_product_and_bigscape_summaries(df, outdir, input_path=None):
    """
    Generate:
      1. All antiSMASH product occurrences.
      2. Top 20 antiSMASH products.
      3. Pie chart of antiSMASH product occurrences.
      4. All BiG-SCAPE category occurrences.
      5. Pie chart of BiG-SCAPE category occurrences.
      6. 100% stacked BiG-SCAPE category plot by dataset,
         with genome size (MB) shown as dots on a secondary axis.

    Comma-separated product values are split and counted separately.
    """

    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------
    # Validate and prepare data
    # ------------------------------------------------------------

    required = ["Database", "product", "BiGSCAPE_Category"]

    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(
            f"Missing required columns: {missing}"
        )

    data = df.copy()

    data["Database"] = (
        data["Database"]
        .fillna("")
        .astype(str)
        .str.strip()
    )

    datasets = data.loc[
        data["Database"].ne(""),
        "Database"
    ].drop_duplicates().tolist()

    # ------------------------------------------------------------
    # Split comma-separated values
    # ------------------------------------------------------------

    def split_values(value):
        if pd.isna(value):
            return []

        value = str(value).strip()

        if value.lower() in {"", "nan", "none", "na"}:
            return []

        return [
            item.strip()
            for item in value.split(",")
            if item.strip()
        ]

    # ------------------------------------------------------------
    # antiSMASH product counts
    # ------------------------------------------------------------

    product_records = []

    for _, row in data.iterrows():
        for product in split_values(row["product"]):
            product_records.append({
                "Database": row["Database"],
                "Product": product
            })

    product_long = pd.DataFrame(
        product_records,
        columns=["Database", "Product"]
    )

    if product_long.empty:
        product_counts = pd.DataFrame(
            columns=["Product", "Occurrences"]
        )
    else:
        product_counts = (
            product_long["Product"]
            .value_counts()
            .rename_axis("Product")
            .reset_index(name="Occurrences")
        )

    product_counts.to_csv(
        outdir / "antismash_product_occurrences.tsv",
        sep="\t",
        index=False
    )

    # Dynamic, consistent colors for each product type.
    product_names = product_counts["Product"].tolist()

    product_palette = dict(zip(
        product_names,
        sns.color_palette(
            "husl",
            n_colors=max(1, len(product_names))
        )
    ))

    # ------------------------------------------------------------
    # BiG-SCAPE category colors
    # ------------------------------------------------------------

    def category_color(category):
        cat = str(category).strip().lower().replace(" ", "")

        # Mixed NRPS/PKS categories take precedence.
        if "nrps" in cat and "pks" in cat:
            return "#F1C40F"  # Yellow

        if "ripp" in cat:
            return "#E74C3C"  # Red

        if "nrps" in cat:
            return "#2ECC71"  # Green

        if "pks" in cat:
            return "#E67E22"  # Orange

        if "terpene" in cat:
            return "#8E44AD"  # Purple

        return "#808080"      # Gray / other

    # ------------------------------------------------------------
    # Shared horizontal bar plotting function
    # ------------------------------------------------------------

    def horizontal_count_plot(
        counts,
        label_col,
        filename,
        title,
        color_func
    ):
        if counts.empty:
            print(f"Skipping empty plot: {filename}")
            return

        counts = counts.sort_values(
            "Occurrences",
            ascending=True
        ).copy()

        fig, ax = plt.subplots(
            figsize=(11, max(4, 0.35 * len(counts) + 1))
        )

        bars = ax.barh(
            counts[label_col],
            counts["Occurrences"],
            color=[
                color_func(value)
                for value in counts[label_col]
            ],
            edgecolor="black",
            linewidth=0.4
        )

        ax.bar_label(
            bars,
            labels=[
                str(int(value))
                for value in counts["Occurrences"]
            ],
            padding=3,
            fontsize=8,
            color="black"
        )

        ax.set_xlabel("Number of occurrences")
        ax.set_ylabel(label_col)
        ax.set_title(title)

        savefig(fig, outdir / filename)

    # ------------------------------------------------------------
    # All products and top 20 products
    # ------------------------------------------------------------

    horizontal_count_plot(
        product_counts,
        "Product",
        "antismash_all_product_occurrences.png",
        "antiSMASH product occurrences",
        lambda value: product_palette[value]
    )

    top20_products = (
        product_counts
        .sort_values("Occurrences", ascending=False)
        .head(20)
    )

    top20_products.to_csv(
        outdir / "antismash_top20_product_occurrences.tsv",
        sep="\t",
        index=False
    )

    horizontal_count_plot(
        top20_products,
        "Product",
        "antismash_top20_product_occurrences.png",
        "Top 20 antiSMASH product occurrences",
        lambda value: product_palette[value]
    )

    # ------------------------------------------------------------
    # Product pie chart
    # ------------------------------------------------------------

    if not product_counts.empty:
        pie_data = product_counts.sort_values(
            "Occurrences",
            ascending=False
        )

        fig, ax = plt.subplots(figsize=(11, 9))

        ax.pie(
            pie_data["Occurrences"],
            labels=pie_data["Product"],
            colors=[
                product_palette[value]
                for value in pie_data["Product"]
            ],
            autopct=lambda p: f"{p:.1f}%" if p >= 2 else "",
            startangle=90,
            pctdistance=0.72,
            textprops={"fontsize": 8}
        )

        ax.set_title(
            "Distribution of antiSMASH product occurrences"
        )

        savefig(
            fig,
            outdir / "antismash_product_occurrences_pie.png"
        )

    # ------------------------------------------------------------
    # BiG-SCAPE category counts
    # Each original category value is counted as one category.
    # For example, NRPS.PKS remains a single category.
    # ------------------------------------------------------------

    category_data = data[
        ["Database", "BiGSCAPE_Category"]
    ].copy()

    category_data["BiGSCAPE_Category"] = (
        category_data["BiGSCAPE_Category"]
        .fillna("")
        .astype(str)
        .str.strip()
    )

    category_data = category_data[
        category_data["BiGSCAPE_Category"].ne("")
        & ~category_data["BiGSCAPE_Category"].str.lower().isin(
            {"nan", "none", "na"}
        )
    ]

    category_counts = (
        category_data["BiGSCAPE_Category"]
        .value_counts()
        .rename_axis("BiGSCAPE_Category")
        .reset_index(name="Occurrences")
    )

    category_counts.to_csv(
        outdir / "bigscape_category_occurrences.tsv",
        sep="\t",
        index=False
    )

    horizontal_count_plot(
        category_counts,
        "BiGSCAPE_Category",
        "bigscape_all_category_occurrences.png",
        "BiG-SCAPE category occurrences",
        category_color
    )

    # ------------------------------------------------------------
    # BiG-SCAPE category pie chart
    # ------------------------------------------------------------

    if not category_counts.empty:
        pie_data = category_counts.sort_values(
            "Occurrences",
            ascending=False
        )

        fig, ax = plt.subplots(figsize=(10, 8))

        ax.pie(
            pie_data["Occurrences"],
            labels=pie_data["BiGSCAPE_Category"],
            colors=[
                category_color(value)
                for value in pie_data["BiGSCAPE_Category"]
            ],
            autopct=lambda p: f"{p:.1f}%" if p >= 2 else "",
            startangle=90,
            pctdistance=0.72,
            textprops={"fontsize": 9}
        )

        ax.set_title(
            "Distribution of BiG-SCAPE category occurrences"
        )

        savefig(
            fig,
            outdir / "bigscape_category_occurrences_pie.png"
        )

    # ------------------------------------------------------------
    # 100% stacked category percentages by dataset
    # ------------------------------------------------------------

    category_by_dataset = (category_data.groupby(["Database", "BiGSCAPE_Category"]).size().unstack(fill_value=0))

    category_by_dataset = category_by_dataset.reindex(datasets, fill_value=0)

    categories = sorted(category_by_dataset.columns.tolist())

    category_by_dataset = category_by_dataset[categories]

    totals = category_by_dataset.sum(axis=1)

    percentages = (category_by_dataset.div(totals.replace(0, np.nan), axis=0).mul(100).fillna(0))

    percentages.to_csv(outdir / "bigscape_category_percentages_by_dataset.tsv", sep="\t", index_label="Database")


    # ------------------------------------------------------------
    # Load dataset sizes from the existing KrillPlots summary
    # ------------------------------------------------------------

    size_path = outdir / "database_size_bgc_summary.csv"

    if not size_path.exists():
        raise FileNotFoundError(
            f"Could not find dataset size summary: {size_path}"
        )

    info = pd.read_csv(size_path, low_memory=False)

    required_size_columns = {"Database", "NT_MB"}
    missing = required_size_columns - set(info.columns)

    if missing:
        raise ValueError(
            f"{size_path} is missing columns: {sorted(missing)}"
        )

    info["Database"] = (
        info["Database"].astype(str).str.strip()
    )
    info["NT_MB"] = pd.to_numeric(
        info["NT_MB"], errors="coerce"
    )

    info = (
        info.dropna(subset=["Database", "NT_MB"])
        .groupby("Database", as_index=False)["NT_MB"]
        .sum()
    )

    genome_sizes = (
        info.set_index("Database")["NT_MB"]
        .reindex(datasets)
    )

    # Save summary values for checking.
    pd.DataFrame({
        "Database": datasets,
        "NT_MB": genome_sizes.to_numpy(),
        "BGCs_with_category": totals.reindex(datasets).to_numpy()
    }).to_csv(
        outdir / "bigscape_category_dataset_summary.tsv",
        sep="\t",
        index=False
    )

    # ------------------------------------------------------------
    # Stacked percentage plot + genome-size dots
    # ------------------------------------------------------------

    x = np.arange(len(datasets))

    fig, ax1 = plt.subplots(
        figsize=(max(10, 1.2 * len(datasets)), 7)
    )

    bottom = np.zeros(len(datasets))

    for category in categories:
        values = percentages[category].reindex(
            datasets,
            fill_value=0
        ).to_numpy()

        ax1.bar(
            x,
            values,
            bottom=bottom,
            color=category_color(category),
            edgecolor="white",
            linewidth=0.5,
            label=category
        )

        bottom += values
    
    ax1.grid(False)
    ax1.set_ylim(0, 100)
    ax1.set_ylabel("BGCs by BiG-SCAPE category (%)")
    ax1.set_xlabel("Dataset")
    ax1.set_xticks(x)
    ax1.set_xticklabels(
        datasets,
        rotation=45,
        ha="right"
    )

    ax1.set_title(
        "BiG-SCAPE category composition and dataset size"
    )

    ax2 = ax1.twinx()

    size_values = genome_sizes.to_numpy()
    valid = np.isfinite(size_values)

    ax2.scatter(
        x[valid],
        size_values[valid],
        color="black",
        marker="o",
        s=65,
        edgecolor="white",
        linewidth=0.7,
        zorder=10
    )

    ax2.set_ylabel("Dataset size (MB)")

    ax2.grid(False)

    category_handles = [
        Patch(
            facecolor=category_color(category),
            edgecolor="white",
            label=category
        )
        for category in categories
    ]

    size_handle = Line2D(
        [0], [0],
        marker="o",
        color="black",
        linestyle="None",
        markersize=8,
        label="Dataset size (MB)"
    )

    ax1.legend(
        handles=category_handles + [size_handle],
        title="BiG-SCAPE category / dataset size",
        bbox_to_anchor=(1.20, 1),
        loc="upper left",
        frameon=False
    )

    savefig(
        fig,
        outdir / "bigscape_category_percentage_and_dataset_size.png"
    )

    print("Product and BiG-SCAPE summary plots generated.")



# ============================================================
# Datasets antiSMASH product heatmap
# ============================================================

def plot_product_heatmap(df, outdir):
    """
    Heatmap of antiSMASH product annotations by dataset.

    Y-axis: original product cell values, kept intact.
    X-axis: datasets.
    Cell values: number of BGCs with each exact product annotation.
    Comma-separated values are NOT split.
    """

    data = df.copy()

    data["Database"] = (
        data["Database"].fillna("").astype(str).str.strip()
    )

    data["product"] = (
        data["product"].fillna("").astype(str).str.strip()
    )

    # Exclude empty or missing product annotations.
    data = data[
        data["Database"].ne("")
        & data["product"].ne("")
        & ~data["product"].str.lower().isin(
            {"nan", "none", "na"}
        )
    ].copy()

    if data.empty:
        print("No antiSMASH product annotations found; skipping heatmap.")
        return

    # Keep each complete cell value as one category.
    counts = pd.crosstab(
        data["product"],
        data["Database"]
    )

    # Preserve dataset order used in the input table.
    datasets = (
        df["Database"]
        .dropna()
        .astype(str)
        .str.strip()
        .loc[lambda s: s.ne("")]
        .drop_duplicates()
        .tolist()
    )

    counts = counts.reindex(
        columns=datasets,
        fill_value=0
    )

    # Sort product annotations by total BGC count, descending.
    counts = counts.loc[
        counts.sum(axis=1).sort_values(ascending=False).index
    ]

    # Save the count matrix.
    counts.to_csv(
        outdir / "antismash_product_heatmap_counts.tsv",
        sep="\t",
        index_label="Product"
    )

    # Dynamic figure size based on datasets and annotations.
    fig, ax = plt.subplots(
        figsize=(
            max(8, 1.2 * len(datasets)),
            max(6, 0.32 * len(counts))
        )
    )

    sns.heatmap(
        counts,
        annot=True,
        fmt="d",
        cmap="YlGnBu",
        linewidths=0.3,
        linecolor="white",
        cbar_kws={"label": "Number of BGCs"},
        ax=ax
    )

    ax.set_xlabel("Dataset")
    ax.set_ylabel("antiSMASH product annotation")
    ax.set_title("antiSMASH product annotations by dataset")

    ax.tick_params(axis="x", labelrotation=45)
    ax.tick_params(axis="y", labelrotation=0)

    savefig(fig, outdir / "antismash_product_heatmap.png")


# ============================================================
# antiSMASH product heatmap - percentages by dataset
# ============================================================

def plot_product_heatmap_percentages(df, outdir):
    """
    Heatmap of antiSMASH product annotations by dataset.

    Y-axis: original product cell values, kept intact.
    X-axis: datasets.
    Cell values: percentage of BGCs with each exact product
                 annotation within each dataset.

    Percentages are calculated independently for each dataset.

    Example:
        LlaimaVolcano has 20 product occurrences.
        RiPP-like occurs 2 times.
        Percentage = (2 / 20) * 100 = 10%.

    Comma-separated product values are NOT split.
    """

    data = df.copy()

    data["Database"] = (
        data["Database"]
        .fillna("")
        .astype(str)
        .str.strip()
    )

    data["product"] = (
        data["product"]
        .fillna("")
        .astype(str)
        .str.strip()
    )

    # Exclude empty or missing product annotations.
    data = data[
        data["Database"].ne("")
        & data["product"].ne("")
        & ~data["product"].str.lower().isin(
            {"nan", "none", "na"}
        )
    ].copy()

    if data.empty:
        print("No antiSMASH product annotations found; skipping heatmap.")
        return

    # --------------------------------------------------------
    # Keep each complete product cell as ONE category.
    # Comma-separated values are NOT split.
    # --------------------------------------------------------

    counts = pd.crosstab(
        data["product"],
        data["Database"]
    )

    # Preserve dataset order used in the input table.
    datasets = (
        df["Database"]
        .dropna()
        .astype(str)
        .str.strip()
        .loc[lambda s: s.ne("")]
        .drop_duplicates()
        .tolist()
    )

    counts = counts.reindex(
        columns=datasets,
        fill_value=0
    )

    # --------------------------------------------------------
    # Convert counts to percentages independently for
    # each dataset.
    #
    # Each column sums to 100%.
    # --------------------------------------------------------

    totals = counts.sum(axis=0)

    percentages = counts.div(
        totals.replace(0, np.nan),
        axis=1
    ) * 100

    percentages = percentages.fillna(0)

    # --------------------------------------------------------
    # Sort product annotations by their total percentage
    # across datasets, descending.
    # --------------------------------------------------------

    percentages = percentages.loc[
        percentages.sum(axis=1)
        .sort_values(ascending=False)
        .index
    ]

    # --------------------------------------------------------
    # Save the percentage matrix.
    # --------------------------------------------------------

    percentages.to_csv(
        outdir / "antismash_product_heatmap_percentages.tsv",
        sep="\t",
        index_label="Product"
    )

    # --------------------------------------------------------
    # Dynamic figure size.
    # --------------------------------------------------------

    fig, ax = plt.subplots(
        figsize=(
            max(8, 1.2 * len(datasets)),
            max(6, 0.32 * len(percentages))
        )
    )

    # --------------------------------------------------------
    # Heatmap
    # --------------------------------------------------------

    max_percentage = percentages.to_numpy().max()

    sns.heatmap(
        percentages,
        annot=True,
        fmt=".1f",
        cmap="YlGnBu",
        linewidths=0.3,
        linecolor="white",
        vmin=0,
        vmax=max_percentage,
        cbar_kws={
            "label": "Percentage of BGCs (%)"
        },
        ax=ax
    )

    ax.set_xlabel("Dataset")
    ax.set_ylabel("antiSMASH product annotation")
    ax.set_title(
        "antiSMASH product annotations by dataset (%)"
    )

    ax.tick_params(
        axis="x",
        labelrotation=45
    )

    ax.tick_params(
        axis="y",
        labelrotation=0
    )

    savefig(
        fig,
        outdir / "antismash_product_heatmap_percentages.png"
    )


# ============================================================
# Known resistance hit product heatmap
# ============================================================

def plot_resistance_product_heatmap(df, outdir):
    """
    Heatmap of known resistance-hit products by dataset.

    Y-axis: individual resistance-hit products.
    X-axis: datasets.
    Cell values: number of resistance hits.

    Comma-separated values in KnownResistanceHit_product
    are SPLIT into individual hits.

    Example:
        MFS_1,ABC_tran,Acetyltransf_1,MarR

    is counted as:

        MFS_1             -> 1 hit
        ABC_tran          -> 1 hit
        Acetyltransf_1    -> 1 hit
        MarR              -> 1 hit
    """

    data = df.copy()

    data["Database"] = (
        data["Database"]
        .fillna("")
        .astype(str)
        .str.strip()
    )

    data["KnownResistanceHit_product"] = (
        data["KnownResistanceHit_product"]
        .fillna("")
        .astype(str)
        .str.strip()
    )

    # Exclude empty or missing resistance-hit annotations.
    data = data[
        data["Database"].ne("")
        & data["KnownResistanceHit_product"].ne("")
        & ~data["KnownResistanceHit_product"].str.lower().isin(
            {"nan", "none", "na"}
        )
    ].copy()

    if data.empty:
        print(
            "No known resistance-hit products found; "
            "skipping heatmap."
        )
        return

    # --------------------------------------------------------
    # Split comma-separated resistance hits.
    #
    # Example:
    # MFS_1,ABC_tran,Acetyltransf_1,MarR
    #
    # becomes four separate rows/hits.
    # --------------------------------------------------------

    data["KnownResistanceHit_product"] = (
        data["KnownResistanceHit_product"]
        .str.split(",")
    )

    data = data.explode(
        "KnownResistanceHit_product"
    )

    # Remove whitespace introduced around commas.
    data["KnownResistanceHit_product"] = (
        data["KnownResistanceHit_product"]
        .astype(str)
        .str.strip()
    )

    # Remove empty values after splitting.
    data = data[
        data["KnownResistanceHit_product"].ne("")
        & ~data["KnownResistanceHit_product"]
            .str.lower()
            .isin({"nan", "none", "na"})
    ].copy()

    if data.empty:
        print(
            "No valid known resistance-hit products found "
            "after splitting; skipping heatmap."
        )
        return

    # --------------------------------------------------------
    # Count individual resistance hits by dataset.
    # --------------------------------------------------------

    counts = pd.crosstab(
        data["KnownResistanceHit_product"],
        data["Database"]
    )

    # Preserve dataset order used in the input table.
    datasets = (
        df["Database"]
        .dropna()
        .astype(str)
        .str.strip()
        .loc[lambda s: s.ne("")]
        .drop_duplicates()
        .tolist()
    )

    counts = counts.reindex(
        columns=datasets,
        fill_value=0
    )

    # --------------------------------------------------------
    # Sort resistance hits by total occurrence,
    # descending.
    # --------------------------------------------------------

    counts = counts.loc[
        counts.sum(axis=1)
        .sort_values(ascending=False)
        .index
    ]

    # --------------------------------------------------------
    # Save absolute-count matrix.
    # --------------------------------------------------------

    counts.to_csv(
        outdir / "known_resistance_product_heatmap_counts.tsv",
        sep="\t",
        index_label="Resistance_hit_product"
    )

    # --------------------------------------------------------
    # Dynamic figure size.
    # --------------------------------------------------------

    fig, ax = plt.subplots(
        figsize=(
            max(8, 1.2 * len(datasets)),
            max(6, 0.32 * len(counts))
        )
    )

    # --------------------------------------------------------
    # Heatmap
    # --------------------------------------------------------

    sns.heatmap(
        counts,
        annot=True,
        fmt="d",
        cmap="YlGnBu",
        linewidths=0.3,
        linecolor="white",
        cbar_kws={
            "label": "Number of resistance hits"
        },
        ax=ax
    )

    ax.set_xlabel("Dataset")
    ax.set_ylabel("Known resistance-hit product")
    ax.set_title(
        "Known resistance-hit products by dataset"
    )

    ax.tick_params(
        axis="x",
        labelrotation=45
    )

    ax.tick_params(
        axis="y",
        labelrotation=0
    )

    savefig(
        fig,
        outdir / "known_resistance_product_heatmap.png"
    )


# ============================================================
# Known resistance hit product heatmap - percentages
# ============================================================

def plot_resistance_product_heatmap_percentages(df, outdir):
    """
    Heatmap of known resistance-hit products by dataset.

    Y-axis: individual resistance-hit products.
    X-axis: datasets.
    Cell values: percentage of resistance hits represented
                 by each product within each dataset.

    Percentages are calculated independently for each dataset.

    Comma-separated values in KnownResistanceHit_product
    are SPLIT into individual hits.

    Example:
        MFS_1,ABC_tran,Acetyltransf_1,MarR

    is counted as four different resistance hits.
    """

    data = df.copy()

    data["Database"] = (
        data["Database"]
        .fillna("")
        .astype(str)
        .str.strip()
    )

    data["KnownResistanceHit_product"] = (
        data["KnownResistanceHit_product"]
        .fillna("")
        .astype(str)
        .str.strip()
    )

    # Exclude empty or missing resistance-hit annotations.
    data = data[
        data["Database"].ne("")
        & data["KnownResistanceHit_product"].ne("")
        & ~data["KnownResistanceHit_product"].str.lower().isin(
            {"nan", "none", "na"}
        )
    ].copy()

    if data.empty:
        print(
            "No known resistance-hit products found; "
            "skipping percentage heatmap."
        )
        return

    # --------------------------------------------------------
    # Split comma-separated resistance hits.
    # --------------------------------------------------------

    data["KnownResistanceHit_product"] = (
        data["KnownResistanceHit_product"]
        .str.split(",")
    )

    data = data.explode(
        "KnownResistanceHit_product"
    )

    # Remove whitespace around individual hits.
    data["KnownResistanceHit_product"] = (
        data["KnownResistanceHit_product"]
        .astype(str)
        .str.strip()
    )

    # Remove empty values after splitting.
    data = data[
        data["KnownResistanceHit_product"].ne("")
        & ~data["KnownResistanceHit_product"]
            .str.lower()
            .isin({"nan", "none", "na"})
    ].copy()

    if data.empty:
        print(
            "No valid known resistance-hit products found "
            "after splitting; skipping percentage heatmap."
        )
        return

    # --------------------------------------------------------
    # Count individual resistance hits by dataset.
    # --------------------------------------------------------

    counts = pd.crosstab(
        data["KnownResistanceHit_product"],
        data["Database"]
    )

    # Preserve dataset order used in the input table.
    datasets = (
        df["Database"]
        .dropna()
        .astype(str)
        .str.strip()
        .loc[lambda s: s.ne("")]
        .drop_duplicates()
        .tolist()
    )

    counts = counts.reindex(
        columns=datasets,
        fill_value=0
    )

    # --------------------------------------------------------
    # Convert counts to percentages independently for
    # each dataset.
    #
    # Each column sums to 100%.
    # --------------------------------------------------------

    totals = counts.sum(axis=0)

    percentages = counts.div(
        totals.replace(0, np.nan),
        axis=1
    ) * 100

    percentages = percentages.fillna(0)

    # --------------------------------------------------------
    # Sort resistance hits by total percentage across
    # datasets, descending.
    # --------------------------------------------------------

    percentages = percentages.loc[
        percentages.sum(axis=1)
        .sort_values(ascending=False)
        .index
    ]

    # --------------------------------------------------------
    # Save percentage matrix.
    # --------------------------------------------------------

    percentages.to_csv(
        outdir / "known_resistance_product_heatmap_percentages.tsv",
        sep="\t",
        index_label="Resistance_hit_product"
    )

    # --------------------------------------------------------
    # Dynamic figure size.
    # --------------------------------------------------------

    fig, ax = plt.subplots(
        figsize=(
            max(8, 1.2 * len(datasets)),
            max(6, 0.32 * len(percentages))
        )
    )

    # --------------------------------------------------------
    # Use the maximum OBSERVED percentage as the maximum
    # value of the heatmap color scale.
    # --------------------------------------------------------

    max_percentage = percentages.to_numpy().max()

    sns.heatmap(
        percentages,
        annot=True,
        fmt=".1f",
        cmap="YlGnBu",
        linewidths=0.3,
        linecolor="white",
        vmin=0,
        vmax=max_percentage,
        cbar_kws={
            "label": "Percentage of resistance hits (%)"
        },
        ax=ax
    )

    ax.set_xlabel("Dataset")
    ax.set_ylabel("Known resistance-hit product")
    ax.set_title(
        "Known resistance-hit products by dataset (%)"
    )

    ax.tick_params(
        axis="x",
        labelrotation=45
    )

    ax.tick_params(
        axis="y",
        labelrotation=0
    )

    savefig(
        fig,
        outdir / "known_resistance_product_heatmap_percentages.png"
    )



# ============================================================
# Known resistance-hit product -> resistance class mapping
# ============================================================

RESISTANCE_CLASS_MAPPING_TEXT = """
16S ribosomal RNA methyltransferase|16S_rRNA_methyltrans
23S ribosomal RNA methyltransferase|Cfr23_rRNA_methyltrans,Erm23S_rRNA_methyltrans,Erm38,ErmA,ErmB,ErmC
Acetyltransferase|Acetyltransf_4,Acetyltransf_1,Acetyltransf_7,Acetyltransf_8,Acetyltransf_3,Acetyltransf_9
Acyltransferase|Acyltransferase
Aminoglycoside acetyltransferase|AAC3,AAC3-I,AAC6-Ib,AAC6-I,AAC6-II,Antibiotic_NAT
Aminoglycoside nucleotidyltransferase|ANT2,ANT3,ANT4,ANT6,ANT9,ANT
Aminoglycoside phosphotransferase|APH3',APH3,APH6
Aminotransferase|Aminotran_1_2,Aminotran_4
Antibiotic efflux pump|MFS_1,efflux_Bcr_CflA,efflux_EmrB,MATE_efflux,ABC_efflux,adeA-adeI,adeB,adeC-adeK-oprM,adeR,adeS,baeR,baeS,Chlor_Efflux_Pump,emrB,emrE,macA,macB,marA,MexA,MexC,MexE,MexH,MexW-MexI,MexX,MFS_efflux,msbA,norA,phoQ,ramA,RND_efflux,robA,soxR,TetA-B,TetA-G,TetA,TetD,TetE,TetH-TetJ,tet_MFS_efflux,TetY,tolC,RND_mfp,ABC_tran,ABC2_membrane,drrA,ABC1,MFS_3,ACR_tran,Small_Multi_Drug_Res
Beta-lactam resistance|blaI,blaR1,mecR1,TE_Inactivator
Beta-lactamase|Lactamase_B,BCII,BJP,BlaB,CARB-PSE,ClassA,ClassB,ClassC-AmpC,ClassD,CMY-LAT-MOX-ACT-MIR-FOX,CTXM,DHA,DIM-GIM-SIM,Exo,GES,GOB,IMP,IND,KHM,KPC,L1,LRA,MoxA,NDM-CcrA,PC1,Sfh,SHV-LEN,SME,SPM,SubclassB1,SubclassB2,SubclassB3,TEM,VEB-PER,VIM,Beta-lactamase2,Beta-lactamase,Lactamase_B_2,CepA
Cephalosporin resistance|CblA,CfxA
Chloramphenicol acetyltransferase|Chlor_Acetyltrans_CAT,CAT
Chloramphenicol phosphotransferase|Chlor_Phospho_CPT,CPT
Dihydrofolate reductase|DHFR_1
Dioxygenase|Glyoxalase
Fluoroquinolone Resistant DNA Topoisomerase|Fluor_Res_DNA_Topo
Glycopeptide resistance|D_ala_D_ala,Dala_Dala_lig_C,Dala_Dala_lig_N,vanA,vanB,vanC,vanD,vanH,vanR,vanS,vanT,vanW,vanX,vanY,vanZ
Macrolide glycosyltransferase|macrolide_glycosyl
MAR regulator|MarR,MarR_2
Methyltransferase|Methyltransf_18,ArmA_Rmt
PBP transpeptidase|Transpeptidase
Peptide resistance|mprF
Phosphotransferase|APH
Quninolone resistance|Qnr
Ribosomal RNA methyltransferase|FmrO
Tetracycline resistance|TetM-TetW-TetO-TetS,tet_ribosomoal_protect,TetX
Thymidylate synthase|Thymidylat_synt,thym_sym
Transcription factor|Whib,SoxR,romA
Transcription regulator|HTH_AraC
Dihydropteroate synthase|Dihydropteroate
Streptomycin phosphotransferase|APH3''
Aminoglycoside resistance|Aminoglyc_resit
Universal stress protein|Usp
Proteasome subunit|Proteasome
PBP transglycosylase|Transgly
G3P dehydrogenase|Gp_dh_N,Gp_dh_C
Hsp90 protein|HSP90
Carbamoyltransferase|OTCace
DNA gyrase|DNA_gyraseB,DNA_topoisoIV
Biotin-requiring enzyme|Biotin_lipoyl
RNA polymerase|RNA_pol,TIGR02013
Carboxyltransferase|Carboxyl_trans,ACCA
Pentapeptide repeats|Pentapeptide_4
DNA polymerase|TIGR00663
"""

def build_resistance_product_class_mapping():
    """
    Convert the embedded mapping text into a dictionary:
        resistance product -> resistance class
    """
    mapping = {}

    for line in RESISTANCE_CLASS_MAPPING_TEXT.strip().splitlines():
        if not line.strip():
            continue

        resistance_class, products_text = line.split("|", 1)

        resistance_class = resistance_class.strip()

        for product in products_text.split(","):
            product = product.strip()

            if product:
                mapping[product] = resistance_class

    return mapping


RESISTANCE_PRODUCT_TO_CLASS = build_resistance_product_class_mapping()


# ============================================================
# Prepare resistance-class counts
# ============================================================

def prepare_resistance_class_counts(df, outdir):
    """
    Split KnownResistanceHit_product into individual hits,
    map each product to its resistance class, and count
    individual hits per class and dataset.

    Returns a DataFrame:
        Rows    = resistance classes
        Columns = datasets
        Values  = number of hits
    """

    data = df.copy()

    # Clean dataset names.
    data["Database"] = (
        data["Database"]
        .fillna("")
        .astype(str)
        .str.strip()
    )

    # Clean resistance product annotations.
    data["KnownResistanceHit_product"] = (
        data["KnownResistanceHit_product"]
        .fillna("")
        .astype(str)
        .str.strip()
    )

    # Remove missing annotations.
    data = data[
        data["Database"].ne("")
        & data["KnownResistanceHit_product"].ne("")
        & ~data["KnownResistanceHit_product"].str.lower().isin(
            {"nan", "none", "na"}
        )
    ].copy()

    if data.empty:
        print("No known resistance-hit products found.")
        return None

    # Split comma-separated annotations into individual hits.
    data["KnownResistanceHit_product"] = (
        data["KnownResistanceHit_product"].str.split(",")
    )

    data = data.explode("KnownResistanceHit_product")

    data["KnownResistanceHit_product"] = (
        data["KnownResistanceHit_product"]
        .astype(str)
        .str.strip()
    )

    # Remove empty values created during splitting.
    data = data[
        data["KnownResistanceHit_product"].ne("")
        & ~data["KnownResistanceHit_product"].str.lower().isin(
            {"nan", "none", "na"}
        )
    ].copy()

    if data.empty:
        print("No valid resistance hits found after splitting.")
        return None

    # Map each resistance product to its class.
    data["Resistance_class"] = data["KnownResistanceHit_product"].map(
        RESISTANCE_PRODUCT_TO_CLASS
    )

    # Report products that are absent from the embedded mapping.
    unmapped = sorted(
        data.loc[
            data["Resistance_class"].isna(),
            "KnownResistanceHit_product"
        ].unique()
    )

    if unmapped:
        print(
            "Warning: the following resistance products were not "
            "found in the mapping and will retain their original names:"
        )

        for product in unmapped:
            print(f"  - {product}")

        # Keep unmapped products rather than silently discarding them.
        data["Resistance_class"] = data["Resistance_class"].fillna(
            data["KnownResistanceHit_product"].map(
                lambda product: f"Unmapped: {product}"
            )
        )

    # Count hits by resistance class and dataset.
    counts = pd.crosstab(
        data["Resistance_class"],
        data["Database"]
    )

    # Preserve dataset order from the original input table.
    datasets = (
        df["Database"]
        .dropna()
        .astype(str)
        .str.strip()
        .loc[lambda s: s.ne("")]
        .drop_duplicates()
        .tolist()
    )

    counts = counts.reindex(
        columns=datasets,
        fill_value=0
    )

    # Sort classes by total number of hits, descending.
    counts = counts.loc[
        counts.sum(axis=1)
        .sort_values(ascending=False)
        .index
    ]

    return counts


# ============================================================
# Resistance-class heatmap: absolute counts
# ============================================================

def plot_resistance_product_heatmap_class(df, outdir):

    counts = prepare_resistance_class_counts(df, outdir)

    if counts is None or counts.empty:
        print("Skipping resistance-class count heatmap.")
        return

    # Save the count matrix.
    counts.to_csv(
        outdir / "known_resistance_class_heatmap_counts.tsv",
        sep="\t",
        index_label="Resistance_class"
    )

    datasets = counts.columns.tolist()

    fig, ax = plt.subplots(
        figsize=(
            max(8, 1.2 * len(datasets)),
            max(6, 0.32 * len(counts))
        )
    )

    sns.heatmap(
        counts,
        annot=True,
        fmt="d",
        cmap="YlGnBu",
        linewidths=0.3,
        linecolor="white",
        cbar_kws={"label": "Number of resistance hits"},
        ax=ax
    )

    ax.set_xlabel("Dataset")
    ax.set_ylabel("Known resistance-hit class")
    ax.set_title("Known resistance-hit classes by dataset")

    ax.tick_params(axis="x", labelrotation=45)
    ax.tick_params(axis="y", labelrotation=0)

    fig.tight_layout()

    savefig(
        fig,
        outdir / "known_resistance_class_heatmap.png"
    )


# ============================================================
# Resistance-class heatmap: percentages
# ============================================================

def plot_resistance_product_heatmap_percentages_class(df, outdir):

    counts = prepare_resistance_class_counts(df, outdir)

    if counts is None or counts.empty:
        print("Skipping resistance-class percentage heatmap.")
        return

    # Calculate percentages independently for each dataset.
    # Each non-empty dataset column sums to 100%.
    totals = counts.sum(axis=0)

    percentages = counts.div(
        totals.replace(0, np.nan),
        axis=1
    ) * 100

    percentages = percentages.fillna(0)

    # Sort classes by total percentage across datasets.
    percentages = percentages.loc[
        percentages.sum(axis=1)
        .sort_values(ascending=False)
        .index
    ]

    # Save the percentage matrix.
    percentages.to_csv(
        outdir / "known_resistance_class_heatmap_percentages.tsv",
        sep="\t",
        index_label="Resistance_class"
    )

    datasets = percentages.columns.tolist()

    fig, ax = plt.subplots(
        figsize=(
            max(8, 1.2 * len(datasets)),
            max(6, 0.32 * len(percentages))
        )
    )

    # Use the maximum observed percentage as the color-scale maximum.
    max_percentage = percentages.to_numpy().max()

    sns.heatmap(
        percentages,
        annot=True,
        fmt=".1f",
        cmap="YlGnBu",
        linewidths=0.3,
        linecolor="white",
        vmin=0,
        vmax=max_percentage if max_percentage > 0 else 1,
        cbar_kws={"label": "Percentage of resistance hits (%)"},
        ax=ax
    )

    ax.set_xlabel("Dataset")
    ax.set_ylabel("Known resistance-hit class")
    ax.set_title("Known resistance-hit classes by dataset (%)")

    ax.tick_params(axis="x", labelrotation=45)
    ax.tick_params(axis="y", labelrotation=0)

    fig.tight_layout()

    savefig(
        fig,
        outdir / "known_resistance_class_heatmap_percentages.png"
    )




# ============================================================
# UpSet plot: shared resistance-hit classes between datasets
# ============================================================

def plot_resistance_upset(df, outdir, colors):

    from upsetplot import UpSet, from_contents

    # Reuse the existing resistance-class preparation function.
    counts = prepare_resistance_class_counts(df, outdir)

    if counts is None or counts.empty:
        print("Skipping resistance UpSet plot: no resistance hits found.")
        return

    # Build a set of resistance classes for each dataset.
    # A class is included only once per dataset, even if that
    # dataset contains multiple hits belonging to the same class.
    contents = {
        dataset: set(
            counts.index[counts[dataset] > 0].tolist()
        )
        for dataset in counts.columns
    }

    # Remove datasets without resistance classes.
    contents = {
        dataset: classes
        for dataset, classes in contents.items()
        if classes
    }

    if not contents:
        print("Skipping resistance UpSet plot: no classes to compare.")
        return

    # Save the unique resistance classes used in the UpSet plot.
    pd.DataFrame([
        {
            "Database": dataset,
            "Resistance_class": resistance_class
        }
        for dataset, classes in contents.items()
        for resistance_class in sorted(classes)
    ]).to_csv(
        outdir / "resistance_upset_classes.tsv",
        sep="\t",
        index=False
    )

    # Convert dataset -> set of classes into UpSet input.
    upset_data = from_contents(contents)

    # Calculate dynamic figure dimensions.
    n_datasets = len(contents)
    n_intersections = len(upset_data.index.unique())

    figure_width = max(10, 0.8 * n_intersections + 4)
    figure_height = max(6, 0.7 * n_datasets + 3)

    fig = plt.figure(
        figsize=(figure_width, figure_height)
    )

    upset = UpSet(
        upset_data,
        subset_size="count",
        show_counts=True,
        sort_by="cardinality",
        sort_categories_by="-input",
        facecolor="black"
    )

    axes = upset.plot(fig=fig)

    # --------------------------------------------------------
    # Color dataset bars using the shared dataset palette.
    # Match bars to dataset labels using their vertical positions.
    # --------------------------------------------------------
    dataset_ax = axes["totals"]

    # Force Matplotlib to calculate tick positions and labels.
    fig.canvas.draw()

    tick_positions = dataset_ax.get_yticks()
    tick_labels = [
        tick.get_text()
        for tick in dataset_ax.get_yticklabels()
    ]

    # Map each dataset label to its y-axis position.
    label_by_position = {
        position: label
        for position, label in zip(tick_positions, tick_labels)
    }

    for bar in dataset_ax.patches:
        bar_center = bar.get_y() + bar.get_height() / 2

        if not label_by_position:
            continue

        # Find the nearest y-axis tick to this bar's center.
        nearest_position = min(
            label_by_position,
            key=lambda position: abs(position - bar_center)
        )

        dataset = label_by_position[nearest_position]

        if dataset in colors:
            bar.set_facecolor(colors[dataset])

    # Title and output.
    fig.suptitle(
        "Resistance-hit classes shared between datasets",
        fontsize=16,
        y=1.02
    )

    fig.savefig(
        outdir / "resistance_upset.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.close(fig)

    print(
        f"Resistance UpSet plot saved to: "
        f"{outdir / 'resistance_upset.png'}"
    )

    # --------------------------------------------------------
    # Build a table of resistance classes for each intersection.
    # --------------------------------------------------------

    # Each UpSet index entry is a Boolean tuple indicating
    # which datasets participate in that intersection.
    intersection_rows = []

    for intersection in upset_data.index.unique():
        if not isinstance(intersection, tuple):
            intersection = (intersection,)

        participating_datasets = [
            dataset
            for dataset, included in zip(
                upset_data.index.names,
                intersection
            )
            if included
        ]

        # Classes must occur in every participating dataset.
        shared_classes = set.intersection(
            *[
                contents[dataset]
                for dataset in participating_datasets
            ]
        )

        # Exclude classes also found in datasets outside this
        # intersection, because UpSet intersections are exact.
        other_datasets = [
            dataset
            for dataset in contents
            if dataset not in participating_datasets
        ]

        if other_datasets:
            classes_in_other_datasets = set.union(
                *[
                    contents[dataset]
                    for dataset in other_datasets
                ]
            )
            shared_classes -= classes_in_other_datasets

        intersection_rows.append({
            "Datasets": " + ".join(participating_datasets),
            "Number_of_classes": len(shared_classes),
            "Resistance_classes": "; ".join(sorted(shared_classes))
        })

    # Save the intersection-class details to a TSV file.
    intersection_df = pd.DataFrame(intersection_rows)

    intersection_df.to_csv(
        outdir / "resistance_upset_intersection_classes.tsv",
        sep="\t",
        index=False
    )



# ============================================================
# Parse the genes column for regulatory genes
# ============================================================

def parse_gene_records(value):
    """
    Parse the genes column into individual gene records.

    Supports:
      - Multiple records without an outer list.
      - Multiple records enclosed in an outer list.
      - A single gene record.
    """

    if not isinstance(value, str) or not value.strip():
        return []

    value = value.strip()

    # First, try parsing the value as it appears in the TSV.
    try:
        parsed = ast.literal_eval(value)
    except (ValueError, SyntaxError):
        parsed = None

    # Multiple comma-separated records without outer brackets
    # are parsed by Python as a tuple.
    if isinstance(parsed, (list, tuple)):
        # A single gene record looks like:
        # ['contig', 'CDS', start, end, ...]
        if (
            len(parsed) > 1
            and parsed[1] == "CDS"
        ):
            return [list(parsed)]

        # A collection of gene records.
        if all(
            isinstance(record, (list, tuple))
            and len(record) > 1
            and record[1] == "CDS"
            for record in parsed
        ):
            return [list(record) for record in parsed]

    # Fallback: wrap multiple records in an outer list.
    try:
        parsed = ast.literal_eval("[" + value.rstrip(",") + "]")
    except (ValueError, SyntaxError):
        return []

    if not isinstance(parsed, list):
        return []

    if all(
        isinstance(record, (list, tuple))
        and len(record) > 1
        and record[1] == "CDS"
        for record in parsed
    ):
        return [list(record) for record in parsed]

    return []


# ============================================================
# Extract a readable regulatory gene description
# ============================================================

def clean_regulatory_description(annotation):
    """
    Convert a regulatory annotation into its gene class.

    Example:
    regulatory (smcogs) SMCOG1057:TetR family transcriptional
    regulator (Score: 69.8; E-value: 6.5e-21)

    becomes:
    TetR family transcriptional regulator
    """
    description = str(annotation).strip()

    # Remove the regulatory category and annotation source.
    description = re.sub(
        r"^regulatory\s+\([^)]*\)\s*",
        "",
        description,
        flags=re.IGNORECASE
    )

    # Remove the SMCOG identifier, if present.
    description = re.sub(
        r"^SMCOG\d+\s*:\s*",
        "",
        description,
        flags=re.IGNORECASE
    )

    # Remove score and E-value information.
    description = re.sub(
        r"\s*\(Score:.*$",
        "",
        description,
        flags=re.IGNORECASE
    )

    return description.strip().rstrip(".")


# ============================================================
# Extract all regulatory gene occurrences
# ============================================================

def extract_regulatory_gene_occurrences(df, outdir):
    """
    Extract regulatory genes from the genes column.

    Returns one row per regulatory gene occurrence, including
    its dataset, BGC information, locus tag, coordinates,
    original annotation, and cleaned gene class.

    Saves:
        regulatory_gene_occurrences.tsv
    """

    occurrences = []

    for _, row in df.iterrows():

        dataset = row.get("Database", "Unknown")
        gene_records = parse_gene_records(row.get("genes", ""))

        for gene in gene_records:

            # Gene record fields in the supplied genes column:
            # 0 = contig
            # 1 = feature type
            # 2 = start
            # 3 = end
            # 5 = annotations
            # 6 = gene categories
            # 7 = locus tags

            if len(gene) < 7:
                continue

            annotations = gene[5] if isinstance(gene[5], list) else []
            categories = gene[6] if isinstance(gene[6], list) else []

            is_regulatory = any(
                str(category).strip().lower() == "regulatory"
                for category in categories
            )

            if not is_regulatory:
                continue

            # Use regulatory annotations only.
            regulatory_annotations = [
                str(annotation).strip()
                for annotation in annotations
                if str(annotation).strip().lower().startswith(
                    "regulatory "
                )
            ]

            # Keep the gene even if it has a regulatory category
            # but lacks a corresponding detailed annotation.
            if not regulatory_annotations:
                regulatory_annotations = ["regulatory (description unavailable)"]

            locus_tags = (
                gene[7]
                if len(gene) > 7 and isinstance(gene[7], list)
                else []
            )

            locus_tag = ", ".join(map(str, locus_tags)) if locus_tags else ""

            # Avoid counting the same description twice for one gene.
            unique_annotations = list(dict.fromkeys(regulatory_annotations))

            for annotation in unique_annotations:

                gene_class = clean_regulatory_description(annotation)

                if not gene_class:
                    gene_class = "Unknown regulatory gene"

                occurrences.append({
                    "Database": dataset,
                    "OriginalName": row.get("OriginalName", ""),
                    "BGC_contig": row.get("contig", ""),
                    "cluster_number": row.get("cluster_number", ""),
                    "gene_contig": gene[0],
                    "gene_start": gene[2],
                    "gene_end": gene[3],
                    "locus_tag": locus_tag,
                    "regulatory_gene_class": gene_class,
                    "original_annotation": annotation
                })

    occurrence_df = pd.DataFrame(occurrences)

    # Save the complete list of regulatory gene occurrences.
    occurrence_df.to_csv(
        outdir / "regulatory_gene_occurrences.tsv",
        sep="\t",
        index=False
    )

    if occurrence_df.empty:
        print("No regulatory genes were found in the genes column.")
        return None

    print(
        f"Extracted {len(occurrence_df)} regulatory gene occurrences "
        f"from {occurrence_df['Database'].nunique()} datasets."
    )

    return occurrence_df


# ============================================================
# Prepare regulatory gene counts per dataset
# ============================================================

def prepare_regulatory_gene_counts(df, outdir):

    occurrences = extract_regulatory_gene_occurrences(df, outdir)

    if occurrences is None or occurrences.empty:
        return None

    # Count regulatory gene occurrences by class and dataset.
    counts = pd.crosstab(
        occurrences["regulatory_gene_class"],
        occurrences["Database"]
    )

    # Preserve the dataset order from the original DataFrame.
    datasets = (
        df["Database"]
        .dropna()
        .astype(str)
        .str.strip()
        .loc[lambda s: s.ne("")]
        .drop_duplicates()
        .tolist()
    )

    counts = counts.reindex(columns=datasets, fill_value=0)

    # Sort classes by total number of occurrences.
    counts = counts.loc[
        counts.sum(axis=1).sort_values(ascending=False).index
    ]

    return counts


# ============================================================
# Heatmap 1: Regulatory gene counts
# ============================================================

def plot_regulatory_gene_heatmap(df, outdir):

    counts = prepare_regulatory_gene_counts(df, outdir)

    if counts is None or counts.empty:
        print("Skipping regulatory gene count heatmap.")
        return

    counts.to_csv(
        outdir / "regulatory_gene_heatmap_counts.tsv",
        sep="\t",
        index_label="Regulatory_gene_class"
    )

    datasets = counts.columns.tolist()

    fig, ax = plt.subplots(
        figsize=(
            max(8, 1.2 * len(datasets)),
            max(6, 0.35 * len(counts))
        )
    )

    sns.heatmap(
        counts,
        annot=True,
        fmt="d",
        cmap="YlGnBu",
        linewidths=0.3,
        linecolor="white",
        cbar_kws={"label": "Number of regulatory genes"},
        ax=ax
    )

    ax.set_xlabel("Dataset")
    ax.set_ylabel("Regulatory gene class")
    ax.set_title("Regulatory gene classes by dataset")

    ax.tick_params(axis="x", labelrotation=45)
    ax.tick_params(axis="y", labelrotation=0)

    fig.tight_layout()

    savefig(
        fig,
        outdir / "regulatory_gene_heatmap.png"
    )


# ============================================================
# Heatmap 2: Regulatory gene percentages
# ============================================================

def plot_regulatory_gene_heatmap_percentages(df, outdir):

    counts = prepare_regulatory_gene_counts(df, outdir)

    if counts is None or counts.empty:
        print("Skipping regulatory gene percentage heatmap.")
        return

    # Calculate each class as a percentage of all regulatory
    # gene occurrences in its dataset.
    totals = counts.sum(axis=0)

    percentages = counts.div(
        totals.replace(0, np.nan),
        axis=1
    ) * 100

    percentages = percentages.fillna(0)

    # Sort by total percentage across datasets.
    percentages = percentages.loc[
        percentages.sum(axis=1)
        .sort_values(ascending=False)
        .index
    ]

    percentages.to_csv(
        outdir / "regulatory_gene_heatmap_percentages.tsv",
        sep="\t",
        index_label="Regulatory_gene_class"
    )

    datasets = percentages.columns.tolist()
    max_percentage = percentages.to_numpy().max()

    fig, ax = plt.subplots(
        figsize=(
            max(8, 1.2 * len(datasets)),
            max(6, 0.35 * len(percentages))
        )
    )

    sns.heatmap(
        percentages,
        annot=True,
        fmt=".1f",
        cmap="YlGnBu",
        linewidths=0.3,
        linecolor="white",
        vmin=0,
        vmax=max_percentage if max_percentage > 0 else 1,
        cbar_kws={"label": "Regulatory genes (%)"},
        ax=ax
    )

    ax.set_xlabel("Dataset")
    ax.set_ylabel("Regulatory gene class")
    ax.set_title("Regulatory gene classes by dataset (%)")

    ax.tick_params(axis="x", labelrotation=45)
    ax.tick_params(axis="y", labelrotation=0)

    fig.tight_layout()

    savefig(
        fig,
        outdir / "regulatory_gene_heatmap_percentages.png"
    )



# ============================================================
# UpSet plot: shared regulatory gene classes between datasets
# ============================================================

def plot_regulatory_upset(df, outdir, colors):

    from upsetplot import UpSet, from_contents
    from matplotlib.patches import Patch

    # Reuse the existing extraction function.
    # This also saves regulatory_gene_occurrences.tsv.
    occurrences = extract_regulatory_gene_occurrences(df, outdir)

    if occurrences is None or occurrences.empty:
        print("Skipping regulatory UpSet plot: no regulatory genes found.")
        return

    # Build a set of unique regulator classes for each dataset.
    # count only once per dataset for the intersection analysis.
    contents = (
        occurrences
        .groupby("Database")["regulatory_gene_class"]
        .apply(lambda values: set(values.dropna()))
        .to_dict()
    )

    # Remove datasets with no regulator classes.
    contents = {
        dataset: classes
        for dataset, classes in contents.items()
        if classes
    }

    if not contents:
        print("Skipping regulatory UpSet plot: no classes to compare.")
        return

    # Save the unique regulator classes used in the UpSet plot.
    pd.DataFrame([
        {
            "Database": dataset,
            "regulatory_gene_class": gene_class
        }
        for dataset, classes in contents.items()
        for gene_class in sorted(classes)
    ]).to_csv(
        outdir / "regulatory_upset_classes.tsv",
        sep="\t",
        index=False
    )



    # Convert dataset -> set of classes into UpSet input.
    upset_data = from_contents(contents)
    # Count the number of datasets.
    n_datasets = len(contents)
    # Count the number of intersections represented in the UpSet data.
    n_intersections = len(upset_data.index.unique())
    # Dynamic figure dimensions. Height increases with the number of datasets. Width increases with the number of intersections.
    figure_width = max(10, 0.8 * n_intersections + 4)
    figure_height = max(6, 0.7 * n_datasets + 3)

    fig = plt.figure(figsize=(figure_width, figure_height))

    upset = UpSet(
        upset_data,
        subset_size="count",
        show_counts=True,
        sort_by="cardinality",
        sort_categories_by="-input",
        facecolor="black"
    )

    axes = upset.plot(fig=fig)

    # --------------------------------------------------------
    # Color the dataset bars using the shared dataset palette.
    # --------------------------------------------------------
    dataset_ax = axes["totals"]

    # UpSet orders categories according to its internal sorting.
    # Match each bar to its dataset using the y-axis tick labels.
    tick_labels = [
        tick.get_text()
        for tick in dataset_ax.get_yticklabels()
    ]

    for bar, label in zip(dataset_ax.patches, tick_labels):
        if label in colors:
            bar.set_facecolor(colors[label])

    fig.suptitle(
        "Regulatory gene classes shared between datasets",
        fontsize=16,
        y=1.02
    )

    fig.savefig(
        outdir / "regulatory_upset.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.close(fig)

    print(
        f"Regulatory UpSet plot saved to: "
        f"{outdir / 'regulatory_upset.png'}"
    )

    # --------------------------------------------------------
    # Build a table of regulatory classes for each intersection.
    # --------------------------------------------------------

    # Each UpSet index entry is a Boolean tuple indicating
    # which datasets participate in that intersection.
    intersection_rows = []

    for intersection in upset_data.index.unique():
        if not isinstance(intersection, tuple):
            intersection = (intersection,)

        participating_datasets = [
            dataset
            for dataset, included in zip(
                upset_data.index.names,
                intersection
            )
            if included
        ]

        # Classes must occur in every participating dataset.
        shared_classes = set.intersection(
            *[
                contents[dataset]
                for dataset in participating_datasets
            ]
        )

        # Exclude classes also found in datasets outside this
        # intersection, because UpSet intersections are exact.
        other_datasets = [
            dataset
            for dataset in contents
            if dataset not in participating_datasets
        ]

        if other_datasets:
            classes_in_other_datasets = set.union(
                *[
                    contents[dataset]
                    for dataset in other_datasets
                ]
            )
            shared_classes -= classes_in_other_datasets

        intersection_rows.append({
            "Datasets": " + ".join(participating_datasets),
            "Number_of_classes": len(shared_classes),
            "Regulatory_classes": "; ".join(sorted(shared_classes))
        })

    # Save the intersection-class details to a TSV file.
    intersection_df = pd.DataFrame(intersection_rows)

    intersection_df.to_csv(
        outdir / "regulatory_upset_intersection_classes.tsv",
        sep="\t",
        index=False
    )


# ============================================================
# General summary
# ============================================================

def write_summary(df, outdir):
    rows = []

    for db, g in df.groupby("Database"):
        rows.append({
            "Database": db,
            "BGC_count": len(g),
            "Complete_BGCs": (g["completeness"].astype(str).str.lower() == "complete").sum(),
            "Fragmented_BGCs": (g["completeness"].astype(str).str.lower() == "fragmented").sum(),
            "Resistance_BGCs": g["KnownResistanceHit_product"].notna().sum(),
            "Regulatory_BGCs": g["regulatory_genes"].sum(),
            "Similarity_values": g[SIM_COL].notna().sum(),
            "Mean_similarity": g[SIM_COL].mean(),
            "Median_similarity": g[SIM_COL].median()
        })

    pd.DataFrame(rows).to_csv(
        outdir / "dataset_summary.csv",
        index=False
    )


# ============================================================
# BiG-SCAPE
# ============================================================

def plot_bigscape(df, outdir):
    for col, stem, title in [
        ("BiGSCAPE_Category", "bigscape_category_percentage", "BiG-SCAPE category composition"),
        ("BiGSCAPE_Class", "bigscape_class_percentage", "BiG-SCAPE class composition")
    ]:
        d = df.dropna(subset=[col])

        if d.empty:
            continue

        counts = pd.crosstab(d["Database"], d[col])
        pct = counts.div(counts.sum(axis=1), axis=0) * 100

        ax = pct.plot(kind="bar", stacked=True, figsize=(12, 7))
        ax.set_xlabel("Database")
        ax.set_ylabel("BGCs (%)")
        ax.set_title(title)
        ax.legend(title=col, bbox_to_anchor=(1.02, 1), loc="upper left")

        savefig(ax.figure, outdir / f"{stem}.png")
        pct.round(3).to_csv(outdir / f"{stem}.csv")


def plot_product_prediction(df, outdir):
    d = df.dropna(subset=["product"])

    if d.empty:
        return

    rows = []

    for _, r in d.iterrows():
        for product in split_cell(r["product"]):
            rows.append({
                "Database": r["Database"],
                "Product": product
            })

    x = pd.DataFrame(rows)

    if x.empty:
        return

    counts = pd.crosstab(x["Database"], x["Product"])
    pct = counts.div(counts.sum(axis=1), axis=0) * 100

    ax = pct.plot(kind="bar", stacked=True, figsize=(12, 7))
    ax.set_xlabel("Database")
    ax.set_ylabel("Product predictions (%)")
    ax.set_title("Product Prediction by database")
    ax.legend(title="Product", bbox_to_anchor=(1.02, 1), loc="upper left")

    savefig(
        ax.figure,
        outdir / "product_prediction_percentage.png"
    )

    pct.round(3).to_csv(
        outdir / "product_prediction_percentage.csv"
    )


# ============================================================
# Resistance
# ============================================================

def resistance_long(df):
    rows = []

    for idx, r in df.iterrows():
        vals = split_cell(r["KnownResistanceHit_product"])

        if not vals:
            vals = split_cell(r["KnownResistanceHit_Resfam"])

        for value in vals:
            rows.append({
                "row_id": idx,
                "Database": r["Database"],
                "Resistance": value
            })

    return pd.DataFrame(rows)


def plot_resistance(df, outdir):
    long = resistance_long(df)

    if long.empty:
        return

    counts = pd.crosstab(long["Database"], long["Resistance"])
    pct = counts.div(counts.sum(axis=1), axis=0) * 100

    fig, ax = plt.subplots(figsize=(12, max(5, 0.35 * len(pct.columns))))
    sns.heatmap(pct, annot=True, fmt=".1f", cmap="viridis", ax=ax)
    ax.set_xlabel("Resistance-associated hit")
    ax.set_ylabel("Database")
    ax.set_title("Resistance hit composition by database")

    savefig(fig, outdir / "resistance_heatmap.png")

    ax = pct.plot(kind="bar", stacked=True, figsize=(12, 7))
    ax.set_xlabel("Database")
    ax.set_ylabel("Resistance hits (%)")
    ax.set_title("Resistance hit composition")
    ax.legend(title="Resistance", bbox_to_anchor=(1.02, 1), loc="upper left")

    savefig(ax.figure, outdir / "resistance_stacked_bar.png")

    pct.round(3).to_csv(outdir / "resistance_percentage.csv")

    binary = pd.crosstab(long["row_id"], long["Resistance"]).clip(upper=1)

    binary = binary.reindex(df.index, fill_value=0)

    binary.insert(0, "Database", df["Database"].values)

    binary.insert(1, "BGC_index", df.index)

    binary.to_csv(outdir / "resistance_bgc_binary_matrix.csv",index=False)

    plot_upset(binary.drop(columns=["Database", "BGC_index"]), "Resistance hit combinations", outdir / "resistance_upset.png")


# ============================================================
# Regulatory genes
# ============================================================

def extract_regulatory_families(value):
    if pd.isna(value):
        return []

    try:
        data = ast.literal_eval(str(value))
    except Exception:
        return []

    if not isinstance(data, list):
        return []

    found = []

    for gene in data:
        if not isinstance(gene, list) or len(gene) < 7:
            continue

        categories = gene[6]

        if not isinstance(categories, list):
            continue

        if not any("regulatory" in str(x).lower() for x in categories):
            continue

        annotations = (
            gene[5]
            if len(gene) > 5 and isinstance(gene[5], list)
            else []
        )

        for annotation in annotations:
            s = str(annotation)

            m = re.search(
                r"SMCOG\d+:([^\(]+)",
                s
            )

            if m:
                found.append(m.group(1).strip())
            else:
                s = re.sub(
                    r"^.*?regulatory\s*\([^)]*\)\s*",
                    "",
                    s,
                    flags=re.I
                )

                if s:
                    found.append(
                        s.split(";")[0].strip()
                    )

    return list(dict.fromkeys(found))


def regulatory_long(df):
    rows = []

    for idx, r in df.iterrows():
        families = extract_regulatory_families(r["genes"])

        if not families and r["regulatory_genes"]:
            families = ["Regulatory gene"]

        for family in families:
            rows.append({
                "row_id": idx,
                "Database": r["Database"],
                "Regulatory_gene": family
            })

    return pd.DataFrame(rows)


def plot_regulatory(df, outdir):
    long = regulatory_long(df)

    if long.empty:
        return

    counts = pd.crosstab(
        long["Database"],
        long["Regulatory_gene"]
    )

    pct = counts.div(
        counts.sum(axis=1),
        axis=0
    ) * 100

    fig, ax = plt.subplots(
        figsize=(12, max(5, 0.35 * len(pct.columns)))
    )

    sns.heatmap(
        pct,
        annot=True,
        fmt=".1f",
        cmap="magma",
        ax=ax
    )

    ax.set_xlabel("Regulatory gene family")
    ax.set_ylabel("Database")
    ax.set_title("Regulatory gene families by database")

    savefig(
        fig,
        outdir / "regulatory_heatmap.png"
    )

    ax = pct.plot(
        kind="bar",
        stacked=True,
        figsize=(12, 7)
    )

    ax.set_xlabel("Database")
    ax.set_ylabel("Regulatory gene families (%)")
    ax.set_title("Regulatory gene family composition")

    ax.legend(
        title="Regulatory gene",
        bbox_to_anchor=(1.02, 1),
        loc="upper left"
    )

    savefig(
        ax.figure,
        outdir / "regulatory_stacked_bar.png"
    )

    pct.round(3).to_csv(
        outdir / "regulatory_percentage.csv"
    )

    binary = pd.crosstab(
        long["row_id"],
        long["Regulatory_gene"]
    ).clip(upper=1)

    binary = binary.reindex(
        df.index,
        fill_value=0
    )

    binary.insert(
        0,
        "Database",
        df["Database"].values
    )

    binary.insert(
        1,
        "BGC_index",
        df.index
    )

    binary.to_csv(
        outdir / "regulatory_bgc_binary_matrix.csv",
        index=False
    )

    plot_upset(
        binary.drop(columns=["Database", "BGC_index"]),
        "Regulatory gene family combinations",
        outdir / "regulatory_upset.png"
    )


# ============================================================
# UpSet
# ============================================================

def plot_upset(binary, title, path, max_sets=12):
    if binary.empty or binary.shape[1] == 0:
        return

    cols = (
        binary.sum()
        .sort_values(ascending=False)
        .head(max_sets)
        .index
        .tolist()
    )

    b = binary[cols].astype(int)
    combos = []

    for bits in itertools.product([0, 1], repeat=len(cols)):
        mask = np.ones(len(b), dtype=bool)

        for col, bit in zip(cols, bits):
            mask &= b[col].to_numpy() == bit

        n = int(mask.sum())

        if n and sum(bits):
            combos.append((bits, n))

    combos.sort(
        key=lambda x: x[1],
        reverse=True
    )

    combos = combos[:20]

    if not combos:
        return

    x = np.arange(len(combos))

    fig, (axbar, axmat) = plt.subplots(
        2,
        1,
        figsize=(max(9, 0.55 * len(combos)), 8),
        gridspec_kw={"height_ratios": [2.2, 3]},
        sharex=True
    )

    axbar.bar(
        x,
        [n for _, n in combos]
    )

    axbar.set_ylabel("BGC count")
    axbar.set_title(title)

    for i, (bits, _) in enumerate(combos):
        active = [j for j, bit in enumerate(bits) if bit]

        for j in range(len(cols)):
            axmat.scatter(
                i,
                j,
                s=35 if j in active else 18,
                alpha=1 if j in active else 0.18
            )

        if len(active) > 1:
            axmat.plot(
                [i, i],
                [min(active), max(active)],
                linewidth=1
            )

    axmat.set_yticks(range(len(cols)))
    axmat.set_yticklabels(cols)
    axmat.set_ylabel("Set")
    axmat.set_xlabel("Intersection")

    savefig(
        fig,
        path
    )


# ============================================================
# Similarity
# ============================================================

def plot_similarity(df, outdir, colors):
    d = df.dropna(subset=[SIM_COL]).copy()

    if d.empty:
        return

    order = (
        d["Database"]
        .drop_duplicates()
        .tolist()
    )

    palette = {
        db: colors.get(
            db,
            (0.5, 0.5, 0.5)
        )
        for db in order
    }

    # Boxplot.
    fig, ax = plt.subplots(
        figsize=(max(6, 1.2 * len(order)), 6)
    )

    sns.boxplot(
        data=d,
        x="Database",
        y=SIM_COL,
        order=order,
        hue="Database",
        palette=palette,
        legend=False,
        ax=ax
    )

    sns.stripplot(
        data=d,
        x="Database",
        y=SIM_COL,
        order=order,
        color="black",
        alpha=0.65,
        jitter=True,
        ax=ax
    )

    ax.set_title(
        "BGC hit mean similarity by database"
    )

    ax.set_ylabel(
        "Mean similarity (%)"
    )

    savefig(
        fig,
        outdir / "similarity_boxplot.png"
    )

    # Histogram.
    fig, ax = plt.subplots(
        figsize=(9, 6)
    )

    for db, g in d.groupby(
        "Database",
        sort=False
    ):
        ax.hist(
            g[SIM_COL],
            bins="auto",
            alpha=0.55,
            label=db,
            color=colors.get(
                db,
                (0.5, 0.5, 0.5)
            )
        )

    ax.set_xlabel(
        "Mean similarity (%)"
    )

    ax.set_ylabel(
        "BGC count"
    )

    ax.set_title(
        "Distribution of BGC hit mean similarities"
    )

    ax.legend()

    savefig(
        fig,
        outdir / "similarity_histogram.png"
    )

    # Q-Q plot.
    fig, ax = plt.subplots(
        figsize=(7, 7)
    )

    stats.probplot(
        d[SIM_COL],
        dist="norm",
        plot=ax
    )

    ax.set_title(
        "Q-Q plot: BGC hit mean similarity"
    )

    savefig(
        fig,
        outdir / "similarity_qqplot.png"
    )

    d[
        ["Database", SIM_COL]
    ].to_csv(
        outdir / "similarity_values.csv",
        index=False
    )

    similarity_statistics(
        d,
        outdir
    )


def bh_adjust(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    result = np.empty(n)
    result[order] = np.minimum(ranked, 1)
    return result


def games_howell(groups):
    if len(groups) < 2 or not hasattr(stats, "studentized_range"):
        return pd.DataFrame()

    rows = []

    for a, b in itertools.combinations(groups, 2):
        x = np.asarray(groups[a], dtype=float)
        y = np.asarray(groups[b], dtype=float)

        x = x[np.isfinite(x)]
        y = y[np.isfinite(y)]

        if len(x) < 2 or len(y) < 2:
            continue

        nx, ny = len(x), len(y)
        mx, my = np.mean(x), np.mean(y)
        vx, vy = np.var(x, ddof=1), np.var(y, ddof=1)

        se = math.sqrt(vx / nx + vy / ny)

        if se == 0:
            continue

        q = abs(mx - my) / se

        df_num = (vx / nx + vy / ny) ** 2

        df_den = (
            vx ** 2 / (nx ** 2 * (nx - 1))
            + vy ** 2 / (ny ** 2 * (ny - 1))
        )

        df_welch = df_num / df_den if df_den else np.nan

        p = stats.studentized_range.sf(
            q * math.sqrt(2),
            len(groups),
            df_welch
        )

        rows.append({
            "group1": a,
            "group2": b,
            "n1": nx,
            "n2": ny,
            "mean1": mx,
            "mean2": my,
            "difference": mx - my,
            "q": q,
            "df": df_welch,
            "pvalue": p
        })

    result = pd.DataFrame(rows)

    if not result.empty:
        result["pvalue_BH"] = bh_adjust(
            result["pvalue"]
        )

    return result


def similarity_statistics(df, outdir):
    groups = {
        db: g[SIM_COL].dropna().to_numpy()
        for db, g in df.groupby("Database")
    }

    groups = {
        k: v
        for k, v in groups.items()
        if len(v)
    }

    rows = []

    if len(groups) >= 2:
        h, p = stats.kruskal(*groups.values())

        rows.append({
            "test": "Kruskal-Wallis",
            "statistic": h,
            "pvalue": p,
            "comparison": ""
        })

        for a, b in itertools.combinations(groups, 2):
            u, p = stats.mannwhitneyu(
                groups[a],
                groups[b],
                alternative="two-sided"
            )

            rows.append({
                "test": "Mann-Whitney U",
                "statistic": u,
                "pvalue": p,
                "comparison": f"{a} vs {b}"
            })

    pd.DataFrame(rows).to_csv(
        outdir / "similarity_statistics.csv",
        index=False
    )

    gh = games_howell(groups)

    if not gh.empty:
        gh.to_csv(
            outdir / "games_howell.csv",
            index=False
        )


# ============================================================
# Resistance dendrogram
# ============================================================

def plot_resistance_dendrogram(df, outdir):
    path = outdir / "resistance_bgc_binary_matrix.csv"

    if not path.exists():
        return

    binary = pd.read_csv(path)

    meta = binary[
        ["BGC_index", "Database"]
    ]

    x = binary.drop(
        columns=["BGC_index", "Database"],
        errors="ignore"
    )

    if x.shape[0] < 2 or x.shape[1] < 1:
        return

    x = x.loc[:, x.nunique(dropna=False) > 1]

    if x.shape[1] == 0:
        return

    dist = pdist(
        x.to_numpy(dtype=float),
        metric="jaccard"
    )

    z = linkage(
        dist,
        method="average"
    )

    fig, ax = plt.subplots(
        figsize=(11, max(6, 0.35 * len(x)))
    )

    dendrogram(
        z,
        labels=[
            f"{db} | BGC {idx}"
            for db, idx in zip(
                meta["Database"],
                meta["BGC_index"]
            )
        ],
        orientation="right",
        ax=ax
    )

    ax.set_title(
        "Resistance-profile dendrogram"
    )

    ax.set_xlabel(
        "Jaccard distance"
    )

    savefig(
        fig,
        outdir / "resistance_dendrogram.png"
    )


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description="Generate consolidated plots from Krill results."
    )

    parser.add_argument(
        "input",
        type=Path,
        help="DBs_BGCs_with_Hits_BiGSCAPE.tsv"
    )

    parser.add_argument(
        "--input-path",
        type=Path,
        default=None,
        help="Krill input directory. If omitted, inferred automatically."
    )

    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("KrillPlots"),
        help="Output directory."
    )

    args = parser.parse_args()

    outdir = args.output.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    sns.set_theme(style="whitegrid")

    # ------------------------------------------------------------
    # Main input
    # ------------------------------------------------------------

    df = read_input(
        args.input
    )

    # ------------------------------------------------------------
    # Dynamic dataset colors
    # ------------------------------------------------------------

    DATASET_COLORS = create_dataset_colors(
        df
    )

    print(
        "\nDatasets detected:"
    )

    for db, color in DATASET_COLORS.items():
        print(
            f"  {db}: {matplotlib.colors.to_hex(color)}"
        )

    # ------------------------------------------------------------
    # Output directories
    # ------------------------------------------------------------

    product_dir = outdir / "01_products_bigscape"
    resistance_dir = outdir / "02_resistance"
    regulatory_dir = outdir / "03_regulatory"
    similarity_dir = outdir / "04_similarity"

    for directory in [
        product_dir,
        resistance_dir,
        regulatory_dir,
        similarity_dir
    ]:
        directory.mkdir(
            parents=True,
            exist_ok=True
        )

    # ------------------------------------------------------------
    # General summary
    # ------------------------------------------------------------

    write_summary(
        df,
        outdir
    )

    # ------------------------------------------------------------
    # Database size / BGC / BGC density / MIBIG similarity
    # ------------------------------------------------------------

    krill_input = find_input_path(
        args.input,
        args.input_path
    )

    plot_database_summary(
        args.input,
        krill_input,
        outdir,
        df,
        DATASET_COLORS
    )

    plot_cumulative_percentage_bars(
        df,
        outdir,
        DATASET_COLORS
    )

    plot_size_violin_plots(
        df,
        outdir,
        DATASET_COLORS
    )

    plot_similarity_violin_plots(
        df,
        outdir,
        DATASET_COLORS
    )

    plot_product_and_bigscape_summaries(
        df=df,
        outdir=outdir,
    )

    plot_product_heatmap(df, outdir)

    plot_product_heatmap_percentages(df, outdir)

    plot_resistance_product_heatmap(df, outdir)

    plot_resistance_product_heatmap_percentages(df, outdir)

    plot_resistance_product_heatmap_class(df, outdir)

    plot_resistance_product_heatmap_percentages_class(df, outdir)

    plot_regulatory_gene_heatmap(df, outdir)

    plot_regulatory_gene_heatmap_percentages(df, outdir)

    plot_regulatory_upset(df, outdir, DATASET_COLORS)

    plot_resistance_upset(df, outdir, DATASET_COLORS)

    # ------------------------------------------------------------
    # BiG-SCAPE
    # ------------------------------------------------------------

    plot_bigscape(
        df,
        product_dir
    )

    plot_product_prediction(
        df,
        product_dir
    )

    # ------------------------------------------------------------
    # Resistance
    # ------------------------------------------------------------

    plot_resistance(
        df,
        resistance_dir
    )

    plot_resistance_dendrogram(
        df,
        resistance_dir
    )

    # ------------------------------------------------------------
    # Regulatory genes
    # ------------------------------------------------------------

    plot_regulatory(
        df,
        regulatory_dir
    )

    # ------------------------------------------------------------
    # Similarity
    # ------------------------------------------------------------

    plot_similarity(
        df,
        similarity_dir,
        DATASET_COLORS
    )

    print(
        "\nAnalysis completed."
    )

    print(
        f"Input: {args.input}"
    )

    print(
        f"Databases: {df['Database'].nunique()}"
    )

    print(
        f"BGCs: {len(df)}"
    )

    print(
        f"Output: {outdir}\n"
    )


if __name__ == "__main__":
    main()