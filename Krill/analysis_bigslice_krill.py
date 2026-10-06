#!/usr/bin/env python3
"""
Analyse BiG-SLiCE results for Krill.

This script is a standalone conversion of analysis.ipynb. It extracts BGC
features, taxonomy and metadata from the BiG-SLiCE SQLite database; obtains GCF
models at a selected clustering threshold; clusters GCFs into GCC-like bins
with K-means; assigns BGCs to the nearest centroid; creates a Newick cladogram;
writes the supplementary GCC metadata and five iTOL annotation files; and
performs the notebook's dataset-, taxonomy-, and dataset-description-based
exploratory analyses, followed by a consolidated taxonomy/HMM feature profile.

Run from the Krill project root, or pass --project-root.
"""

from __future__ import annotations

import argparse
import os
import sqlite3
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors
from sklearn.cluster import KMeans, AgglomerativeClustering
from scipy.cluster import hierarchy
from ete3 import PhyloTree
from matplotlib import pyplot as plt
from matplotlib import ticker, rc
from random import randint
import random
import pickle



def render_publication_tree(tree, metadata: pd.DataFrame, output_dir: Path, dpi: int = 600,
                           label_font_size: float = 4.0, show_labels: bool = True) -> None:
    """Render a rectangular GCC tree with metadata tracks as PNG, PDF and SVG."""
    output_dir.mkdir(parents=True, exist_ok=True)
    leaves = tree.get_leaves()
    if not leaves:
        raise ValueError("The GCC tree contains no leaves.")

    # Order leaves as represented by the rooted tree and calculate branch lengths.
    leaf_y = {leaf.name: float(i) for i, leaf in enumerate(leaves)}
    node_y = {}
    node_x = {}
    for node in tree.traverse("postorder"):
        if node.is_leaf():
            node_y[node] = leaf_y[node.name]
        else:
            child_ys = [node_y[ch] for ch in node.children]
            node_y[node] = sum(child_ys) / len(child_ys)
    node_x[tree] = 0.0
    for node in tree.traverse("preorder"):
        for child in node.children:
            node_x[child] = node_x[node] + max(float(child.dist or 0.0), 0.0)

    # Match GCC leaves to metadata rows; the dummy outgroup has no metadata row.
    track_specs = [
        ("bgc_count", "BGC count", "viridis"),
        ("avg_signal_strength", "Mean feature signal", "Reds"),
        ("avg_bgc_length", "Mean BGC length (bp)", "Greys"),
        ("avg_gcf_size", "Mean GCF size", "Blues"),
    ]
    available_tracks = [(c, label, cmap) for c, label, cmap in track_specs if c in metadata.columns]
    class_names = ["Polyketide", "NRP", "RiPP", "Saccharide", "Terpene", "Alkaloid", "Other"]
    class_cols = ["class_" + name for name in class_names if "class_" + name in metadata.columns]
    dataset_cols = [c for c in metadata.columns if c.startswith("dataset_")]

    n = len(leaves)
    fig_height = max(8.0, min(30.0, n * 0.035))
    track_count = len(available_tracks) + (1 if class_cols else 0) + (1 if dataset_cols else 0)
    fig_width = 16 if track_count else 10
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    max_x = max(node_x.values()) or 1.0
    tree_width = max_x * 0.68

    # Draw rectangular branches.
    for node in tree.traverse("preorder"):
        if node.children:
            ys = [node_y[ch] for ch in node.children]
            ax.plot([node_x[node], node_x[node]], [min(ys), max(ys)], color="#333333", lw=0.45)
            for child in node.children:
                ax.plot([node_x[node], node_x[child]], [node_y[child], node_y[child]], color="#333333", lw=0.45)

    if show_labels:
        for leaf in leaves:
            ax.text(node_x[leaf] + max_x * 0.008, node_y[leaf], leaf.name,
                    va="center", ha="left", fontsize=label_font_size, color="#222222")

    track_start = max_x + max_x * (0.25 if show_labels else 0.08)
    track_gap = max_x * 0.025
    track_width = max_x * 0.07
    # Leaf-aligned scalar tracks.
    for i, (column, label, cmap) in enumerate(available_tracks):
        x0 = track_start + i * (track_width + track_gap)
        values = []
        ys = []
        for leaf in leaves:
            if leaf.name.startswith("GCC"):
                idx = int(leaf.name[3:])
                if idx in metadata.index and column in metadata.columns:
                    values.append(float(metadata.loc[idx, column]))
                    ys.append(node_y[leaf])
        if values:
            vmax = max(values) or 1.0
            norm = plt.Normalize(vmin=min(0.0, min(values)), vmax=vmax)
            cmap_obj = plt.get_cmap(cmap)
            for value, y in zip(values, ys):
                ax.add_patch(plt.Rectangle((x0, y - 0.38), track_width, 0.76,
                                           color=cmap_obj(norm(value)), linewidth=0))
            ax.text(x0 + track_width / 2, -1.5, label, ha="center", va="bottom",
                    rotation=35, fontsize=7)

    x0 = track_start + len(available_tracks) * (track_width + track_gap)
    if class_cols:
        colors = {"class_Polyketide": "#ff872b", "class_NRP": "#1cb226", "class_RiPP": "#d64859",
                  "class_Saccharide": "#0101ad", "class_Terpene": "#ae00ae", "class_Alkaloid": "#cccc00",
                  "class_Other": "#7e7e7f"}
        for leaf in leaves:
            if not leaf.name.startswith("GCC"):
                continue
            idx = int(leaf.name[3:])
            if idx not in metadata.index:
                continue
            vals = np.array([float(metadata.loc[idx, c]) for c in class_cols])
            total = vals.sum()
            if total <= 0:
                continue
            left = x0
            for c, value in zip(class_cols, vals / total):
                width = track_width * value
                ax.add_patch(plt.Rectangle((left, node_y[leaf] - 0.38), width, 0.76,
                                           color=colors.get(c, "#999999"), linewidth=0))
                left += width
        ax.text(x0 + track_width / 2, -1.5, "Chemical classes", ha="center", va="bottom",
                rotation=35, fontsize=7)
        x0 += track_width + track_gap

    if dataset_cols:
        dataset_colors = plt.get_cmap("tab20")
        for leaf in leaves:
            if not leaf.name.startswith("GCC"):
                continue
            idx = int(leaf.name[3:])
            if idx not in metadata.index:
                continue
            vals = np.array([float(metadata.loc[idx, c]) for c in dataset_cols])
            total = vals.sum()
            if total <= 0:
                continue
            left = x0
            for j, value in enumerate(vals / total):
                width = track_width * value
                ax.add_patch(plt.Rectangle((left, node_y[leaf] - 0.38), width, 0.76,
                                           color=dataset_colors(j % 20), linewidth=0))
                left += width
        ax.text(x0 + track_width / 2, -1.5, "Datasets", ha="center", va="bottom",
                rotation=35, fontsize=7)

    ax.set_xlim(-max_x * 0.02, track_start + (track_count + 1) * (track_width + track_gap))
    ax.set_ylim(n - 0.5, -3.0)
    ax.set_yticks([])
    ax.set_xlabel("Branch length (Euclidean distance; Ward linkage)")
    ax.set_title("BiG-SLiCE GCC hierarchical tree", loc="left", fontsize=12, weight="bold")
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.spines["bottom"].set_color("#777777")
    fig.tight_layout()
    for ext in ("png", "pdf", "svg"):
        fig.savefig(output_dir / f"figure_6_tree.{ext}", dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Analyse BiG-SLiCE results and generate GCC tree/metadata for Krill."
    )
    parser.add_argument(
        "--project-root", type=Path, default=Path.cwd(),
        help="Krill project directory containing source_files/ (default: current directory).",
    )
    parser.add_argument(
        "--db", type=Path, default=None,
        help="Optional path to BiG-SLiCE data.db; overrides the default under project root.",
    )
    parser.add_argument(
        "--bigslice-analysis-threshold", "--threshold", dest="threshold", type=float, default=None,
        help="BiG-SLiCE clustering threshold to analyse (default: first threshold in the database).",
    )
    parser.add_argument(
        "--bigslice-analysis-bins", "--number-of-bins", dest="number_of_bins", type=int, default=520,
        help="Number of K-means GCC bins (default: 520, as in analysis.ipynb).",
    )
    parser.add_argument(
        "--dataset-keyword", default="DeceptionIsland",
        help="Dataset name/keyword for Analysis 1 (default: DeceptionIsland).",
    )
    parser.add_argument(
        "--target-taxon", default="Actinomycetota",
        help="Taxon substring for Analysis 2 (default: Actinomycetota).",
    )
    parser.add_argument(
        "--taxon-rank", default="Phylum",
        help="Taxonomy column to search (default: Phylum).",
    )
    parser.add_argument(
        "--description-keyword", default="marine",
        help="Keyword searched in dataset descriptions for Analysis 3, e.g., 'soil', 'ice', 'ocean', 'gut', 'forest', 'marine' (default: marine).",
    )
    parser.add_argument(
        "--top-n", type=int, default=100,
        help="Maximum number of GCC bins profiled by each exploratory analysis (default: 100).",
    )
    parser.add_argument(
        "--random-seed", type=int, default=17081945,
        help="Random seed for K-means and dataset annotation colors.",
    )
    parser.add_argument("--tree-dpi", type=int, default=600,
                        help="Resolution for the publication PNG (default: 600 dpi).")
    parser.add_argument("--tree-label-font-size", type=float, default=4.0,
                        help="Font size for GCC leaf labels in the tree image (default: 4).")
    parser.add_argument("--hide-tree-labels", action="store_true",
                        help="Hide individual GCC labels for a less crowded tree figure.")
    return parser


def main():

    # Notebook functions and analysis cells are nested in main so their variables
    # share this run's configuration and data.
    import sqlite3
    import pandas as pd
    import numpy as np
    from sklearn.neighbors import NearestNeighbors
    from sklearn.cluster import Birch, KMeans, AgglomerativeClustering
    import pickle
    import random
    from matplotlib import pyplot as plt
    from matplotlib import ticker, rc
    from random import randint
    from scipy.cluster import hierarchy
    from ete3 import PhyloTree
    from os import path


    global args, selected_threshold
    args = parse_args().parse_args()
    project_root = args.project_root.expanduser().resolve()
    if not project_root.is_dir():
        raise FileNotFoundError(f"Project root does not exist: {project_root}")
    os.chdir(project_root)

    db_path = args.db.expanduser().resolve() if args.db else project_root / "source_files/full_run_result/result/data.db"
    if not db_path.is_file():
        raise FileNotFoundError(
            f"BiG-SLiCE database not found: {db_path}\n"
            "Use --db or --project-root to point to the correct location."
        )

    # Keep the paths used by the notebook, but create output directories first.
    Path("source_files").mkdir(parents=True, exist_ok=True)
    Path("tables").mkdir(parents=True, exist_ok=True)
    Path("plots").mkdir(parents=True, exist_ok=True)
    Path("source_files/full_run_result/result/cache").mkdir(parents=True, exist_ok=True)

    # The notebook uses relative database paths in its functions. Change to a
    # project-root-relative symlink-free path only when the default layout is used.
    # If --db points elsewhere, functions below use DB_PATH via a compatibility
    # substitution of the database path string.
    global DB_PATH
    DB_PATH = str(db_path)
    # Set a deterministic random seed for repeatable K-means and color generation.
    random.seed(args.random_seed)
    np.random.seed(args.random_seed)



    # fetch bgc_features
    def fetch_bgc_features():
        cached_bgc_features = "./source_files/cached_bgc_features.pkl"
        if path.exists(cached_bgc_features):
            print("loading bgc features (cached)...")
            bgc_features = pd.read_pickle(cached_bgc_features)
        else:
            print("loading bgc features...")
            with sqlite3.connect(DB_PATH) as con:
                cur = con.cursor()
                bgc_ids = [row[0] for row in cur.execute("select id from bgc order by id asc").fetchall()]
                print(f"Total number of unique BGCs found in database: {len(set(bgc_ids))}")
                hmm_ids = [row[0] for row in cur.execute("select id, name from hmm where db_id=1 order by id asc").fetchall()]
                bgc_features = pd.DataFrame(
                    np.zeros((len(bgc_ids), len(hmm_ids)), dtype=np.float32),
                    index=bgc_ids,
                    columns=hmm_ids
                )
                for bgc_id, hmm_id, value in cur.execute((
                    "select bgc_id, hmm_id, value"
                    " from bgc_features,bgc,hmm"
                    " where bgc_features.bgc_id=bgc.id"
                    " and hmm.id=bgc_features.hmm_id"
                    " and hmm.db_id=1"
                )).fetchall():
                    bgc_features.at[bgc_id, hmm_id] = value
                bgc_features.to_pickle(cached_bgc_features)
            
        # Give the HMM feature columns their names. Query in a separate
        # connection so this also works when the feature matrix was cached.
        with sqlite3.connect(DB_PATH) as con:
            hmm_rows = con.execute(
                "select id, name from hmm where db_id=1 order by id asc"
            ).fetchall()
        hmm_names = pd.DataFrame(hmm_rows, columns=["hmm_id", "name"]).set_index("hmm_id")
        bgc_features = bgc_features.loc[:, hmm_names.index.intersection(bgc_features.columns)]
        bgc_features.columns = hmm_names.loc[bgc_features.columns, "name"].to_numpy()
    
        return bgc_features

    bgc_features = fetch_bgc_features()



    # fetch bgc_taxonomy
    def fetch_bgc_taxonomy():
        cached_bgc_taxonomy = "./source_files/cached_bgc_taxonomy.pkl"
        if path.exists(cached_bgc_taxonomy):
            print("loading bgc taxonomy (cached)...")
            bgc_taxonomy = pd.read_pickle(cached_bgc_taxonomy)
        else:
            print("loading bgc taxonomy...")
            with sqlite3.connect(DB_PATH) as con:
                cur = con.cursor()
                bgc_ids = [row[0] for row in cur.execute("select id from bgc order by id asc").fetchall()]
                taxon_levels = [row[0] for row in cur.execute("select name from taxon_class order by level asc").fetchall()]
                bgc_taxonomy = pd.DataFrame(
                    np.full((len(bgc_ids), len(taxon_levels)), "", dtype=object),
                    index=bgc_ids,
                    columns=taxon_levels
                )
                for bgc_id, taxon_level, name in cur.execute((
                    "select bgc_id,taxon_class.name,taxon.name"
                    " from bgc,bgc_taxonomy,taxon,taxon_class"
                    " where bgc.id=bgc_taxonomy.bgc_id and taxon.id=bgc_taxonomy.taxon_id"
                    " and taxon.level=taxon_class.level"
                )).fetchall():
                    bgc_taxonomy.at[bgc_id, taxon_level] = name
                bgc_taxonomy.to_pickle(cached_bgc_taxonomy)
    
        return bgc_taxonomy

    bgc_taxonomy = fetch_bgc_taxonomy()


    # fetch bgc_names
    def fetch_bgc_meta():
        cached_bgc_meta = "./source_files/cached_bgc_meta.pkl"
        if path.exists(cached_bgc_meta):
            print("loading bgc meta (cached)...")
            bgc_meta = pd.read_pickle(cached_bgc_meta)
        else:
            print("loading bgc meta...")
            with sqlite3.connect(DB_PATH) as con:
                cur = con.cursor()
                bgc_ids, genome_names, bgc_names, contig_edges, lengths  = tuple(zip(*cur.execute(
                    "select id, orig_folder, name, on_contig_edge, length_nt from bgc order by id asc"
                ).fetchall()))
                bgc_meta = pd.DataFrame(
                    {
                        "genome": genome_names,
                        "name": bgc_names,
                        "fragmented": contig_edges,
                        "length": lengths
                    },
                    index=bgc_ids
                )
                #bgc_meta.to_pickle(cached_bgc_meta)
    
        return bgc_meta

    bgc_meta = fetch_bgc_meta()



    import sqlite3
    Available_thresholds_in_DB = []
    print(Available_thresholds_in_DB)
    with sqlite3.connect(DB_PATH) as con:
        cur = con.cursor()
        available_thresholds = cur.execute("select distinct threshold from clustering").fetchall()
        print("Available thresholds in DB:", [t[0] for t in available_thresholds])
        Available_thresholds_in_DB = [t[0] for t in available_thresholds]
        print(Available_thresholds_in_DB)

    if not Available_thresholds_in_DB:
        raise ValueError("No BiG-SLiCE clustering thresholds were found in the database.")
    if args.threshold is None:
        selected_threshold = Available_thresholds_in_DB[0]
    else:
        matching = [t for t in Available_thresholds_in_DB if np.isclose(float(t), args.threshold)]
        if not matching:
            raise ValueError(
                f"Requested threshold {args.threshold} is unavailable. "
                f"Available thresholds: {Available_thresholds_in_DB}"
            )
        selected_threshold = matching[0]
    print(f"Using BiG-SLiCE threshold: {selected_threshold}")


    # fetch gcf centroids from bigslice
    def fetch_gcf_models(threshold):    
        with sqlite3.connect(DB_PATH) as con:
            cur = con.cursor()

            hmm_ids, hmm_names = list(
                # assume the same hmm_db_id for all data
                zip(*cur.execute("select id, name from hmm where db_id=1 order by id asc").fetchall())
            )
        
            clustering_id = cur.execute((
                "select id"
                " from clustering"
                " where threshold=?"
            ), (threshold, )).fetchall()[0][0]

            # check cache
            cache_path = path.join("./source_files/full_run_result/result/cache/", "clustering_{}.pkl".format(clustering_id))
            if path.exists(cache_path):
                print("fetching gcf models for t={} (using cache)".format(threshold))
                df = pd.read_pickle(cache_path).loc[:, hmm_ids]
            else:
                print("fetching gcf models for t={}".format(threshold))
                gcf_ids = [row[0] for row in cur.execute((
                    "select id from gcf"
                    " where clustering_id=?"
                    " order by id"
                ), (clustering_id, )).fetchall()]
                df = pd.DataFrame(np.zeros((len(gcf_ids), len(hmm_ids)), dtype=np.float32), index=gcf_ids, columns=hmm_ids)

                # fetch data
                for gcf_id, hmm_id, value in cur.execute((
                    "select gcf_id, hmm_id, value"
                    " from gcf,gcf_models"
                    " where gcf.id=gcf_models.gcf_id"
                    " and gcf.clustering_id=?"
                ), (clustering_id, )).fetchall():
                    df.at[gcf_id, hmm_id] = value
                
            return df

    gcf_models = fetch_gcf_models(selected_threshold)
                

    # define the number of bins:

    number_of_bins = args.number_of_bins

    print('Number of Bins: ', number_of_bins)

    # generate gcf bin centroids

    def run_kmeans(gcf_models, k):
        km = KMeans(
            n_clusters=k,
            init='random',
            n_init=1,
            max_iter=300,
            algorithm='lloyd',  # 'full' was renamed to 'lloyd'
            random_state=args.random_seed
        )
        if len(gcf_models) < k:
            raise ValueError(
                f"Cannot create {k} GCC bins from only {len(gcf_models)} GCF models. "
                "Reduce --number-of-bins."
            )
        km.fit(gcf_models)
        kmeans_result = {
            "km": km,
            "gcf_labels": pd.DataFrame({"bin": km.labels_}, index=gcf_models.index),
            "bin_centroids": km.cluster_centers_.astype(np.float32)
        }
        return kmeans_result
        
    kmeans_result = run_kmeans(gcf_models, number_of_bins)

    # bgc membership assignment
    def match_bgcs(bgc_features, kmeans_result):
        nn = NearestNeighbors(metric='euclidean', algorithm='brute', n_jobs=1)
        nn.fit(kmeans_result["bin_centroids"])
        return pd.DataFrame({"bin": nn.kneighbors(bgc_features.values, n_neighbors=1)[1][:, 0]}, index=bgc_features.index)

    bgc_bins =  match_bgcs(bgc_features, kmeans_result)


    ### data for Figure 6A ###

    # construct newick tree
    def construct_tree(gcc_centroids, save_to):
        def get_linkage_matrix(model, **kwargs):
            # from https://scikit-learn.org/stable/auto_examples/cluster/plot_agglomerative_dendrogram.html

            # Create linkage matrix and then plot the dendrogram
            # create the counts of samples under each node
            counts = np.zeros(model.children_.shape[0])
            n_samples = len(model.labels_)
            for i, merge in enumerate(model.children_):
                current_count = 0
                for child_idx in merge:
                    if child_idx < n_samples:
                        current_count += 1  # leaf node
                    else:
                        current_count += counts[child_idx - n_samples]
                counts[i] = current_count
            linkage_matrix = np.column_stack([model.children_, model.distances_,
                                              counts]).astype(float)
            return linkage_matrix

        def getNewick(node, newick, parentdist, leaf_names):
            # from https://stackoverflow.com/questions/28222179/save-dendrogram-to-newick-format
            if node.is_leaf():
                return "%s:%.2f%s" % (leaf_names[node.id], parentdist - node.dist, newick)
            else:
                if len(newick) > 0:
                    newick = "):%.2f%s" % (parentdist - node.dist, newick)
                else:
                    newick = ");"
                newick = getNewick(node.get_left(), newick, node.dist, leaf_names)
                newick = getNewick(node.get_right(), ",%s" % (newick), node.dist, leaf_names)
                newick = "(%s" % (newick)
                return newick

        # append a zero GCC feature for tree rooting
        _gcc_centroids = np.append(
            gcc_centroids,
            np.zeros((1, gcc_centroids.shape[1]), dtype=np.float32), axis=0
        )

        # perform hierarchical clustering
        clusterer = AgglomerativeClustering(
            n_clusters=None,
            distance_threshold=0,
            metric="euclidean",  # Changed from 'affinity' to 'metric'
            linkage="ward"
        )
        clusterer.fit(_gcc_centroids)    

        # generate unrooted tree newick
        scipy_tree = hierarchy.to_tree(get_linkage_matrix(clusterer), False)
        ete_tree = PhyloTree(getNewick(scipy_tree, "", scipy_tree.dist, [
            "GCC{:05d}".format(gcc_id) for gcc_id in np.arange(scipy_tree.count)
        ]))

        # reroot tree to the dummy zeros feature
        ete_tree.set_outgroup("GCC{:05d}".format(scipy_tree.count - 1))

        # write newick to file
        with open(save_to, "w") as txt:
            txt.write(ete_tree.write())
        
        return ete_tree
        
    ete_tree = construct_tree(kmeans_result["bin_centroids"], "./plots/figure_6_tree.newick")

    # construct metadata
    def get_phylo_metadata(threshold):
        bin_members = bgc_bins.groupby("bin").indices
        bin_counts = bgc_bins.groupby(["bin"]).size()
        phylo_metadata = pd.DataFrame(bin_counts, columns=["bgc_count"])
    
        # annotation 1: average "signal strength" (feature values)
        print("annotation 1: average \"signal strength\" (feature values)")
        def calc_signal_strength(bin_id):
            return int(bgc_features.iloc[bin_members[bin_id]].sum(axis=1).mean())

        phylo_metadata["avg_signal_strength"] = phylo_metadata.index.map(lambda bin_id: calc_signal_strength(bin_id))
    
        # annotation 2: average bgc lengths
        print("annotation 2: average bgc lengths")
        with sqlite3.connect(DB_PATH) as con:
            cur = con.cursor()
            bgc_lengths = pd.DataFrame.from_dict({row[0]: row[1] for row in cur.execute((
                "select id, length_nt"
                " from bgc order by id asc"
            )).fetchall()}, orient="index", columns=["length"])

        def calc_bgc_length(bin_id):
            return int(bgc_lengths.iloc[bin_members[bin_id]]["length"].mean())

        phylo_metadata["avg_bgc_length"] = phylo_metadata.index.map(lambda bin_id: calc_bgc_length(bin_id))
    
        # annotation 3: average number of BGCs in a GCF
        print("annotation 3: average number of BGCs in a GCF")
        bin_gcfs = kmeans_result["gcf_labels"].groupby("bin").groups
        clustering_id = cur.execute((
            "select id"
            " from clustering"
            " where threshold=?"
        ), (threshold, )).fetchall()[0][0]
        gcf_sizes = pd.DataFrame.from_dict({row[0]: row[1] for row in cur.execute((
            "select gcf_id, count(bgc_id)"
            " from gcf,gcf_membership"
            " where gcf.id=gcf_membership.gcf_id"
            " and gcf.clustering_id=?"
            " and gcf_membership.rank=0"
            " and gcf_membership.membership_value <= ?"
            " group by gcf_id"
        ), (clustering_id, 2*threshold)).fetchall()}, orient="index", columns=["size"])

        def calc_gcf_size(bin_id):
            return int(gcf_sizes.loc[bin_gcfs[bin_id]]["size"].mean())

        phylo_metadata["avg_gcf_size"] = phylo_metadata.index.map(lambda bin_id: calc_gcf_size(bin_id))
    
        # annotation 4: bgc generic classes distribution
        print("annotation 4: bgc generic classes distribution")
        with sqlite3.connect(DB_PATH) as con:
            cur = con.cursor()
            class_ids, class_names = tuple(zip(*cur.execute((
                "select id, name from chem_class where name not like 'Unknown'"
            )).fetchall()))
            bgc_classes = pd.DataFrame({class_id: [False]*len(bgc_features) for class_id in class_ids}, index=bgc_features.index)
        
            for bgc_id, class_id in cur.execute((
                "select bgc.id, chem_subclass.class_id"
                " from bgc,bgc_class,chem_subclass"
                " where bgc.id=bgc_class.bgc_id"
                " and bgc_class.chem_subclass_id=chem_subclass.id"
            )).fetchall():
                if bgc_id in bgc_classes.index: # SAFETY CHECK
                    bgc_classes.at[bgc_id, class_id] = True

        def calc_bgc_class(bin_id, class_id):
            # Look up exact BGC IDs instead of blind row positions
            bin_bgc_ids = bgc_bins[bgc_bins["bin"] == bin_id].index
            return int(bgc_classes.loc[bin_bgc_ids, class_id].sum())

        # --- BATCH ADDITION ---
        class_cols = {}
        for i, class_id in enumerate(class_ids):
            class_cols["class_" + class_names[i]] = phylo_metadata.index.map(lambda bin_id: calc_bgc_class(bin_id, class_id))
        phylo_metadata = pd.concat([phylo_metadata, pd.DataFrame(class_cols, index=phylo_metadata.index)], axis=1)
    
        # annotation 5: bgc dataset distribution
        print("annotation 5: bgc dataset distribution")
        with sqlite3.connect(DB_PATH) as con:
            cur = con.cursor()
            dataset_ids, dataset_names = tuple(zip(*cur.execute((
                "select id, name from dataset order by name asc"
            )).fetchall()))
            bgc_datasets = pd.DataFrame.from_dict({row[0]: row[1] for row in cur.execute((
                "select id, dataset_id"
                " from bgc order by id asc"
            )).fetchall()}, orient="index", columns=["dataset_id"])

        def calc_bgc_dataset(bin_id, dataset_id):
            # Look up exact BGC IDs
            bin_bgc_ids = bgc_bins[bgc_bins["bin"] == bin_id].index
            return int((bgc_datasets.loc[bin_bgc_ids, "dataset_id"] == dataset_id).sum())

        # --- BATCH ADDITION ---
        dataset_cols = {}
        for i, dataset_id in enumerate(dataset_ids):
            dataset_cols["dataset_" + dataset_names[i]] = phylo_metadata.index.map(lambda bin_id: calc_bgc_dataset(bin_id, dataset_id))
        phylo_metadata = pd.concat([phylo_metadata, pd.DataFrame(dataset_cols, index=phylo_metadata.index)], axis=1)

        # extra annotation: BGC subclasses
        print("extra annotation: BGC subclasses")
        with sqlite3.connect(DB_PATH) as con:
            cur = con.cursor()
            subclass_ids, subclass_names = tuple(zip(*cur.execute((
                "select id, name from chem_subclass"
            )).fetchall()))
            bgc_subclasses = pd.DataFrame({subclass_id: [False]*len(bgc_features) for subclass_id in subclass_ids}, index=bgc_features.index)
        
            for bgc_id, subclass_id in cur.execute((
                "select bgc_id, bgc_class.chem_subclass_id"
                " from bgc_class"
            )).fetchall():
                if bgc_id in bgc_subclasses.index: # CRITICAL SAFETY CHECK
                    bgc_subclasses.at[bgc_id, subclass_id] = True

        def calc_bgc_subclass(bin_id, subclass_id):
            # Look up exact BGC IDs
            bin_bgc_ids = bgc_bins[bgc_bins["bin"] == bin_id].index
            return int(bgc_subclasses.loc[bin_bgc_ids, subclass_id].sum())

        # --- BATCH ADDITION WITH OVERWRITE PROTECTION ---
        subclass_cols = {}
        for i, subclass_id in enumerate(subclass_ids):
            col_name = "subclass_" + subclass_names[i]
            counts = phylo_metadata.index.map(lambda bin_id: calc_bgc_subclass(bin_id, subclass_id))
        
            # PREVENT OVERWRITING: If "subclass_generic" already exists, add the new counts to the old ones!
            if col_name in subclass_cols:
                subclass_cols[col_name] += counts
            else:
                subclass_cols[col_name] = counts
            
        phylo_metadata = pd.concat([phylo_metadata, pd.DataFrame(subclass_cols, index=phylo_metadata.index)], axis=1)
    
        return phylo_metadata


    phylo_metadata = get_phylo_metadata(selected_threshold)

    phylo_metadata.to_csv(
        "./tables/supplementary_table_5.tsv",
        sep="\t",
        index_label="bin"
    )

    print("Rendering publication-style GCC tree...")
    render_publication_tree(
        ete_tree,
        phylo_metadata,
        Path("./plots"),
        dpi=args.tree_dpi,
        label_font_size=args.tree_label_font_size,
        show_labels=not args.hide_tree_labels,
    )
    print("Rendered tree: ./plots/figure_6_tree.png, .pdf, .svg")


    ## save annotation tables for iTOL ##

    # annotation 1: average "signal strength"
    def annotate_signal_strength(save_to):
        annotation_txt = """
    DATASET_GRADIENT
    SEPARATOR SPACE
    DATASET_LABEL signal_strength
    COLOR #ff0000

    #define the gradient colors. Values in the dataset will be mapped onto the corresponding color gradient.
    COLOR_MIN #ffffff
    COLOR_MAX #000000

    DATA
    """
        # --- CHANGED iteritems() TO items() ---
        for idx, value in phylo_metadata["avg_signal_strength"].items(): 
            bin_name = "GCC{:05d}".format(idx)
            annotation_txt += "{} {}\n".format(bin_name, value)
    
        with open(save_to, "w") as fo:
            fo.write(annotation_txt)
    
    annotate_signal_strength("./plots/figure_6-annot-1.txt")

    # annotation 2: average bgc length
    def annotate_bgc_length(save_to):
        annotation_txt = """
    DATASET_GRADIENT
    SEPARATOR SPACE
    DATASET_LABEL bgc_length
    COLOR #000000

    #define the gradient colors. Values in the dataset will be mapped onto the corresponding color gradient.
    COLOR_MIN #ffffff
    COLOR_MAX #000000

    DATA
    """
        for idx, value in phylo_metadata["avg_bgc_length"].items(): 
            bin_name = "GCC{:05d}".format(idx)
            annotation_txt += "{} {}\n".format(bin_name, value)
    
        with open(save_to, "w") as fo:
            fo.write(annotation_txt)
    
    annotate_bgc_length("./plots/figure_6-annot-2.txt")

    # annotation 3: average bgc counts
    def annotate_bgc_counts(save_to):
        annotation_txt = """
    DATASET_GRADIENT
    SEPARATOR SPACE
    DATASET_LABEL bgc_counts
    COLOR #0000ff

    #define the gradient colors. Values in the dataset will be mapped onto the corresponding color gradient.
    COLOR_MIN #ffffff
    COLOR_MAX #0000ff

    DATA
    """
        for idx, value in phylo_metadata["avg_gcf_size"].map(lambda x: min(x, 50)).items(): 
            bin_name = "GCC{:05d}".format(idx)
            annotation_txt += "{} {}\n".format(bin_name, value)
    
        with open(save_to, "w") as fo:
            fo.write(annotation_txt)
    
    annotate_bgc_counts("./plots/figure_6-annot-3.txt")

    # annotation 4: bgc generic classes distribution
    def annotate_chem_class(save_to):
        all_classes = [
            "Polyketide",
            "NRP",
            "RiPP",
            "Saccharide",
            "Terpene",
            "Alkaloid",
            "Other"
        ]
        all_classes_shape = {
            "Polyketide": ("rgba(255,135,43,0.5)", 1),
            "NRP": ("rgba(28,178,38,0.5)", 1),
            "RiPP": ("rgba(214,72,89,0.5)", 1),
            "Saccharide": ("rgba(1,1,173,0.5)", 1),
            "Terpene": ("rgba(174,0,174,0.5)", 1),
            "Alkaloid": ("rgba(255,255,56,0.5)", 1),
            "Other": ("rgba(126,126,127,0.5)", 1),
        }

        itol_annot = """
    DATASET_EXTERNALSHAPE
    SEPARATOR SPACE

    DATASET_LABEL class_ratio
    COLOR #ffffff

    FIELD_COLORS """ + " ".join([all_classes_shape[chem_class][0] for chem_class in all_classes]) + """
    FIELD_LABELS """ + " ".join(all_classes) + """

    DATA
        """

        col_idx = ["class_" + name for name in all_classes]
        for idx, rows in phylo_metadata.iterrows():
            row_classes = rows[col_idx]
            total = row_classes.sum()
            itol_annot += "{} {}\n".format(
                "GCC{:05d}".format(idx),
                " ".join(map("{:.4f}".format, (row_classes.values) / total))
            )

        with open(save_to, "w") as txt:
            txt.write(itol_annot)
        
    annotate_chem_class("./plots/figure_6-annot-4.txt")

    import sqlite3
    import random

    # annotation 5: bgc dataset distribution
    def annotate_dataset_counts(save_to):
        # Automatically fetch datasets from the database
        with sqlite3.connect(DB_PATH) as con:
            cur = con.cursor()
            all_datasets = [row[0] for row in cur.execute("select name from dataset order by name asc").fetchall()]
        
        # Dynamically assign a random color to each dataset
        all_datasets_color = {}
        for dataset in all_datasets:
            r = random.randint(0, 200) # Kept under 200 to avoid colors that are too bright/white
            g = random.randint(0, 200)
            b = random.randint(0, 200)
            all_datasets_color[dataset] = "rgba({}, {}, {}, 0.7)".format(r, g, b)
    
        itol_annot = """
    DATASET_MULTIBAR
    SEPARATOR SPACE

    DATASET_LABEL dataset_distribution
    COLOR #ffffff

    FIELD_COLORS """ + " ".join([all_datasets_color[dataset] for dataset in all_datasets]) + """
    FIELD_LABELS """ + " ".join(all_datasets) + """

    DATA
        """

        col_idx = ["dataset_" + name for name in all_datasets]
        sum_datasets = phylo_metadata.loc[:, col_idx].sum()
    
        for idx, rows in phylo_metadata.iterrows():
            row_datasets = rows[col_idx]
            ratio_within_datasets = (row_datasets.values / sum_datasets)
            total_ratio = ratio_within_datasets.sum()
        
            # Handle cases where total_ratio is 0 to avoid division by zero errors
            if total_ratio > 0:
                final_ratios = ratio_within_datasets / total_ratio
            else:
                final_ratios = ratio_within_datasets
            
            itol_annot += "{} {}\n".format(
                "GCC{:05d}".format(idx),
                " ".join(map("{:.4f}".format, final_ratios))
            )

        with open(save_to, "w") as txt:
            txt.write(itol_annot)
        
    annotate_dataset_counts("./plots/figure_6-annot-5.txt")


    # Print all column names that are datasets
    print([col for col in phylo_metadata.columns if col.startswith('dataset_')])

    # ==========================================
    # ANALYSIS 1: BY ENVIRONMENT / DATASET
    # ==========================================

    # 1. Define your environment keyword (e.g., 'DeceptionIsland', 'TARA', 'WhaleFall', 'BlackSea', or 'LlaimaVolcano')
    dataset_keyword = ('dataset_' + args.dataset_keyword)
    print('Searching for the dataset:', dataset_keyword)

    # 2. Find which columns represent these samples
    env_columns = [col for col in phylo_metadata.columns if dataset_keyword.lower() in col.lower()]
    print(env_columns)

    if not env_columns:
        print(f"Warning: No datasets with '{dataset_keyword}' in the name were found.")
    else:
        # 3. Sort bins by how many hits they contain for this environment
        env_bins_counts = phylo_metadata[env_columns].sum(axis=1).sort_values(ascending=False)
    
        # Grab the top 100 bins with the most hits
        env_bins_counts = env_bins_counts[env_bins_counts > 0]
        env_bin_idxs = env_bins_counts.head(args.top_n).index.tolist()
    
        print(f"--- Top {args.top_n} Bins for Environment: '{dataset_keyword}' ---")
        print(f"Bin IDs: {env_bin_idxs}\n")

        if env_bin_idxs:
            # 4. Run subclass (chemistry) analysis on these bins
            env_subclasses = phylo_metadata.loc[env_bin_idxs, phylo_metadata.columns.str.startswith("subclass_")].sum().sort_values(ascending=False)
            env_subclasses = env_subclasses[env_subclasses > 0]
        
            print(f"Subclass distribution (%) for '{dataset_keyword}' bins:")
            print(round((env_subclasses / env_subclasses.sum()) * 100, 1))

            # 5. Run broad class (chemistry) analysis on these bins
            env_classes = phylo_metadata.loc[env_bin_idxs, phylo_metadata.columns.str.startswith("class_")].sum().sort_values(ascending=False)
            env_classes = env_classes[env_classes > 0]
        
            print(f"\nBroad Class distribution (%) for '{dataset_keyword}' bins:")
            print(round((env_classes / env_classes.sum()) * 100, 1))




    # ==========================================
    # ANALYSIS 2: BY TAXONOMY / STRAIN
    # ==========================================

    # 1. Define your target bacteria (e.g., 'Psychrobacter', 'Steptomyces'), or phyla (change code below in 'Genus')
    target_genus = args.target_taxon

    bgc_taxonomy.to_csv(
        "./tables/taxonomy.tsv",
        sep="\t",
        index_label="bin"
    )

    # 2. Find all BGCs in the taxonomy table that match this genus
    if args.taxon_rank not in bgc_taxonomy.columns:
        raise ValueError(f"Taxonomy rank '{args.taxon_rank}' not found. Available ranks: {list(bgc_taxonomy.columns)}")
    target_bgcs = bgc_taxonomy[bgc_taxonomy[args.taxon_rank].str.contains(target_genus, case=False, na=False)]

    if target_bgcs.empty:
        print(f"Warning: No BGCs found for genus '{target_genus}'.")
    else:
        print(f"--- Profiling Taxonomy: '{target_genus}' ---")
        print(f"Found {len(target_bgcs)} total BGCs assigned to this genus.\n")

        # 3. Find which bins contain these specific BGCs
        target_bin_assignments = bgc_bins.loc[target_bgcs.index]
        bin_counts = target_bin_assignments['bin'].value_counts()
    
        # Select the top 4 bins where this genus clusters
        taxa_bin_idxs = bin_counts.head(args.top_n).index.tolist()
        print(f"Top 100 Bin IDs where {target_genus} clusters: {taxa_bin_idxs}\n")

        # 4. Run subclass (chemistry) analysis on these bins
        taxa_subclasses = phylo_metadata.loc[taxa_bin_idxs, phylo_metadata.columns.str.startswith("subclass_")].sum().sort_values(ascending=False)
        taxa_subclasses = taxa_subclasses[taxa_subclasses > 0]

        print(f"Subclass distribution (%) for '{target_genus}' bins:")
        print(round((taxa_subclasses / taxa_subclasses.sum()) * 100, 1))

        print(" ")
            
        # 5. Run broad class (chemistry) analysis on these bins
        taxa_classes = phylo_metadata.loc[taxa_bin_idxs, phylo_metadata.columns.str.startswith("class_")].sum().sort_values(ascending=False)
        taxa_classes = taxa_classes[taxa_classes > 0]
        print(f"Class distribution (%) for '{target_genus}' environment bins:")
        print(round((taxa_classes / taxa_classes.sum()) * 100, 1))



    
        print("\n---------------------------------------------------")
    
        # 5. Check the dataset origins for these specific bins
        all_dataset_cols = phylo_metadata.columns[phylo_metadata.columns.str.startswith("dataset_")]
        taxa_datasets = phylo_metadata.loc[taxa_bin_idxs, all_dataset_cols].sum()
        taxa_datasets = taxa_datasets[taxa_datasets > 0].sort_values(ascending=False)
    
        print(f"Dataset origins for these '{target_genus}' bins:")
        for dataset_name, count in taxa_datasets.items():
            clean_name = dataset_name.replace("dataset_", "")
            print(f" - {clean_name}: {int(count)} BGCs")

    # ==========================================
    # ANALYSIS 3: BY DATASET DESCRIPTION (METADATA)
    # ==========================================
    import sqlite3

    # 1. Define your ecological keyword to search in descriptions (e.g., 'soil', 'ice', 'ocean', 'gut', 'forest', 'marine')
    description_keyword = args.description_keyword

    # 2. Query the database to find datasets with this keyword in their description
    with sqlite3.connect(DB_PATH) as con:
        cur = con.cursor()
        # Using LOWER() to make the SQL search case-insensitive
        query = "SELECT name, description FROM dataset WHERE LOWER(description) LIKE ?"
        matches = cur.execute(query, ('%' + description_keyword.lower() + '%',)).fetchall()

    if not matches:
        print(f"Warning: No datasets found with '{description_keyword}' in their description.")
    else:
        print(f"Found {len(matches)} dataset(s) matching '{description_keyword}' in the description.")
        matched_dataset_names = [row[0] for row in matches]
    
        # 3. Find the corresponding columns in your metadata table
        desc_columns = [f"dataset_{name}" for name in matched_dataset_names if f"dataset_{name}" in phylo_metadata.columns]
    
        if not desc_columns:
            print("Warning: Matched datasets have no corresponding BGCs at the current clustering threshold.")
        else:
            # 4. Sort bins by how many hits they contain for these environment-specific datasets
            desc_bins_counts = phylo_metadata[desc_columns].sum(axis=1).sort_values(ascending=False)
        
            # Grab the top 4 bins with the most hits
            desc_bins_counts = desc_bins_counts[desc_bins_counts > 0]
            desc_bin_idxs = desc_bins_counts.head(args.top_n).index.tolist()
        
            print(f"\n--- Top {args.top_n} Bins for Environment Description: '{description_keyword}' ---")
            print(f"Bin IDs: {desc_bin_idxs}\n")

            if desc_bin_idxs:
                # 5. Run subclass (chemistry) analysis on these specific bins
                desc_subclasses = phylo_metadata.loc[desc_bin_idxs, phylo_metadata.columns.str.startswith("subclass_")].sum().sort_values(ascending=False)
                desc_subclasses = desc_subclasses[desc_subclasses > 0]
            
                print(f"Subclass distribution (%) for '{description_keyword}' environment bins:")
                print(round((desc_subclasses / desc_subclasses.sum()) * 100, 1))
                       
                print(" ")
            
                # 5. Run broad class (chemistry) analysis on these bins
                desc_classes = phylo_metadata.loc[desc_bin_idxs, phylo_metadata.columns.str.startswith("class_")].sum().sort_values(ascending=False)
                desc_classes = desc_classes[desc_classes > 0]
                print(f"Class distribution (%) for '{description_keyword}' environment bins:")
                print(round((desc_classes / desc_classes.sum()) * 100, 1))
            
            
            
                print("\n---------------------------------------------------")
                print("Specific contributing datasets found in these top bins:")
            
                # Show exactly which datasets from your SQL query are showing up in these top bins
                desc_datasets_present = phylo_metadata.loc[desc_bin_idxs, desc_columns].sum()
                desc_datasets_present = desc_datasets_present[desc_datasets_present > 0].sort_values(ascending=False)
            
                for dataset_col, count in desc_datasets_present.items():
                    clean_name = dataset_col.replace("dataset_", "")
                    print(f" - {clean_name}: {int(count)} BGCs")




    # ==========================================
    # CONSOLIDATED DATA SLICER
    # ==========================================
    sliced_data = {}

    # 1. Slice for Environment Analysis
    if 'env_bin_idxs' in locals() and env_bin_idxs:
        env_bgcs = bgc_bins[bgc_bins["bin"].isin(env_bin_idxs)]
        sliced_data["Environment"] = {
            "taxonomy": bgc_taxonomy.loc[env_bgcs.index],
            "features": bgc_features.loc[env_bgcs.index]
        }

    # 2. Slice for Taxonomy Analysis
    if 'taxa_bin_idxs' in locals() and taxa_bin_idxs:
        taxa_bgcs = bgc_bins[bgc_bins["bin"].isin(taxa_bin_idxs)]
        sliced_data["Taxonomy"] = {
            "taxonomy": bgc_taxonomy.loc[taxa_bgcs.index],
            "features": bgc_features.loc[taxa_bgcs.index]
        }

    # 3. Slice for Description/Metadata Analysis
    if 'desc_bin_idxs' in locals() and desc_bin_idxs:
        desc_bgcs = bgc_bins[bgc_bins["bin"].isin(desc_bin_idxs)]
        sliced_data["Description"] = {
            "taxonomy": bgc_taxonomy.loc[desc_bgcs.index],
            "features": bgc_features.loc[desc_bgcs.index]
        }

    print(f"Successfully prepared data slices for: {list(sliced_data.keys())}")

    # ==========================================
    # UNIVERSAL CLADE PROFILER
    # ==========================================
    def profile_target_bins(analysis_name, subset_taxonomy, subset_features):
        print(f"\n========================================================")
        print(f" PROFILING RESULTS FOR: {analysis_name.upper()}")
        print(f"========================================================")
    
        # 1. Kingdom Distribution
        print("\n--- 1. Kingdom Distribution ---")
        print(subset_taxonomy.groupby(["Kingdom"])["Organism"].count().sort_values(ascending=False))
    
        # 2. Bacterial Taxonomy Breakdown
        bac_taxa = subset_taxonomy[subset_taxonomy["Kingdom"] == "Bacteria"]
        if not bac_taxa.empty:
            print("\n--- 2. Top 5 Bacterial Phyla ---")
            print(bac_taxa.groupby(["Phylum"])["Organism"].count().sort_values(ascending=False).head(5))
        
            print("\n--- 3. Top 5 Bacterial Classes ---")
            print(bac_taxa.groupby(["Class"])["Organism"].count().sort_values(ascending=False).head(5))
        
            print("\n--- 4. Top 10 Bacterial Genera ---")
            print(bac_taxa.groupby(["Genus"])["Organism"].count().sort_values(ascending=False).head(10))
        else:
            print("\n(No Bacterial BGCs found in this slice to profile taxonomy.)")

        # 3. Enriched Biosynthetic Features (Top HMMs)
        print("\n--- 5. Top 15 Enriched HMM Features (Average Score) ---")
        # Calculate the mean score of every feature across these specific BGCs
        mean_features = subset_features.mean().sort_values(ascending=False)
        # Filter out zeros to only show features that are actually present
        top_features = mean_features[mean_features > 0].head(15)
        print(top_features)
        print("========================================================\n")

    # Execute the profiler for all active analyses
    if not sliced_data:
        print("No sliced data available. Please run Analysis 1, 2, or 3 first.")
    else:
        for name, data in sliced_data.items():
            profile_target_bins(name, data["taxonomy"], data["features"])

if __name__ == "__main__":
    try:
        main()
    except (FileNotFoundError, ValueError, RuntimeError, sqlite3.Error) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
