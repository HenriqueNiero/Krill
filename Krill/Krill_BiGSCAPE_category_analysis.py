#!/usr/bin/env python3

import sqlite3
import pandas as pd
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.manifold import MDS
import networkx as nx
from matplotlib.patches import Polygon
from scipy.spatial import ConvexHull


# ============================================================
# Configuration
# ============================================================

MIBIG_IDENTIFIER = "/site-packages/big_scape/MIBiG"




# ============================================================
# Utility functions
# ============================================================

def is_mibig(path):
    """
    Return True if a path belongs to the BiG-SCAPE MIBiG database.
    """
    if pd.isna(path):
        return False

    return MIBIG_IDENTIFIER in str(path)


def is_user_data(path, input_path):
    """
    Return True if the path belongs to the user's Krill input directory.
    """
    if pd.isna(path):
        return False

    path = os.path.abspath(str(path))
    input_path = os.path.abspath(input_path)

    return path.startswith(input_path)


def get_database_from_path(path, input_path):
    """
    Determine whether a BGC belongs to MIBiG or to the user dataset.

    For user data, the first directory below input_path is used as the
    database/sample name when possible.
    """
    if pd.isna(path):
        return "Unknown"

    path = os.path.abspath(str(path))
    input_path = os.path.abspath(input_path)

    if is_mibig(path):
        return "MIBiG"

    if path.startswith(input_path):
        relative = os.path.relpath(path, input_path)

        parts = relative.split(os.sep)

        if len(parts) > 1:
            return parts[0]

        return os.path.basename(path)

    return "Unknown"


# ============================================================
# Read SQLite tables
# ============================================================

def read_tables(db_file):

    print("\n==============================================")
    print("Reading BiG-SCAPE SQLite database")
    print("==============================================")

    if not os.path.isfile(db_file):
        raise FileNotFoundError(
            f"BiG-SCAPE database not found:\n{db_file}"
        )

    conn = sqlite3.connect(db_file)

    distance = pd.read_sql_query(
        "SELECT * FROM distance",
        conn
    )

    gbk = pd.read_sql_query(
        "SELECT * FROM gbk",
        conn
    )

    bgc_record = pd.read_sql_query(
        "SELECT * FROM bgc_record",
        conn
    )

    conn.close()

    print(f"distance rows:    {len(distance)}")
    print(f"gbk rows:         {len(gbk)}")
    print(f"bgc_record rows:  {len(bgc_record)}")

    return distance, gbk, bgc_record


# ============================================================
# Add GBK path to bgc_record
# ============================================================

def add_path_to_bgc_record(bgc_record, gbk):

    print("\n==============================================")
    print("Adding GBK path to bgc_record")
    print("==============================================")

    gbk_path = gbk[
        [
            "id",
            "path"
        ]
    ].copy()

    gbk_path = gbk_path.rename(
        columns={
            "id": "gbk_id",
            "path": "path"
        }
    )

    # Many bgc_record rows can have the same gbk_id.
    # This is expected.
    bgc_record = bgc_record.merge(
        gbk_path,
        on="gbk_id",
        how="left"
    )

    missing = bgc_record["path"].isna().sum()

    print(f"bgc_record rows: {len(bgc_record)}")
    print(f"Missing paths:   {missing}")

    return bgc_record


# ============================================================
# Add BGC information to distance table
# ============================================================

def add_bgc_information_to_distance(
    distance,
    bgc_record
):

    print("\n==============================================")
    print("Adding BGC information to distance")
    print("==============================================")

    lookup = bgc_record[
        [
            "id",
            "path",
            "record_type",
            "gbk_id",
            "record_number",
            "product",
            "category"
        ]
    ].copy()

    # --------------------------------------------------------
    # record_a
    # --------------------------------------------------------

    lookup_a = lookup.rename(
        columns={
            "id": "record_a_id",
            "path": "record_a_id_path",
            "record_type": "record_a_record_type",
            "gbk_id": "record_a_gbk_id",
            "record_number": "record_a_record_number",
            "product": "record_a_product",
            "category": "record_a_category"
        }
    )

    distance = distance.merge(
        lookup_a,
        on="record_a_id",
        how="left"
    )

    # --------------------------------------------------------
    # record_b
    # --------------------------------------------------------

    lookup_b = lookup.rename(
        columns={
            "id": "record_b_id",
            "path": "record_b_id_path",
            "record_type": "record_b_record_type",
            "gbk_id": "record_b_gbk_id",
            "record_number": "record_b_record_number",
            "product": "record_b_product",
            "category": "record_b_category"
        }
    )

    distance = distance.merge(
        lookup_b,
        on="record_b_id",
        how="left"
    )

    print("Added BGC metadata to distance table.")

    return distance


# ============================================================
# Add distance source
# ============================================================

def add_distance_source(distance):

    """
    Label original BiG-SCAPE distances.

    No distance values are generated or modified.
    """

    distance = distance.copy()

    distance["distance_source"] = "BiG-SCAPE"

    return distance


# ============================================================
# Identify user cand_cluster BGCs
# ============================================================

def get_user_candidate_clusters(
    bgc_record,
    input_path
):

    print("\n==============================================")
    print("Finding user cand_cluster BGCs")
    print("==============================================")

    user_mask = bgc_record[
        "path"
    ].apply(
        lambda x: is_user_data(x, input_path)
    )

    candidate_mask = (
        bgc_record["record_type"]
        == "cand_cluster"
    )

    user_candidates = bgc_record[
        user_mask & candidate_mask
    ].copy()

    print(
        f"User cand_cluster BGCs: "
        f"{len(user_candidates)}"
    )

    return user_candidates


# ============================================================
# Identify singleton BGCs
# ============================================================

def identify_singletons(
    distance,
    user_candidates,
    threshold=0.3
):
    """
    Classify user cand_cluster BGCs according to their
    minimum BiG-SCAPE distance.

    For every user BGC:

        min_distance =
            minimum distance to any OTHER cand_cluster

    Singleton definition:

        min_distance > threshold

    Therefore, with threshold = 0.3:

        distance <= 0.3  -> non-singleton
        distance >  0.3  -> singleton

    Important:
    - MIBiG comparisons are included when calculating min_distance.
    - No distance values are created or modified.
    - Missing user-user pairs are NOT interpreted as distance=1.
    """

    print("\n==============================================")
    print(
        f"Identifying singletons using minimum "
        f"distance and threshold <= {threshold}"
    )
    print("==============================================")

    # --------------------------------------------------------
    # User BGC IDs
    # --------------------------------------------------------

    user_ids = set(
        pd.to_numeric(
            user_candidates["id"],
            errors="coerce"
        )
        .dropna()
        .astype(int)
    )

    # --------------------------------------------------------
    # Prepare distance table
    # --------------------------------------------------------

    d = distance.copy()

    d["record_a_id"] = pd.to_numeric(
        d["record_a_id"],
        errors="coerce"
    )

    d["record_b_id"] = pd.to_numeric(
        d["record_b_id"],
        errors="coerce"
    )

    d["distance"] = pd.to_numeric(
        d["distance"],
        errors="coerce"
    )

    d = d.dropna(
        subset=[
            "record_a_id",
            "record_b_id",
            "distance"
        ]
    )

    # Remove self-comparisons
    d = d[
        d["record_a_id"] != d["record_b_id"]
    ]

    # Only cand_cluster × cand_cluster
    d = d[
        (d["record_a_record_type"] == "cand_cluster") &
        (d["record_b_record_type"] == "cand_cluster")
    ]

    # --------------------------------------------------------
    # Store every distance associated with each BGC
    # --------------------------------------------------------

    bgc_distances = {
        bgc_id: []
        for bgc_id in user_ids
    }

    for row in d.itertuples(index=False):

        a = int(row.record_a_id)
        b = int(row.record_b_id)

        value = float(row.distance)

        # If A is one of the user's BGCs,
        # this distance belongs to A.
        if a in bgc_distances:
            bgc_distances[a].append(value)

        # If B is one of the user's BGCs,
        # this distance belongs to B.
        if b in bgc_distances:
            bgc_distances[b].append(value)

    # --------------------------------------------------------
    # Build status table
    # --------------------------------------------------------

    status = user_candidates.copy()

    # --------------------------------------------------------
    # Number of comparisons
    # --------------------------------------------------------

    status["n_comparisons"] = (
        status["id"]
        .apply(
            lambda x:
            len(
                bgc_distances.get(
                    int(x),
                    []
                )
            )
        )
    )

    # --------------------------------------------------------
    # Minimum distance
    #
    # This is the minimum distance across ALL observed
    # cand_cluster comparisons for that BGC.
    # --------------------------------------------------------

    def get_min_distance(bgc_id):

        distances = bgc_distances.get(
            int(bgc_id),
            []
        )

        if len(distances) == 0:
            return np.nan

        return min(distances)

    status["min_distance"] = (
        status["id"]
        .apply(get_min_distance)
    )

    # --------------------------------------------------------
    # Number of neighbors <= threshold
    # --------------------------------------------------------

    status["n_neighbors_at_threshold"] = (
        status["id"]
        .apply(
            lambda x:
            sum(
                d <= threshold
                for d in bgc_distances.get(
                    int(x),
                    []
                )
            )
        )
    )

    # --------------------------------------------------------
    # Singleton classification
    #
    # A BGC is a singleton if its minimum distance is
    # greater than the threshold.
    #
    # If min_distance is missing because there are literally
    # no comparisons, we do NOT silently call it a singleton.
    # --------------------------------------------------------

    status["threshold"] = threshold

    status["singleton"] = (
        status["min_distance"].notna()
        &
        (
            status["min_distance"]
            > threshold
        )
    )

    status["status"] = np.where(
        status["singleton"],
        "singleton",
        "non-singleton"
    )

    # --------------------------------------------------------
    # Report
    # --------------------------------------------------------

    singleton_count = int(
        status["singleton"].sum()
    )

    non_singleton_count = int(
        (~status["singleton"]).sum()
    )

    no_comparison_count = int(
        status["min_distance"].isna().sum()
    )

    print(
        f"Singletons:              {singleton_count}"
    )

    print(
        f"Non-singletons:          {non_singleton_count}"
    )

    print(
        f"No comparisons:          {no_comparison_count}"
    )

    print(
        f"Total user BGCs:         {len(status)}"
    )

    return status


# ============================================================
# Add singleton columns to distance table
# ============================================================

def add_singleton_columns(
    distance,
    status
):

    """
    Add singleton information for record_a and record_b.
    """

    status_lookup = status[
        [
            "id",
            "singleton",
            "status",
            "n_neighbors_at_threshold",
            "min_distance"
        ]
    ].copy()

    # --------------------------------------------------------
    # A
    # --------------------------------------------------------

    lookup_a = status_lookup.rename(
        columns={
            "id": "record_a_id",
            "singleton": "record_a_singleton",
            "status": "record_a_status",
            "n_neighbors_at_threshold":
                "record_a_n_neighbors_at_threshold",
            "min_distance":
                "record_a_min_distance"
        }
    )

    distance = distance.merge(
        lookup_a,
        on="record_a_id",
        how="left"
    )

    # --------------------------------------------------------
    # B
    # --------------------------------------------------------

    lookup_b = status_lookup.rename(
        columns={
            "id": "record_b_id",
            "singleton": "record_b_singleton",
            "status": "record_b_status",
            "n_neighbors_at_threshold":
                "record_b_n_neighbors_at_threshold",
            "min_distance":
                "record_b_min_distance"
        }
    )

    distance = distance.merge(
        lookup_b,
        on="record_b_id",
        how="left"
    )

    return distance


# ============================================================
# Remove MIBiG × MIBiG
# ============================================================

def remove_mibig_mibig(distance):

    print("\n==============================================")
    print("Removing MIBiG × MIBiG comparisons")
    print("==============================================")

    a_is_mibig = (
        distance[
            "record_a_id_path"
        ]
        .apply(is_mibig)
    )

    b_is_mibig = (
        distance[
            "record_b_id_path"
        ]
        .apply(is_mibig)
    )

    remove_mask = (
        a_is_mibig &
        b_is_mibig
    )

    removed = int(
        remove_mask.sum()
    )

    distance = distance[
        ~remove_mask
    ].copy()

    print(
        f"Removed MIBiG × MIBiG: {removed}"
    )

    print(
        f"Rows remaining:       {len(distance)}"
    )

    return distance


# ============================================================
# Keep only cand_cluster × cand_cluster
# ============================================================

def keep_cand_clusters(distance):

    print("\n==============================================")
    print("Keeping cand_cluster × cand_cluster comparisons")
    print("==============================================")

    a_is_candidate = (
        distance[
            "record_a_record_type"
        ]
        == "cand_cluster"
    )

    b_is_candidate = (
        distance[
            "record_b_record_type"
        ]
        == "cand_cluster"
    )

    keep_mask = (
        a_is_candidate &
        b_is_candidate
    )

    removed = int(
        (~keep_mask).sum()
    )

    distance = distance[
        keep_mask
    ].copy()

    print(
        f"Non-cand_cluster rows removed: "
        f"{removed}"
    )

    print(
        f"Rows remaining: "
        f"{len(distance)}"
    )

    return distance


# ============================================================
# Remove all MIBiG comparisons
# ============================================================

def remove_all_mibig(distance):

    """
    Keep only user-user comparisons.

    Any row where either side belongs to MIBiG is removed.
    """

    a_is_mibig = (
        distance[
            "record_a_id_path"
        ]
        .apply(is_mibig)
    )

    b_is_mibig = (
        distance[
            "record_b_id_path"
        ]
        .apply(is_mibig)
    )

    remove_mask = (
        a_is_mibig |
        b_is_mibig
    )

    return distance[
        ~remove_mask
    ].copy()


# ============================================================
# Add database labels
# ============================================================

def add_database(
    input_path,
    distance
):

    print("\n==============================================")
    print("Adding Database information")
    print("==============================================")

    distance = distance.copy()

    distance["record_a_database"] = (
        distance[
            "record_a_id_path"
        ]
        .apply(
            lambda x:
            get_database_from_path(
                x,
                input_path
            )
        )
    )

    distance["record_b_database"] = (
        distance[
            "record_b_id_path"
        ]
        .apply(
            lambda x:
            get_database_from_path(
                x,
                input_path
            )
        )
    )

    print(
        "Database names added."
    )

    return distance


# ============================================================
# Identify missing user-user pairs
# ============================================================

def identify_missing_user_pairs(
    distance,
    user_candidates
):

    """
    Compare all possible user-user BGC pairs with the pairs
    actually present in the BiG-SCAPE distance table.

    Missing pairs are reported but NOT assigned a distance.
    """

    print("\n==============================================")
    print("Checking user-user pair completeness")
    print("==============================================")

    user_ids = sorted(
        pd.to_numeric(
            user_candidates["id"],
            errors="coerce"
        )
        .dropna()
        .astype(int)
        .unique()
    )

    # --------------------------------------------------------
    # Observed pairs
    # --------------------------------------------------------

    observed_pairs = set()

    d = distance.copy()

    d["record_a_id"] = pd.to_numeric(
        d["record_a_id"],
        errors="coerce"
    )

    d["record_b_id"] = pd.to_numeric(
        d["record_b_id"],
        errors="coerce"
    )

    d = d.dropna(
        subset=[
            "record_a_id",
            "record_b_id"
        ]
    )

    for row in d.itertuples(index=False):

        a = int(row.record_a_id)
        b = int(row.record_b_id)

        if a == b:
            continue

        if (
            a in user_ids and
            b in user_ids
        ):
            observed_pairs.add(
                tuple(sorted((a, b)))
            )

    # --------------------------------------------------------
    # All possible pairs
    # --------------------------------------------------------

    possible_pairs = set()

    for i in range(len(user_ids)):

        for j in range(i + 1, len(user_ids)):

            possible_pairs.add(
                (
                    user_ids[i],
                    user_ids[j]
                )
            )

    missing_pairs = (
        possible_pairs -
        observed_pairs
    )

    print(
        f"Possible user-user pairs: "
        f"{len(possible_pairs)}"
    )

    print(
        f"Observed user-user pairs: "
        f"{len(observed_pairs)}"
    )

    print(
        f"Missing user-user pairs: "
        f"{len(missing_pairs)}"
    )

    # --------------------------------------------------------
    # Create output table
    # --------------------------------------------------------

    # Map BGC IDs to GBK IDs.  BiG-SCAPE intentionally does not
    # generate distances between records belonging to the same GBK.
    gbk_lookup = {}
    if "gbk_id" in user_candidates.columns:
        tmp = user_candidates[["id", "gbk_id"]].copy()
        tmp["id"] = pd.to_numeric(tmp["id"], errors="coerce")
        for row in tmp.dropna(subset=["id"]).itertuples(index=False):
            gbk_lookup[int(row.id)] = row.gbk_id

    rows = []

    for a, b in sorted(missing_pairs):

        gbk_a = gbk_lookup.get(a, np.nan)
        gbk_b = gbk_lookup.get(b, np.nan)

        rows.append(
            {
                "record_a_id": a,
                "record_b_id": b,
                "record_a_gbk_id": gbk_a,
                "record_b_gbk_id": gbk_b,
                "same_gbk_id": (
                    pd.notna(gbk_a) and
                    pd.notna(gbk_b) and
                    str(gbk_a) == str(gbk_b)
                ),
                "pair_status": "missing"
            }
        )

    return pd.DataFrame(rows)


# ============================================================
# Select one representative BGC per affected GBK
# ============================================================

def select_gbk_representatives_from_missing_pairs(
    missing_pairs,
    user_candidates
):
    """
    Select one BGC representative for each GBK implicated in
    same-GBK missing comparisons.

    BiG-SCAPE 2.0.3 intentionally skips comparisons between BGC
    records whose ``parent_gbk`` is the same. Therefore these
    missing pairs do not represent failed distance calculations.

    For every GBK that appears in a same-GBK missing pair, all user
    cand_cluster records from that GBK are considered and the record
    with the largest genomic span (End - Start) is retained. All
    other records from that GBK are excluded from matrix/statistical
    analyses.

    The function is deliberately applied only to GBKs implicated by
    the named same-GBK missing pairs. BGCs from unaffected GBKs are
    retained.

    BiG-SCAPE 2 stores BGC coordinates in ``bgc_record.nt_start``
    and ``bgc_record.nt_stop``; those are accepted directly. Common
    ``Start``/``End`` aliases are also supported for compatibility.
    """

    if missing_pairs is None or missing_pairs.empty:
        return set(pd.to_numeric(user_candidates["id"], errors="coerce")
                   .dropna().astype(int)), pd.DataFrame()

    candidates = user_candidates.copy()
    candidates["id"] = pd.to_numeric(candidates["id"], errors="coerce")

    # --------------------------------------------------------
    # Identify same-GBK missing pairs
    # --------------------------------------------------------
    if "same_gbk_id" in missing_pairs.columns:
        same = missing_pairs[
            missing_pairs["same_gbk_id"] == True
        ].copy()
    else:
        # Backward-compatible fallback.
        same = missing_pairs.copy()
        if {"record_a_gbk_id", "record_b_gbk_id"}.issubset(same.columns):
            same = same[
                same["record_a_gbk_id"].notna() &
                same["record_b_gbk_id"].notna() &
                (same["record_a_gbk_id"].astype(str) ==
                 same["record_b_gbk_id"].astype(str))
            ]
        else:
            same = same.iloc[0:0]

    if same.empty:
        return set(candidates["id"].dropna().astype(int)), pd.DataFrame()

    affected_gbks = set()
    for col in ["record_a_gbk_id", "record_b_gbk_id"]:
        if col in same.columns:
            affected_gbks.update(
                same[col].dropna().tolist()
            )

    if not affected_gbks or "gbk_id" not in candidates.columns:
        raise ValueError(
            "Same-GBK missing pairs were detected, but gbk_id information "
            "is unavailable in user_candidates."
        )

    # --------------------------------------------------------
    # Find coordinate columns
    # --------------------------------------------------------
    def find_column(frame, names):
        normalized = {str(c).strip().lower(): c for c in frame.columns}
        for name in names:
            if name.lower() in normalized:
                return normalized[name.lower()]
        return None

    start_col = find_column(
        candidates,
        [
            "Start", "start",
            "nt_start",
            "record_start",
            "bgc_start",
            "start_position"
        ]
    )
    end_col = find_column(
        candidates,
        [
            "End", "end",
            "nt_stop",
            "record_end",
            "bgc_end",
            "end_position"
        ]
    )

    if start_col is None or end_col is None:
        raise ValueError(
            "Cannot select the largest BGC representative because Start/End "
            "coordinates were not found in bgc_record/user_candidates. "
            f"Available columns: {list(candidates.columns)}"
        )

    candidates[start_col] = pd.to_numeric(
        candidates[start_col], errors="coerce"
    )
    candidates[end_col] = pd.to_numeric(
        candidates[end_col], errors="coerce"
    )

    candidates["_Krill_BGC_size"] = (
        candidates[end_col] - candidates[start_col]
    ).abs()

    # --------------------------------------------------------
    # Choose the largest BGC from each affected GBK
    # --------------------------------------------------------
    affected = candidates[
        candidates["gbk_id"].isin(affected_gbks)
    ].copy()

    representatives = []
    excluded = []

    for gbk_id, group in affected.groupby("gbk_id", dropna=False):
        group = group.sort_values(
            ["_Krill_BGC_size", "id"],
            ascending=[False, True],
            na_position="last"
        )

        if group.empty:
            continue

        if group["_Krill_BGC_size"].notna().sum() == 0:
            raise ValueError(
                f"Cannot select a largest BGC for GBK {gbk_id}: "
                "all affected cand_cluster records have missing nt_start/nt_stop."
            )

        rep = group.iloc[0]
        rep_id = int(rep["id"])
        representatives.append(rep_id)

        for _, row in group.iloc[1:].iterrows():
            excluded.append({
                "record_id": int(row["id"]),
                "gbk_id": gbk_id,
                "size": row["_Krill_BGC_size"],
                "representative_record_id": rep_id,
                "reason": "same_gbk_smaller_than_representative"
            })

        print(
            f"GBK {gbk_id}: representative BGC {rep_id} "
            f"(size={rep['_Krill_BGC_size']})"
        )

    # All unaffected GBKs retain all their BGCs. Affected GBKs retain
    # only the largest representative selected above.
    affected_gbk_ids = set(affected["gbk_id"].tolist())
    all_ids = set(candidates["id"].dropna().astype(int))

    unaffected_ids = set(
        candidates.loc[
            ~candidates["gbk_id"].isin(affected_gbk_ids),
            "id"
        ].dropna().astype(int)
    )

    representative_ids = unaffected_ids | set(representatives)

    excluded_df = pd.DataFrame(excluded)
    if not excluded_df.empty:
        excluded_df = excluded_df.sort_values(
            ["gbk_id", "record_id"]
        )

    print(
        f"\nAffected GBKs with same-GBK missing pairs: "
        f"{len(affected_gbk_ids)}"
    )
    print(
        f"Representative BGCs retained from affected GBKs: "
        f"{len(representatives)}"
    )
    print(
        f"Smaller same-GBK BGCs excluded: "
        f"{len(excluded_df)}"
    )
    print(
        f"Total matrix BGC representatives retained: "
        f"{len(representative_ids)}"
    )

    return representative_ids, excluded_df


# ============================================================
# Complete-case IDs
# ============================================================

def get_complete_case_ids(
    distance,
    user_candidates
):

    """
    Iteratively remove BGCs that have missing pairwise
    comparisons.

    The result is the largest set obtained by repeatedly
    removing BGCs involved in missing user-user pairs.

    No missing distances are imputed.
    """

    ids = set(
        pd.to_numeric(
            user_candidates["id"],
            errors="coerce"
        )
        .dropna()
        .astype(int)
    )

    d = distance.copy()

    d["record_a_id"] = pd.to_numeric(
        d["record_a_id"],
        errors="coerce"
    )

    d["record_b_id"] = pd.to_numeric(
        d["record_b_id"],
        errors="coerce"
    )

    d = d.dropna(
        subset=[
            "record_a_id",
            "record_b_id"
        ]
    )

    # --------------------------------------------------------
    # Observed unordered pairs
    # --------------------------------------------------------

    observed_pairs = set()

    for row in d.itertuples(index=False):

        a = int(row.record_a_id)
        b = int(row.record_b_id)

        if a == b:
            continue

        if a in ids and b in ids:

            observed_pairs.add(
                tuple(sorted((a, b)))
            )

    # --------------------------------------------------------
    # Iteratively remove IDs involved in missing pairs
    # --------------------------------------------------------

    changed = True

    while changed:

        changed = False

        missing_for_current = []

        current_ids = sorted(ids)

        for i in range(len(current_ids)):

            for j in range(i + 1, len(current_ids)):

                pair = (
                    current_ids[i],
                    current_ids[j]
                )

                if pair not in observed_pairs:

                    missing_for_current.append(
                        pair
                    )

        if missing_for_current:

            to_remove = set()

            for a, b in missing_for_current:

                to_remove.add(a)
                to_remove.add(b)

            new_ids = ids - to_remove

            if new_ids != ids:

                ids = new_ids
                changed = True

    return sorted(ids)


# ============================================================
# Subset complete-case distance table
# ============================================================

def subset_complete_case_distance(
    distance,
    user_candidates
):

    

    complete_ids = get_complete_case_ids(
        distance,
        user_candidates
    )

    complete_set = set(
        complete_ids
    )

    d = distance.copy()

    print("\n distance d table from subset_complete_case_distance before function: \n ", d )


    '''
    d = d[
        d["record_a_id"].isin(
            complete_set
        ) &
        d["record_b_id"].isin(
            complete_set
        )
    ].copy()

    print("\n distance d table from subset_complete_case_distance : \n ", d )
    '''

    return d, complete_ids



# ============================================================
# Build distance matrix
# ============================================================

def build_distance_matrix(
    distance,
    user_only=True,
    ids=None
):

    """
    Build a symmetric observed-distance matrix.

    Missing pairwise distances remain NaN. They are never imputed.

    Parameters
    ----------
    distance : pandas.DataFrame
        Long-format BiG-SCAPE distance table.
    user_only : bool
        If True, remove all MIBiG comparisons before building the matrix.
    ids : iterable, optional
        Explicit BGC IDs to retain in the matrix. Supplying this is important
        for pairwise-available analyses because a BGC must remain in the
        matrix even when one or more of its pairwise distances are missing.
    """

    d = distance.copy()

    if user_only:
        d = remove_all_mibig(d)

    d["record_a_id"] = pd.to_numeric(
        d["record_a_id"],
        errors="coerce"
    )

    d["record_b_id"] = pd.to_numeric(
        d["record_b_id"],
        errors="coerce"
    )

    d["distance"] = pd.to_numeric(
        d["distance"],
        errors="coerce"
    )

    d = d.dropna(
        subset=[
            "record_a_id",
            "record_b_id",
            "distance"
        ]
    )

    d = d[d["record_a_id"] != d["record_b_id"]]

    if ids is None:
        ids = sorted(
            set(d["record_a_id"].astype(int))
            |
            set(d["record_b_id"].astype(int))
        )
    else:
        ids = sorted(
            set(
                pd.to_numeric(
                    pd.Series(list(ids)),
                    errors="coerce"
                ).dropna().astype(int)
            )
        )

    matrix = pd.DataFrame(
        np.nan,
        index=ids,
        columns=ids,
        dtype=float
    )

    # Use a writable copy. Newer pandas/NumPy combinations can expose
    # DataFrame.values as read-only.
    matrix_array = matrix.to_numpy(copy=True)
    np.fill_diagonal(matrix_array, 0.0)
    matrix = pd.DataFrame(
        matrix_array,
        index=ids,
        columns=ids
    )

    for row in d.itertuples(index=False):

        a = int(row.record_a_id)
        b = int(row.record_b_id)
        value = float(row.distance)

        if a not in matrix.index or b not in matrix.index:
            continue

        # Keep the minimum if duplicated.
        if pd.isna(matrix.loc[a, b]):
            matrix.loc[a, b] = value
            matrix.loc[b, a] = value
        else:
            matrix.loc[a, b] = min(
                matrix.loc[a, b],
                value
            )
            matrix.loc[b, a] = min(
                matrix.loc[b, a],
                value
            )

    return matrix


# ============================================================
# PCoA
# ============================================================

def pcoa(distance_matrix):

    """
    Principal Coordinates Analysis using Gower's centered matrix.
    """

    D = np.asarray(
        distance_matrix,
        dtype=float
    )

    n = D.shape[0]

    if n < 2:
        raise ValueError(
            "At least two BGCs are required for PCoA."
        )

    if np.isnan(D).any():
        raise ValueError(
            "Distance matrix contains missing values."
        )

    # Square distances
    D2 = D ** 2

    # Centering matrix
    J = (
        np.eye(n)
        -
        np.ones((n, n)) / n
    )

    B = (
        -0.5
        *
        J
        @ D2
        @ J
    )

    eigenvalues, eigenvectors = np.linalg.eigh(B)

    # Sort descending
    order = np.argsort(
        eigenvalues
    )[::-1]

    eigenvalues = eigenvalues[
        order
    ]

    eigenvectors = eigenvectors[
        :,
        order
    ]

    # Keep positive axes
    positive = (
        eigenvalues > 1e-10
    )

    eigenvalues = eigenvalues[
        positive
    ]

    eigenvectors = eigenvectors[
        :,
        positive
    ]

    if len(eigenvalues) == 0:

        raise ValueError(
            "No positive eigenvalues were found."
        )

    coordinates = (
        eigenvectors
        *
        np.sqrt(eigenvalues)
    )

    total = np.sum(
        eigenvalues
    )

    explained = (
        eigenvalues /
        total
    )

    return (
        coordinates,
        eigenvalues,
        explained
    )


# ============================================================
# PERMANOVA
# ============================================================

def permanova(
    distance_matrix,
    groups,
    permutations=999
):

    """
    One-factor PERMANOVA.

    Uses pseudo-F based on sums of squared distances.

    Permutations shuffle group labels.
    """

    D = np.asarray(
        distance_matrix,
        dtype=float
    )

    groups = np.asarray(
        groups
    )

    n = len(groups)

    if D.shape != (n, n):

        raise ValueError(
            "Distance matrix and groups have incompatible dimensions."
        )

    unique_groups = np.unique(
        groups
    )

    if len(unique_groups) < 2:

        raise ValueError(
            "PERMANOVA requires at least two groups."
        )

    # --------------------------------------------------------
    # Total sum of squares
    # --------------------------------------------------------

    total_ss = (
        np.sum(D ** 2)
        /
        n
    )

    # --------------------------------------------------------
    # Within-group sum of squares
    # --------------------------------------------------------

    within_ss = 0.0

    for group in unique_groups:

        idx = np.where(
            groups == group
        )[0]

        if len(idx) == 0:
            continue

        Dg = D[
            np.ix_(idx, idx)
        ]

        within_ss += (
            np.sum(Dg ** 2)
            /
            len(idx)
        )

    between_ss = (
        total_ss -
        within_ss
    )

    df_between = (
        len(unique_groups) -
        1
    )

    df_within = (
        n -
        len(unique_groups)
    )

    if df_within <= 0:

        raise ValueError(
            "Not enough observations within groups for PERMANOVA."
        )

    ms_between = (
        between_ss /
        df_between
    )

    ms_within = (
        within_ss /
        df_within
    )

    pseudo_F = (
        ms_between /
        ms_within
        if ms_within > 0
        else np.inf
    )

    R2 = (
        between_ss /
        total_ss
        if total_ss > 0
        else np.nan
    )

    # --------------------------------------------------------
    # Permutation test
    # --------------------------------------------------------

    rng = np.random.default_rng(
        12345
    )

    permuted_F = []

    for _ in range(
        permutations
    ):

        shuffled = rng.permutation(
            groups
        )

        perm_within = 0.0

        for group in unique_groups:

            idx = np.where(
                shuffled == group
            )[0]

            if len(idx) == 0:
                continue

            Dg = D[
                np.ix_(idx, idx)
            ]

            perm_within += (
                np.sum(Dg ** 2)
                /
                len(idx)
            )

        perm_between = (
            total_ss -
            perm_within
        )

        perm_ms_between = (
            perm_between /
            df_between
        )

        perm_ms_within = (
            perm_within /
            df_within
        )

        if perm_ms_within > 0:

            permuted_F.append(
                perm_ms_between /
                perm_ms_within
            )

    permuted_F = np.asarray(
        permuted_F
    )

    if len(permuted_F) > 0:

        p_value = (
            np.sum(
                permuted_F >= pseudo_F
            )
            + 1
        ) / (
            len(permuted_F)
            + 1
        )

    else:

        p_value = np.nan

    return {
        "pseudo_F": pseudo_F,
        "R2": R2,
        "p_value": p_value,
        "permutations": permutations
    }


# ============================================================
# PERMDISP
# ============================================================

def permdisp(
    distance_matrix,
    groups,
    permutations=999
):

    """
    PERMDISP-like test based on distances to group centroids.

    This is a companion to PERMANOVA because PERMANOVA can be
    influenced by differences in within-group dispersion.
    """

    D = np.asarray(
        distance_matrix,
        dtype=float
    )

    groups = np.asarray(
        groups
    )

    unique_groups = np.unique(
        groups
    )

    if len(unique_groups) < 2:

        raise ValueError(
            "PERMDISP requires at least two groups."
        )

    # --------------------------------------------------------
    # Convert distance matrix to Euclidean coordinates
    # --------------------------------------------------------

    coords, _, _ = pcoa(
        D
    )

    # --------------------------------------------------------
    # Distance to centroid
    # --------------------------------------------------------

    distances_to_centroid = []

    labels = []

    for group in unique_groups:

        idx = np.where(
            groups == group
        )[0]

        if len(idx) < 2:
            continue

        centroid = np.mean(
            coords[idx],
            axis=0
        )

        for i in idx:

            dist = np.linalg.norm(
                coords[i] -
                centroid
            )

            distances_to_centroid.append(
                dist
            )

            labels.append(
                group
            )

    distances_to_centroid = np.asarray(
        distances_to_centroid
    )

    labels = np.asarray(
        labels
    )

    # --------------------------------------------------------
    # Observed one-way ANOVA on centroid distances
    # --------------------------------------------------------

    group_values = [
        distances_to_centroid[
            labels == group
        ]
        for group in unique_groups
        if np.sum(
            labels == group
        ) > 0
    ]

    grand_mean = np.mean(
        distances_to_centroid
    )

    ss_between = 0.0
    ss_within = 0.0

    for values in group_values:

        group_mean = np.mean(
            values
        )

        ss_between += (
            len(values)
            *
            (group_mean - grand_mean) ** 2
        )

        ss_within += np.sum(
            (
                values -
                group_mean
            ) ** 2
        )

    df_between = (
        len(group_values) -
        1
    )

    df_within = (
        len(distances_to_centroid) -
        len(group_values)
    )

    if df_within <= 0:

        raise ValueError(
            "Not enough observations for PERMDISP."
        )

    ms_between = (
        ss_between /
        df_between
    )

    ms_within = (
        ss_within /
        df_within
    )

    F = (
        ms_between /
        ms_within
        if ms_within > 0
        else np.inf
    )

    # --------------------------------------------------------
    # Permutation
    # --------------------------------------------------------

    rng = np.random.default_rng(
        12345
    )

    permuted_F = []

    for _ in range(
        permutations
    ):

        shuffled = rng.permutation(
            labels
        )

        perm_group_values = [
            distances_to_centroid[
                shuffled == group
            ]
            for group in unique_groups
            if np.sum(
                shuffled == group
            ) > 0
        ]

        if len(perm_group_values) < 2:
            continue

        grand_mean_perm = np.mean(
            distances_to_centroid
        )

        ssb = 0.0
        ssw = 0.0

        for values in perm_group_values:

            group_mean = np.mean(
                values
            )

            ssb += (
                len(values)
                *
                (
                    group_mean -
                    grand_mean_perm
                ) ** 2
            )

            ssw += np.sum(
                (
                    values -
                    group_mean
                ) ** 2
            )

        dfb = (
            len(perm_group_values) -
            1
        )

        dfw = (
            len(distances_to_centroid) -
            len(perm_group_values)
        )

        if dfw > 0 and ssw > 0:

            permuted_F.append(
                (
                    ssb / dfb
                )
                /
                (
                    ssw / dfw
                )
            )

    if len(permuted_F) > 0:

        p_value = (
            np.sum(
                np.asarray(permuted_F) >= F
            )
            + 1
        ) / (
            len(permuted_F)
            + 1
        )

    else:

        p_value = np.nan

    return {
        "F": F,
        "p_value": p_value,
        "permutations": permutations
    }


# ============================================================
# Pairwise PERMANOVA
# ============================================================

def pairwise_permanova(
    distance_matrix,
    metadata,
    permutations=999
):

    """
    Pairwise PERMANOVA between databases.

    P-values are corrected using Benjamini-Hochberg FDR.
    """

    results = []

    groups = sorted(
        metadata.unique()
    )

    for i in range(
        len(groups)
    ):

        for j in range(
            i + 1,
            len(groups)
        ):

            group_a = groups[i]
            group_b = groups[j]

            selected = metadata.isin(
                [
                    group_a,
                    group_b
                ]
            )

            ids = metadata.index[
                selected
            ]

            if len(ids) < 2:
                continue

            submatrix = distance_matrix.loc[
                ids,
                ids
            ]

            labels = metadata.loc[
                ids
            ].values

            if len(
                np.unique(labels)
            ) < 2:

                continue

            try:

                result = permanova(
                    submatrix.values,
                    labels,
                    permutations=permutations
                )

                results.append(
                    {
                        "group_a": group_a,
                        "group_b": group_b,
                        "pseudo_F":
                            result["pseudo_F"],
                        "R2":
                            result["R2"],
                        "p_value":
                            result["p_value"]
                    }
                )

            except ValueError:

                continue

    result_df = pd.DataFrame(
        results
    )

    if result_df.empty:

        return result_df

    # --------------------------------------------------------
    # Benjamini-Hochberg correction
    # --------------------------------------------------------

    pvals = result_df[
        "p_value"
    ].values

    order = np.argsort(
        pvals
    )

    adjusted = np.empty_like(
        pvals,
        dtype=float
    )

    m = len(
        pvals
    )

    running_min = 1.0

    for rank_reverse, idx in enumerate(
        order[::-1],
        start=1
    ):

        rank = m - rank_reverse + 1

        value = (
            pvals[idx]
            *
            m
            /
            rank
        )

        running_min = min(
            running_min,
            value
        )

        adjusted[idx] = min(
            running_min,
            1.0
        )

    result_df[
        "p_adjust_BH"
    ] = adjusted

    return result_df


# ============================================================
# Plot PCoA / MDS
# ============================================================

def plot_ordination(
    matrix,
    metadata,
    output_prefix,
    method="PCoA"
):

    if len(matrix) < 2:
        return

    if method == "PCoA":

        coordinates, eigenvalues, explained = pcoa(
            matrix.values
        )

        if coordinates.shape[1] < 2:
            raise ValueError(
                "PCoA produced fewer than two axes."
            )

        x = coordinates[:, 0]
        y = coordinates[:, 1]

        xlabel = (
            f"PCoA1 "
            f"({explained[0] * 100:.1f}%)"
        )

        ylabel = (
            f"PCoA2 "
            f"({explained[1] * 100:.1f}%)"
            if len(explained) > 1
            else "PCoA2"
        )

    elif method == "MDS":

        model = MDS(
            n_components=2,
            dissimilarity="precomputed",
            random_state=12345
        )

        coordinates = model.fit_transform(
            matrix.values
        )

        x = coordinates[:, 0]
        y = coordinates[:, 1]

        xlabel = "MDS1"
        ylabel = "MDS2"

    else:

        raise ValueError(
            f"Unknown ordination method: {method}"
        )

    plot_df = pd.DataFrame(
        {
            "x": x,
            "y": y,
            "database": metadata.values
        },
        index=metadata.index
    )

    plt.figure(
        figsize=(9, 7)
    )

    sns.scatterplot(
        data=plot_df,
        x="x",
        y="y",
        hue="database",
        s=80
    )

    plt.xlabel(
        xlabel
    )

    plt.ylabel(
        ylabel
    )

    plt.title(
        f"{method} of user BGCs"
    )

    plt.tight_layout()

    output_file = (
        output_prefix
        +
        f"_{method}.png"
    )

    plt.savefig(
        output_file,
        dpi=300
    )

    plt.close()

    print(
        f"Saved {method} plot: "
        f"{output_file}"
    )


# ============================================================
# Complete statistical analysis
# ============================================================

def run_complete_statistics(
    distance,
    output_prefix,
    label="analysis",
    permutations=999
):

    print("\n==============================================")
    print(
        f"Running statistical analysis: {label}"
    )
    print("==============================================")

    print("\n arquivo distance utilizado: \n ", distance)

    matrix = build_distance_matrix(
        distance,
        user_only=True
    )

    # Remove incomplete rows/columns just in case
    complete_mask = ~matrix.isna().any(
        axis=1
    )

    matrix = matrix.loc[
        complete_mask,
        complete_mask
    ]

    if len(matrix) < 2:

        raise ValueError(
            "Not enough BGCs for statistical analysis."
        )

    metadata_a = distance[
        ["record_a_id", "record_a_database"]
    ].rename(
        columns={
            "record_a_id": "record_id",
            "record_a_database": "database"
        }
    )

    metadata_b = distance[
        ["record_b_id", "record_b_database"]
    ].rename(
        columns={
            "record_b_id": "record_id",
            "record_b_database": "database"
        }
    )

    metadata = pd.concat(
        [metadata_a, metadata_b],
        ignore_index=True
    ).drop_duplicates(
        subset=["record_id"]
    ).set_index("record_id")["database"]


    '''
    metadata = pd.Series(
        distance[
            [
                "record_a_id",
                "record_a_database"
            ]
        ]
        .drop_duplicates()
        .set_index(
            "record_a_id"
        )[
            "record_a_database"
        ],
        name="database"
    )'''

    #print("\n statistics metadata: \n  ", metadata)

    metadata = metadata.reindex(
        matrix.index
    )

    valid = ~metadata.isna()

    matrix = matrix.loc[
        valid,
        valid
    ]

    metadata = metadata.loc[
        valid
    ]

    if metadata.nunique() < 2:

        raise ValueError(
            "At least two databases/groups are required."
        )

    #print("\n statistics metadata second: \n  ", metadata)

    # --------------------------------------------------------
    # PCoA
    # --------------------------------------------------------

    plot_ordination(
        matrix,
        metadata,
        output_prefix + "_" + label,
        method="PCoA"
    )

    # --------------------------------------------------------
    # MDS
    # --------------------------------------------------------

    plot_ordination(
        matrix,
        metadata,
        output_prefix + "_" + label,
        method="MDS"
    )

    # --------------------------------------------------------
    # PERMANOVA
    # --------------------------------------------------------

    perm = permanova(
        matrix.values,
        metadata.values,
        permutations=permutations
    )

    # --------------------------------------------------------
    # PERMDISP
    # --------------------------------------------------------

    disp = permdisp(
        matrix.values,
        metadata.values,
        permutations=permutations
    )

    # --------------------------------------------------------
    # Pairwise PERMANOVA
    # --------------------------------------------------------

    pairwise = pairwise_permanova(
        matrix,
        metadata,
        permutations=permutations
    )

    pairwise_file = (
        output_prefix
        +
        f"_{label}_pairwise_PERMANOVA.tsv"
    )

    pairwise.to_csv(
        pairwise_file,
        sep="\t",
        index=False
    )

    # --------------------------------------------------------
    # Summary
    # --------------------------------------------------------

    stats_summary = pd.DataFrame(
        [
            {
                "analysis": label,
                "n_BGCs": len(matrix),
                "n_databases":
                    metadata.nunique(),
                "PERMANOVA_p":
                    perm["p_value"],
                "PERMANOVA_R2":
                    perm["R2"],
                "PERMANOVA_pseudo_F":
                    perm["pseudo_F"],
                "PERMDISP_p":
                    disp["p_value"],
                "PERMDISP_F":
                    disp["F"]
            }
        ]
    )

    stats_file = (
        output_prefix
        +
        f"_{label}_statistics.tsv"
    )

    stats_summary.to_csv(
        stats_file,
        sep="\t",
        index=False
    )

    print("\nStatistical results")
    print("----------------------------------------------")

    print(
        f"PERMANOVA pseudo-F: "
        f"{perm['pseudo_F']:.4f}"
    )

    print(
        f"PERMANOVA R2:       "
        f"{perm['R2']:.4f}"
    )

    print(
        f"PERMANOVA p:        "
        f"{perm['p_value']:.5f}"
    )

    print(
        f"PERMDISP F:         "
        f"{disp['F']:.4f}"
    )

    print(
        f"PERMDISP p:         "
        f"{disp['p_value']:.5f}"
    )

    print(
        f"Statistics file:    "
        f"{stats_file}"
    )

    print(
        f"Pairwise file:      "
        f"{pairwise_file}"
    )

    return (
        matrix,
        metadata,
        perm,
        disp,
        pairwise
    )


# ============================================================
# Legacy plotting wrapper
# ============================================================

def plot_bgc_distribution(
    distance_noMIBIG,
    output_path=None,
    distance_threshold=0.95
):

    """
    Backward-compatible plotting function.
    """

    if output_path is None:
        output_path = "bgc_distribution"

    prefix = output_path

    matrix = build_distance_matrix(
        distance_noMIBIG,
        user_only=True
    )

    metadata = pd.Series(
        distance_noMIBIG[
            [
                "record_a_id",
                "record_a_database"
            ]
        ]
        .drop_duplicates()
        .set_index(
            "record_a_id"
        )[
            "record_a_database"
        ]
    )

    metadata = metadata.reindex(
        matrix.index
    )

    plot_ordination(
        matrix,
        metadata,
        prefix,
        method="PCoA"
    )

    plot_ordination(
        matrix,
        metadata,
        prefix,
        method="MDS"
    )


# ============================================================
# Save outputs
# ============================================================

def save_outputs(
    distance,
    bgc_record,
    status,
    output_prefix,
    distance_filtered,
    observed_user=None
):

    distance_file = (
        output_prefix +
        "_distance_annotated.tsv"
    )

    bgc_record_file = (
        output_prefix +
        "_bgc_record.tsv"
    )

    status_file = (
        output_prefix +
        "_singleton_classification.tsv"
    )

    distance_filtered_file = (
        output_prefix +
        "_distance_cand_cluster.tsv"
    )

    distance.to_csv(
        distance_file,
        sep="\t",
        index=False
    )

    bgc_record.to_csv(
        bgc_record_file,
        sep="\t",
        index=False
    )

    status.to_csv(
        status_file,
        sep="\t",
        index=False
    )

    distance_filtered.to_csv(
        distance_filtered_file,
        sep="\t",
        index=False
    )

    if observed_user is not None:

        observed_file = (
            output_prefix +
            "_distance_observed_user_only.tsv"
        )

        observed_user.to_csv(
            observed_file,
            sep="\t",
            index=False
        )

    print("\n==============================================")
    print("Core output files")
    print("==============================================")

    print(
        distance_file
    )

    print(
        bgc_record_file
    )

    print(
        status_file
    )

    print(
        distance_filtered_file
    )

    if observed_user is not None:

        print(
            observed_file
        )



# ============================================================
# BGC CATEGORY-SPECIFIC DATASET COMPARISONS
# ============================================================

# These are the seven broad BGC classes requested for the
# category-specific analyses.
#
# IMPORTANT:
# BiG-SCAPE can store combined categories such as:
#   NRPS.PKS
#   NRPS.PKS.RiPP
#   PKS.other.terpene
#
# A BGC with a combined category is included in EVERY broad
# category represented in that string. For example:
#
#   NRPS.PKS.RiPP
#
# is included in:
#   NRPS
#   PKS
#   RiPP
#
# The category matrices therefore intentionally overlap.
CATEGORY_NAMES = [
    "RiPP",
    "terpene",
    "PKS",
    "NRPS",
    "other",
    "NRPS.PKS",
    "saccharide"
]


def category_contains(category_value, requested_category):
    """
    Return True when requested_category is represented in a
    BiG-SCAPE category string.

    Examples
    --------
    NRPS.PKS.RiPP contains:
        NRPS
        PKS
        RiPP

    It does NOT treat "NRPS.PKS" as a substring category only;
    the category string is parsed on '.' and every token is
    considered a category.

    The requested category "NRPS.PKS" is therefore handled as
    a special combined class below.
    """

    if pd.isna(category_value):
        return False

    value = str(category_value).strip()

    if value == "":
        return False

    tokens = set(
        token.strip()
        for token in value.split(".")
        if token.strip()
    )

    # "NRPS.PKS" is explicitly requested as its own class.
    # A BGC belongs to this class when BOTH NRPS and PKS are
    # represented in its BiG-SCAPE category.
    if requested_category == "NRPS.PKS":
        return (
            "NRPS" in tokens and
            "PKS" in tokens
        )

    return requested_category in tokens


def get_category_bgc_ids(
    distance,
    user_candidates,
    category
):
    """
    Identify user BGC IDs belonging to one broad BGC category.

    Category membership is determined from BOTH:
        record_a_category
        record_b_category

    This is important because a BGC can occur on either side of
    a distance comparison.

    The user_candidates table is also used as a fallback so that
    a user BGC that has no observed distance row is not silently
    assigned to the wrong category. No distance is created for
    such a BGC.
    """

    category_ids = set()

    # --------------------------------------------------------
    # A-side category membership
    # --------------------------------------------------------

    if (
        "record_a_category" in distance.columns and
        "record_a_id" in distance.columns
    ):

        a_mask = distance[
            "record_a_category"
        ].apply(
            lambda x:
            category_contains(
                x,
                category
            )
        )

        a_ids = pd.to_numeric(
            distance.loc[
                a_mask,
                "record_a_id"
            ],
            errors="coerce"
        ).dropna().astype(int)

        category_ids.update(
            a_ids.tolist()
        )

    # --------------------------------------------------------
    # B-side category membership
    # --------------------------------------------------------

    if (
        "record_b_category" in distance.columns and
        "record_b_id" in distance.columns
    ):

        b_mask = distance[
            "record_b_category"
        ].apply(
            lambda x:
            category_contains(
                x,
                category
            )
        )

        b_ids = pd.to_numeric(
            distance.loc[
                b_mask,
                "record_b_id"
            ],
            errors="coerce"
        ).dropna().astype(int)

        category_ids.update(
            b_ids.tolist()
        )

    # --------------------------------------------------------
    # Restrict to user BGCs.
    #
    # This prevents MIBiG IDs from entering category matrices.
    # --------------------------------------------------------

    user_ids = set(
        pd.to_numeric(
            user_candidates["id"],
            errors="coerce"
        ).dropna().astype(int)
    )

    category_ids &= user_ids

    # --------------------------------------------------------
    # Fallback from user_candidates.category.
    #
    # This is useful for a BGC with no observed distance row.
    # It does not create any distance value.
    # --------------------------------------------------------

    if (
        "category" in user_candidates.columns
    ):

        candidate_mask = user_candidates[
            "category"
        ].apply(
            lambda x:
            category_contains(
                x,
                category
            )
        )

        fallback_ids = pd.to_numeric(
            user_candidates.loc[
                candidate_mask,
                "id"
            ],
            errors="coerce"
        ).dropna().astype(int)

        category_ids.update(
            fallback_ids.tolist()
        )

    return sorted(category_ids)


def build_category_distance_table(
    distance,
    category_ids
):
    """
    Keep only observed user-user distance comparisons in which
    BOTH BGCs belong to the requested category.

    No missing pair is created and no distance is imputed.

    If a BGC belongs to multiple categories, its observed
    comparisons are consequently present in multiple category
    tables/matrices.
    """

    if len(category_ids) == 0:
        return distance.iloc[0:0].copy()

    category_set = set(
        int(x)
        for x in category_ids
    )

    d = distance.copy()

    d["record_a_id"] = pd.to_numeric(
        d["record_a_id"],
        errors="coerce"
    )

    d["record_b_id"] = pd.to_numeric(
        d["record_b_id"],
        errors="coerce"
    )

    d["distance"] = pd.to_numeric(
        d["distance"],
        errors="coerce"
    )

    d = d.dropna(
        subset=[
            "record_a_id",
            "record_b_id",
            "distance"
        ]
    ).copy()

    d["record_a_id"] = d[
        "record_a_id"
    ].astype(int)

    d["record_b_id"] = d[
        "record_b_id"
    ].astype(int)

    # Only comparisons between two BGCs that belong to this
    # category are retained.
    d = d[
        d["record_a_id"].isin(category_set) &
        d["record_b_id"].isin(category_set)
    ].copy()

    return d


def build_category_matrix(
    distance,
    category_ids
):
    """
    Build a symmetric matrix for one BGC category.

    All category BGC IDs are included, including IDs that have
    missing comparisons. Missing pairwise distances remain NaN.

    Diagonal values are 0.0 because self-distance is defined as
    zero for the matrix representation; this does NOT create a
    missing pair or modify BiG-SCAPE's original distance table.
    """

    ids = sorted(
        int(x)
        for x in category_ids
    )

    matrix = pd.DataFrame(
        np.nan,
        index=ids,
        columns=ids,
        dtype=float
    )

    if len(ids) == 0:
        return matrix

    matrix_array = matrix.to_numpy(copy=True)
    np.fill_diagonal(
        matrix_array,
        0.0
    )
    matrix = pd.DataFrame(
        matrix_array,
        index=ids,
        columns=ids
    )

    d = build_category_distance_table(
        distance,
        ids
    )

    for row in d.itertuples(index=False):

        a = int(row.record_a_id)
        b = int(row.record_b_id)

        if a == b:
            continue

        value = float(row.distance)

        # If duplicate rows exist, retain the smallest observed
        # BiG-SCAPE distance, exactly as build_distance_matrix().
        if pd.isna(matrix.loc[a, b]):

            matrix.loc[a, b] = value
            matrix.loc[b, a] = value

        else:

            matrix.loc[a, b] = min(
                matrix.loc[a, b],
                value
            )

            matrix.loc[b, a] = min(
                matrix.loc[b, a],
                value
            )

    return matrix


def build_category_metadata(
    distance,
    category_ids
):
    """
    Build BGC -> database metadata using BOTH record_a and
    record_b information.

    This avoids the A-side-only metadata issue that can occur
    when the first occurrence of a BGC happens to be on side B.
    """

    category_set = set(
        int(x)
        for x in category_ids
    )

    metadata_a = distance[
        [
            "record_a_id",
            "record_a_database"
        ]
    ].rename(
        columns={
            "record_a_id": "record_id",
            "record_a_database": "database"
        }
    )

    metadata_b = distance[
        [
            "record_b_id",
            "record_b_database"
        ]
    ].rename(
        columns={
            "record_b_id": "record_id",
            "record_b_database": "database"
        }
    )

    metadata = pd.concat(
        [
            metadata_a,
            metadata_b
        ],
        ignore_index=True
    )

    metadata["record_id"] = pd.to_numeric(
        metadata["record_id"],
        errors="coerce"
    )

    metadata = metadata.dropna(
        subset=["record_id"]
    )

    metadata["record_id"] = metadata[
        "record_id"
    ].astype(int)

    metadata = metadata[
        metadata["record_id"].isin(
            category_set
        )
    ]

    # Prefer a non-missing database label when duplicate
    # records occur.
    metadata = (
        metadata
        .dropna(subset=["database"])
        .drop_duplicates(
            subset=["record_id"],
            keep="first"
        )
        .set_index("record_id")["database"]
    )

    return metadata


def save_category_matrix(
    matrix,
    category,
    output_prefix
):
    """
    Save one category distance matrix.
    """

    safe_category = (
        category
        .replace(".", "_")
        .replace("/", "_")
        .replace(" ", "_")
    )

    matrix_file = (
        output_prefix +
        f"_category_{safe_category}_distance_matrix.tsv"
    )

    matrix.to_csv(
        matrix_file,
        sep="\t",
        index=True
    )

    return matrix_file


def save_category_distance_table(
    distance,
    category,
    output_prefix
):
    """
    Save the observed long-format distance table for one
    category.
    """

    safe_category = (
        category
        .replace(".", "_")
        .replace("/", "_")
        .replace(" ", "_")
    )

    distance_file = (
        output_prefix +
        f"_category_{safe_category}_distance.tsv"
    )

    distance.to_csv(
        distance_file,
        sep="\t",
        index=False
    )

    return distance_file


def save_category_membership(
    user_candidates,
    category_ids,
    category,
    output_prefix
):
    """
    Save the BGC membership list for one category.

    This makes the intentional overlap between categories
    explicit.
    """

    category_set = set(
        int(x)
        for x in category_ids
    )

    membership = user_candidates.copy()

    membership["record_id"] = pd.to_numeric(
        membership["id"],
        errors="coerce"
    )

    membership = membership[
        membership["record_id"].isin(
            category_set
        )
    ].copy()

    membership["analysis_category"] = category

    safe_category = (
        category
        .replace(".", "_")
        .replace("/", "_")
        .replace(" ", "_")
    )

    membership_file = (
        output_prefix +
        f"_category_{safe_category}_BGC_membership.tsv"
    )

    membership.to_csv(
        membership_file,
        sep="\t",
        index=False
    )

    return membership_file


def run_category_comparisons(
    observed_user,
    user_candidates,
    output_prefix,
    permutations=999,
    representative_ids=None
):
    """
    Run dataset-comparison analyses separately for each broad BGC category.

    Category BGCs are first restricted to the representative set selected
    from same-GBK missing pairs. This guarantees that intentionally undefined
    same-GBK BiG-SCAPE comparisons do not create missing values in the
    statistical matrix.

    If unexpected missing distances remain between the retained representatives,
    a complete-case subset is used for the classical statistical analysis and
    the missing pairs are reported rather than imputed.
    """

    print("\n\n")
    print("============================================================")
    print("BGC CATEGORY-SPECIFIC DATASET COMPARISONS")
    print("============================================================")

    category_summary_rows = []

    for category in CATEGORY_NAMES:

        print("\n------------------------------------------------------------")
        print(f"Category analysis: {category}")
        print("------------------------------------------------------------")

        category_ids = get_category_bgc_ids(
            observed_user,
            user_candidates,
            category
        )

        # Apply the same one-representative-per-affected-GBK rule
        # used by the global matrix.
        if representative_ids is not None:
            category_ids = sorted(
                set(category_ids) & set(representative_ids)
            )

        print(
            f"User BGC representatives in {category}: "
            f"{len(category_ids)}"
        )

        membership_file = save_category_membership(
            user_candidates,
            category_ids,
            category,
            output_prefix
        )

        category_distance = build_category_distance_table(
            observed_user,
            category_ids
        )

        print(
            f"Observed {category} distance rows: "
            f"{len(category_distance)}"
        )

        distance_file = save_category_distance_table(
            category_distance,
            category,
            output_prefix
        )

        category_matrix = build_category_matrix(
            observed_user,
            category_ids
        )

        matrix_file = save_category_matrix(
            category_matrix,
            category,
            output_prefix
        )

        possible_pairs = (
            len(category_ids)
            * (len(category_ids) - 1)
            // 2
        )

        if len(category_matrix) > 1:
            n_observed_pairs = int(
                np.triu(
                    np.isfinite(category_matrix.values),
                    k=1
                ).sum()
            )
        else:
            n_observed_pairs = 0

        n_missing_pairs = possible_pairs - n_observed_pairs

        print(
            f"Possible {category} pairs: {possible_pairs}"
        )
        print(
            f"Observed {category} pairs: {n_observed_pairs}"
        )
        print(
            f"Missing {category} pairs after representative filtering: "
            f"{n_missing_pairs}"
        )

        safe_category = (
            category
            .replace(".", "_")
            .replace("/", "_")
            .replace(" ", "_")
        )

        # ----------------------------------------------------
        # Complete matrix / classical analysis
        # ----------------------------------------------------

        if len(category_matrix) >= 2:
            complete_mask = ~category_matrix.isna().any(axis=1)
            complete_matrix = category_matrix.loc[
                complete_mask,
                complete_mask
            ]
        else:
            complete_matrix = category_matrix.copy()

        complete_matrix_file = None

        if len(complete_matrix) > 0:
            complete_matrix_file = (
                output_prefix
                + f"_category_{safe_category}_complete_case_matrix.tsv"
            )

            complete_matrix.to_csv(
                complete_matrix_file,
                sep="\t",
                index=True
            )

        print(
            f"Complete-case {category} BGCs: "
            f"{len(complete_matrix)}"
        )

        category_results = None

        if len(complete_matrix) >= 2:

            complete_ids = set(
                complete_matrix.index.astype(int)
            )

            category_complete_distance = category_distance[
                category_distance["record_a_id"].isin(complete_ids)
                & category_distance["record_b_id"].isin(complete_ids)
            ].copy()

            try:
                category_results = run_complete_statistics(
                    category_complete_distance,
                    output_prefix,
                    label=(
                        "category_"
                        + safe_category
                        + "_complete_case"
                    ),
                    permutations=permutations
                )
            except (
                ValueError,
                np.linalg.LinAlgError
            ) as e:
                print(
                    f"{category} complete-case statistical "
                    f"analysis could not be completed:"
                )
                print(f"  {e}")

        # ----------------------------------------------------
        # Category summary
        # ----------------------------------------------------

        category_row = {
            "category": category,
            "n_category_BGCs": len(category_ids),
            "n_possible_pairs": possible_pairs,
            "n_observed_pairs": n_observed_pairs,
            "n_missing_pairs": n_missing_pairs,
            "n_complete_case_BGCs": len(complete_matrix),
            "n_complete_case_databases": np.nan,
            "complete_case_PERMANOVA_p": np.nan,
            "complete_case_PERMANOVA_R2": np.nan,
            "complete_case_PERMANOVA_pseudo_F": np.nan,
            "complete_case_PERMDISP_p": np.nan,
            "complete_case_PERMDISP_F": np.nan
        }

        if category_results is not None:
            (
                category_stat_matrix,
                category_metadata,
                category_perm,
                category_disp,
                category_pairwise
            ) = category_results

            category_row.update(
                {
                    "n_complete_case_databases":
                        category_metadata.nunique(),
                    "complete_case_PERMANOVA_p":
                        category_perm["p_value"],
                    "complete_case_PERMANOVA_R2":
                        category_perm["R2"],
                    "complete_case_PERMANOVA_pseudo_F":
                        category_perm["pseudo_F"],
                    "complete_case_PERMDISP_p":
                        category_disp["p_value"],
                    "complete_case_PERMDISP_F":
                        category_disp["F"]
                }
            )

        category_summary_rows.append(category_row)

        print(f"Saved category membership: {membership_file}")
        print(f"Saved category distances: {distance_file}")
        print(f"Saved category matrix: {matrix_file}")

        if complete_matrix_file is not None:
            print(
                f"Saved complete-case matrix: "
                f"{complete_matrix_file}"
            )

    category_summary = pd.DataFrame(
        category_summary_rows
    )

    category_summary_file = (
        output_prefix
        + "_category_comparison_summary.tsv"
    )

    category_summary.to_csv(
        category_summary_file,
        sep="\t",
        index=False
    )

    print("\n============================================================")
    print("Category-specific analysis finished.")
    print("============================================================")
    print(f"Category summary: {category_summary_file}")

    return category_summary


# ============================================================
# Main
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description=(
            "Analyze BiG-SCAPE SQLite distances for Krill.\n\n"
            "Singleton definition:\n"
            "A user cand_cluster BGC is a singleton when it has "
            "no other cand_cluster with distance <= threshold.\n\n"
            "Same-GBK missing comparisons are handled by retaining "
            "the largest BGC as the GBK representative."
        )
    )

    parser.add_argument(
        "db_file",
        help="Path to BiGSCAPE.db"
    )

    parser.add_argument(
        "input_path",
        help=(
            "Krill input directory containing "
            "user BGCs"
        )
    )

    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="Output prefix"
    )

    parser.add_argument(
        "--threshold",
        type=float,
        default=0.3,
        help=(
            "Distance threshold for singleton "
            "classification. Default: 0.3"
        )
    )

    parser.add_argument(
        "--permutations",
        type=int,
        default=999,
        help=(
            "Number of permutations for "
            "PERMANOVA/PERMDISP. Default: 999"
        )
    )

    args = parser.parse_args()

    db_file = os.path.abspath(
        args.db_file
    )

    input_path = os.path.abspath(
        args.input_path
    )

    if args.output is None:

        output_dir = os.path.join(os.path.dirname(db_file), "Krill_BiGSCAPE")
        print("\n\n output_directory: ", output_dir)

        os.makedirs(output_dir, exist_ok=True)

        output_prefix = os.path.join(output_dir, "Krill_BiGSCAPE")
        print("\n\n output_prefix: ", output_prefix)

    else:

        output_dir = os.path.abspath(args.output)
        print("\n\n output_args: ", output_dir)

        output_prefix = os.path.join(output_dir, "Krill_BiGSCAPE")
        print("\n\n output_prefix: ", output_prefix)

        if output_dir:
            os.makedirs(output_dir,exist_ok=True)



    # ========================================================
    # 1. Read SQLite
    # ========================================================

    distance, gbk, bgc_record = read_tables(
        db_file
    )

    original_distance = distance.copy()

    # ========================================================
    # 2. Add GBK paths
    # ========================================================

    bgc_record = add_path_to_bgc_record(
        bgc_record,
        gbk
    )

    # ========================================================
    # 3. Annotate distance table
    # ========================================================

    distance = add_bgc_information_to_distance(
        distance,
        bgc_record
    )

    distance = add_distance_source(
        distance
    )

    # ========================================================
    # 4. Identify user cand_clusters
    # ========================================================

    user_candidates = get_user_candidate_clusters(
        bgc_record,
        input_path
    )

    if user_candidates.empty:

        raise ValueError(
            "No user cand_cluster BGCs were found "
            "under input_path."
        )

    # ========================================================
    # 5. SINGLETON CLASSIFICATION
    #
    # IMPORTANT:
    #
    # This happens BEFORE removing MIBiG.
    #
    # A user BGC is a singleton when it has no other
    # cand_cluster with distance <= threshold.
    #
    # MIBiG neighbors therefore count for this classification.
    # ========================================================

    status = identify_singletons(
        distance,
        user_candidates,
        threshold=args.threshold
    )

    status["database"] = (
        status["path"]
        .apply(
            lambda x:
            get_database_from_path(
                x,
                input_path
            )
        )
    )

    # Add singleton information to the complete
    # annotated distance table.
    distance = add_singleton_columns(
        distance,
        status
    )

    # ========================================================
    # 6. Remove MIBiG × MIBiG
    # ========================================================

    distance_filtered = remove_mibig_mibig(
        distance
    )

    # ========================================================
    # 7. Keep cand_cluster × cand_cluster
    # ========================================================

    distance_filtered = keep_cand_clusters(
        distance_filtered
    )

    # ========================================================
    # 8. Add database labels
    # ========================================================

    distance_filtered = add_database(
        input_path,
        distance_filtered
    )

    # ========================================================
    # 9. USER-ONLY observed distances
    #
    # Remove every comparison involving MIBiG.
    #
    # This is separate from singleton classification.
    # ========================================================

    observed_user = remove_all_mibig(
        distance_filtered
    )

    observed_user_file = (
        output_prefix +
        "_distance_observed_user_only.tsv"
    )

    observed_user.to_csv(
        observed_user_file,
        sep="\t",
        index=False
    )

    print(
        f"\nUser-only distance table: "
        f"{observed_user_file}"
    )

    # ========================================================
    # 10. Missing user-user pairs
    # ========================================================

    missing_pairs = identify_missing_user_pairs(
        observed_user,
        user_candidates
    )

    missing_pairs_file = (
        output_prefix +
        "_missing_user_user_pairs.tsv"
    )

    missing_pairs.to_csv(
        missing_pairs_file,
        sep="\t",
        index=False
    )

    print(
        f"Missing-pair table: "
        f"{missing_pairs_file}"
    )

    # ========================================================
    # 11. Select one representative per affected GBK
    #
    # BiG-SCAPE intentionally does not calculate distances between
    # BGCs from the same GBK. Therefore same-GBK missing pairs do
    # not justify removing every BGC involved. Instead, for each
    # affected GBK, retain the largest BGC (End - Start) and use
    # that BGC as the representative for matrix/statistical analysis.
    #
    # Singleton classification above remains unchanged and still
    # uses the complete BiG-SCAPE distance table.
    # ========================================================

    matrix_user_ids, same_gbk_excluded = (
        select_gbk_representatives_from_missing_pairs(
            missing_pairs,
            user_candidates
        )
    )

    same_gbk_excluded_file = (
        output_prefix +
        "_same_gbk_smaller_BGCs_excluded.tsv"
    )

    same_gbk_excluded.to_csv(
        same_gbk_excluded_file,
        sep="\t",
        index=False
    )

    # Candidate table used by matrix/category analyses.
    representative_candidates = user_candidates[
        pd.to_numeric(
            user_candidates["id"],
            errors="coerce"
        ).isin(matrix_user_ids)
    ].copy()

    # ========================================================
    # 12. Complete-case observed dataset
    #
    # The complete-case diagnostic is now evaluated only after
    # the one-representative-per-affected-GBK rule. This prevents
    # intentionally undefined same-GBK comparisons from causing
    # all records from an affected GBK to be discarded.
    # ========================================================

    observed_complete, complete_ids = (
        subset_complete_case_distance(
            observed_user,
            representative_candidates
        )
    )


    excluded_ids = sorted(
        set(matrix_user_ids) -
        set(complete_ids)
    )

    excluded_file = (
        output_prefix +
        "_observed_complete_case_excluded_BGCs.tsv"
    )

    pd.DataFrame(
        {
            "record_id":
                excluded_ids
        }
    ).to_csv(
        excluded_file,
        sep="\t",
        index=False
    )

    print(
        f"\nComplete-case BGCs: "
        f"{len(complete_ids)}"
    )

    print(
        f"Excluded BGCs:      "
        f"{len(excluded_ids)}"
    )

    # ========================================================
    # 12. Save complete-case observed distances
    # ========================================================

    observed_complete_file = (
        output_prefix +
        "_distance_observed_complete_case.tsv"
    )

    observed_complete.to_csv(
        observed_complete_file,
        sep="\t",
        index=False
    )

    # ========================================================
    # 13. Save core tables
    # ========================================================

    distance_annotated_file = (
        output_prefix +
        "_distance_annotated.tsv"
    )

    distance_filtered_file = (
        output_prefix +
        "_distance_cand_cluster.tsv"
    )

    bgc_record_file = (
        output_prefix +
        "_bgc_record.tsv"
    )

    status_file = (
        output_prefix +
        "_singleton_classification.tsv"
    )

    distance.to_csv(
        distance_annotated_file,
        sep="\t",
        index=False
    )

    distance_filtered.to_csv(
        distance_filtered_file,
        sep="\t",
        index=False
    )

    bgc_record.to_csv(
        bgc_record_file,
        sep="\t",
        index=False
    )

    status.to_csv(
        status_file,
        sep="\t",
        index=False
    )

    # ========================================================
    # 14. Summary
    # ========================================================

    possible_pairs = (
        len(user_candidates)
        *
        (
            len(user_candidates) -
            1
        )
        //
        2
    )

    missing_count = len(
        missing_pairs
    )

    observed_pairs = (
        possible_pairs -
        missing_count
    )

    summary = {
        "threshold": args.threshold,
        "original_distance_rows": len(original_distance),
        "annotated_distance_rows": len(distance),
        "cand_cluster_distance_rows": len(distance_filtered),
        "user_cand_clusters": len(user_candidates),
        "singleton_BGCs": int(status["singleton"].sum()),
        "non_singleton_BGCs": int((~status["singleton"]).sum()),
        "possible_user_user_pairs": possible_pairs,
        "observed_user_user_pairs": observed_pairs,
        "missing_user_user_pairs": missing_count,
        "matrix_representative_BGCs": len(matrix_user_ids),
        "same_gbk_smaller_BGCs_excluded": len(same_gbk_excluded),
        "observed_complete_case_BGCs": len(complete_ids),
        "observed_excluded_BGCs": len(excluded_ids)
    }

    summary_file = (
        output_prefix +
        "_analysis_summary.tsv"
    )

    pd.DataFrame(
        [summary]
    ).to_csv(
        summary_file,
        sep="\t",
        index=False
    )

    # ========================================================
    # 15. ADDITIONAL BGC CATEGORY-SPECIFIC ANALYSES
    #
    # This is additive and does NOT replace the original
    # complete-case dataset comparison above.
    #
    # Each requested broad category receives its own matrix
    # and its own dataset-comparison statistics.
    # ========================================================

    category_results_summary = (
        run_category_comparisons(
            observed_user,
            user_candidates,
            output_prefix,
            permutations=args.permutations,
            representative_ids=matrix_user_ids
        )
    )

    # ========================================================
    # 16. Final report
    # ========================================================

    print("\n==============================================")
    print("Finished successfully")
    print("==============================================")

    print(
        f"Singleton threshold: "
        f"distance <= {args.threshold}"
    )

    print(
        f"Singletons:          "
        f"{int(status['singleton'].sum())}"
    )

    print(
        f"Non-singletons:      "
        f"{int((~status['singleton']).sum())}"
    )

    print(
        f"User BGCs:           "
        f"{len(user_candidates)}"
    )

    print(
        f"Possible user pairs: "
        f"{possible_pairs}"
    )

    print(
        f"Observed user pairs: "
        f"{observed_pairs}"
    )

    print(
        f"Missing user pairs:  "
        f"{missing_count}"
    )

    print(
        f"BGC representatives retained for matrix: "
        f"{len(matrix_user_ids)}"
    )

    print(
        f"Smaller same-GBK BGCs excluded: "
        f"{len(same_gbk_excluded)}"
    )

    print(
        f"Complete-case BGCs for classical methods: "
        f"{len(complete_ids)}"
    )

    print(
        f"BGCs excluded from complete-case subset: "
        f"{len(excluded_ids)}"
    )

    print(
        "\nBiGSCAPE.db was NOT modified."
    )

    print(
        "\nMain output:"
    )

    print(
        f"  {status_file}"
    )

    print(
        f"  {missing_pairs_file}"
    )

    print(
        f"  {summary_file}"
    )

    category_summary_file = (
        output_prefix +
        "_category_comparison_summary.tsv"
    )

    print(
        f"  {category_summary_file}"
    )


# ============================================================
# Run
# ============================================================

if __name__ == "__main__":

    main()