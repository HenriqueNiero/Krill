#!/usr/bin/env python3
"""
Krill BiG-SCAPE integration.

This module combines the three post-processing steps previously performed by
separate scripts:

1. Extract cand_cluster coordinates from antiSMASH region GBKs.
2. Add those coordinates to the BiG-SCAPE record/clustering tables and create
   BiGSCAPE/big_scape_extraction.tsv.
3. Integrate BiG-SCAPE GCF/CC/class/category information into
   DBsReportOutput/DBs_BGCs_with_Hits.tsv.

The coordinate convention intentionally follows Krill's AntiSMASH_extractor:
Biopython feature.location.start is kept as-is (0-based) and
feature.location.end is kept as-is (end-exclusive).
"""

import glob
import os
import re

import pandas as pd
from Bio import SeqIO


def _clean(value):
    """Return a stripped string, or an empty string for NaN/None."""
    if value is None or pd.isna(value):
        return ""
    return str(value).strip()


def _normalise_gbk_name(value):
    """
    Normalize a BiG-SCAPE GBK identifier to the original antiSMASH region
    filename stem.

    Examples:
        000000001_27.region001
        000000001_27.region001.gbk
        000000001_27.region001.gbk_region_1
    -> 000000001_27.region001
    """
    value = _clean(value)
    value = re.sub(r"_gbk_region_\d+$", "", value)
    value = re.sub(r"\.gbk$", "", value)
    return os.path.basename(value)


def _contig_from_gbk(value):
    """
    Derive the Krill contig name from an antiSMASH region GBK name.

    Example:
        000000001_27.region001.gbk -> 000000001_27
    """
    value = _normalise_gbk_name(value)
    value = re.sub(r"\.region\d+$", "", value)
    return value


def _extract_region_number(value):
    """Extract antiSMASH region number from a GBK identifier."""
    match = re.search(r"\.region(\d+)", _clean(value))
    return int(match.group(1)) if match else None


def _extract_original_coordinates(record):
    """
    Read antiSMASH original coordinates from the structured comment.

    Returns:
        (Orig_Start, Orig_End)

    The values are None when antiSMASH original-coordinate metadata is absent.
    """
    structured_comment = record.annotations.get("structured_comment", {})
    antismash_data = structured_comment.get("antiSMASH-Data", {})

    orig_start = antismash_data.get("Orig. start")
    orig_end = antismash_data.get("Orig. end")

    try:
        orig_start = int(orig_start) if orig_start is not None else None
    except (TypeError, ValueError):
        orig_start = None

    try:
        orig_end = int(orig_end) if orig_end is not None else None
    except (TypeError, ValueError):
        orig_end = None

    return orig_start, orig_end


def _get_candidate_number(feature):
    """
    Extract a candidate-cluster number.

    The normal antiSMASH qualifier is candidate_cluster_number. The fallback
    qualifiers make the parser tolerant of older antiSMASH files.
    """
    qualifiers = feature.qualifiers

    for key in (
        "candidate_cluster_number",
        "protocluster_number",
        "number",
        "label",
        "note",
    ):
        values = qualifiers.get(key)
        if not values:
            continue

        if isinstance(values, list):
            value = " ".join(map(str, values))
        else:
            value = str(values)

        match = re.search(r"\d+", value)
        if match:
            return int(match.group())

    return None


def _find_bigscape_run(path, bigscape_dir, cutoff):
    """Locate the newest BiG-SCAPE2 run for the requested cutoff."""
    output_dir = os.path.join(bigscape_dir, "output_files")

    if not os.path.isdir(output_dir):
        raise FileNotFoundError(
            "BiG-SCAPE output directory not found:\n{}".format(output_dir)
        )

    pattern = os.path.join(output_dir, "*c{}".format(cutoff))
    run_dirs = sorted(
        d for d in glob.glob(pattern) if os.path.isdir(d)
    )

    if not run_dirs:
        # Be slightly more tolerant of 0.3 vs 0.30 naming.
        all_runs = sorted(
            d for d in glob.glob(os.path.join(output_dir, "*"))
            if os.path.isdir(d)
        )
        run_dirs = [
            d for d in all_runs
            if re.search(r"_c0*{}$".format(re.escape(str(cutoff).replace(".", r"\."))), os.path.basename(d))
        ]

    if not run_dirs:
        raise FileNotFoundError(
            "No BiG-SCAPE run found for cutoff {} in:\n{}".format(
                cutoff, output_dir
            )
        )

    run_dir = run_dirs[-1]

    record_file = os.path.join(run_dir, "record_annotations.tsv")
    mix_file = os.path.join(
        run_dir, "mix", "mix_clustering_c{}.tsv".format(cutoff)
    )

    if not os.path.isfile(record_file):
        raise FileNotFoundError(
            "BiG-SCAPE record_annotations.tsv not found:\n{}".format(
                record_file
            )
        )

    if not os.path.isfile(mix_file):
        # Fall back to the only mix clustering file if the cutoff was written
        # as 0.30 rather than 0.3.
        candidates = glob.glob(
            os.path.join(run_dir, "mix", "mix_clustering_c*.tsv")
        )
        if len(candidates) == 1:
            mix_file = candidates[0]
        else:
            raise FileNotFoundError(
                "BiG-SCAPE mix clustering file not found:\n{}".format(
                    mix_file
                )
            )

    return run_dir, record_file, mix_file


def extract_candidate_clusters_coordinates(input_dir, output_file=None):
    """
    Extract antiSMASH cand_cluster coordinates from region GBKs.

    Only *.region*.gbk files are processed. BiG-SCAPE output is skipped.

    For antiSMASH region files, candidate-cluster coordinates are converted
    back to original-record coordinates using:

        Start = Orig. start + candidate_relative_start
        End   = Orig. start + candidate_relative_end

    Coordinates use the same convention as Krill's existing extractor:
    0-based start and end-exclusive end.

    Returns a pandas DataFrame.
    """
    rows = []

    for root, dirs, files in os.walk(input_dir):
        dirs[:] = [
            d for d in dirs
            if d != "BiGSCAPE"
        ]

        for filename in files:
            if not filename.endswith(".gbk"):
                continue

            if ".region" not in filename:
                continue

            gbk_path = os.path.join(root, filename)

            try:
                for record in SeqIO.parse(gbk_path, "genbank"):
                    orig_start, orig_end = _extract_original_coordinates(record)

                    if orig_start is None or orig_end is None:
                        print(
                            "WARNING: Missing antiSMASH original coordinates: "
                            "{}".format(gbk_path)
                        )
                        continue

                    for feature in record.features:
                        if feature.type != "cand_cluster":
                            continue

                        candidate_number = _get_candidate_number(feature)

                        if candidate_number is None:
                            print(
                                "WARNING: cand_cluster without a candidate "
                                "cluster number: {}".format(gbk_path)
                            )
                            continue

                        candidate_start = int(feature.location.start)
                        candidate_end = int(feature.location.end)

                        original_start = orig_start + candidate_start
                        original_end = orig_start + candidate_end

                        rows.append(
                            {
                                "GBK": _normalise_gbk_name(filename),
                                "GBK_File": filename,
                                "GBK_Path": gbk_path,
                                "Record_ID": record.id,
                                "Record_Type": "cand_cluster",
                                "Record_Number": int(candidate_number),
                                "candidate_cluster_number": int(candidate_number),
                                "Orig_Start": orig_start,
                                "Orig_End": orig_end,
                                "Region_Start": candidate_start,
                                "Region_End": candidate_end,
                                "Start": original_start,
                                "End": original_end,
                                "contig": _contig_from_gbk(filename),
                                "product": ",".join(
                                    feature.qualifiers.get("product", [])
                                ),
                                "kind": ",".join(
                                    feature.qualifiers.get("kind", [])
                                ),
                                "protoclusters": ",".join(
                                    feature.qualifiers.get("protoclusters", [])
                                ),
                            }
                        )

            except Exception as exc:
                print(
                    "WARNING: Could not process {}: {}".format(
                        gbk_path, exc
                    )
                )

    coordinates = pd.DataFrame(rows)

    if coordinates.empty:
        raise RuntimeError(
            "No antiSMASH cand_cluster coordinates were extracted from:\n{}".format(
                input_dir
            )
        )

    coordinates = coordinates.drop_duplicates(
        subset=["GBK", "Record_Type", "Record_Number", "Start", "End"],
        keep="first",
    )

    coordinates["Record_Number"] = pd.to_numeric(
        coordinates["Record_Number"],
        errors="coerce",
    ).astype("Int64")

    coordinates = coordinates.sort_values(
        ["GBK", "Record_Number"]
    ).reset_index(drop=True)

    if output_file:
        os.makedirs(os.path.dirname(os.path.abspath(output_file)), exist_ok=True)
        coordinates.to_csv(
            output_file,
            sep="\t",
            index=False,
        )
        print(
            "Candidate-cluster coordinate table written to:\n{}".format(
                output_file
            )
        )

    print(
        "Candidate clusters extracted: {}".format(len(coordinates))
    )

    return coordinates


def build_bigscape_extraction_table(
    path,
    bigscape_dir=None,
    cutoff="0.3",
    coordinates=None,
):
    """
    Build the coordinate-aware BiG-SCAPE extraction table.

    Reads:
        BiGSCAPE/output_files/<run>/record_annotations.tsv
        BiGSCAPE/output_files/<run>/mix/mix_clustering_c<cutoff>.tsv

    Adds antiSMASH cand_cluster coordinates and creates:
        BiGSCAPE/big_scape_extraction.tsv

    Returns the path to the generated extraction table.
    """
    if bigscape_dir is None:
        bigscape_dir = os.path.join(path, "BiGSCAPE")

    if coordinates is None:
        coordinate_file = os.path.join(
            bigscape_dir,
            "candidate_clusters_coordinates.tsv",
        )
        coordinates = extract_candidate_clusters_coordinates(
            path,
            output_file=coordinate_file,
        )

    run_dir, record_file, mix_file = _find_bigscape_run(
        path,
        bigscape_dir,
        cutoff,
    )

    print("\n============================================")
    print("Building BiG-SCAPE extraction table")
    print("============================================")
    print("BiG-SCAPE run: {}".format(run_dir))
    print("Record annotations: {}".format(record_file))
    print("Mix clustering: {}".format(mix_file))

    records = pd.read_csv(
        record_file,
        sep="\t",
        dtype=str,
    )

    mix = pd.read_csv(
        mix_file,
        sep="\t",
        dtype=str,
    )

    required_record_columns = [
        "Record",
        "GBK",
        "Record_Type",
        "Record_Number",
        "Class",
        "Category",
    ]

    required_mix_columns = [
        "Record",
        "GBK",
        "Record_Type",
        "Record_Number",
        "CC",
        "Family",
    ]

    for column in required_record_columns:
        if column not in records.columns:
            raise ValueError(
                "Column '{}' missing from {}".format(
                    column, record_file
                )
            )

    for column in required_mix_columns:
        if column not in mix.columns:
            raise ValueError(
                "Column '{}' missing from {}".format(
                    column, mix_file
                )
            )

    # Keep all BiG-SCAPE annotation information, but normalize the keys.
    records = records.copy()
    mix = mix.copy()

    for df in (records, mix):
        df["GBK"] = df["GBK"].map(_normalise_gbk_name)
        df["Record_Type"] = df["Record_Type"].astype(str).str.strip()
        df["Record_Number"] = pd.to_numeric(
            df["Record_Number"],
            errors="coerce",
        ).astype("Int64")

    coordinates = coordinates.copy()
    coordinates["GBK"] = coordinates["GBK"].map(_normalise_gbk_name)
    coordinates["Record_Type"] = (
        coordinates["Record_Type"].astype(str).str.strip()
    )
    coordinates["Record_Number"] = pd.to_numeric(
        coordinates["Record_Number"],
        errors="coerce",
    ).astype("Int64")

    coordinate_columns = [
        "GBK",
        "Record_Type",
        "Record_Number",
        "contig",
        "Orig_Start",
        "Orig_End",
        "Region_Start",
        "Region_End",
        "Start",
        "End",
        "candidate_cluster_number",
        "GBK_File",
        "GBK_Path",
        "Record_ID",
    ]

    coordinate_columns = [
        c for c in coordinate_columns
        if c in coordinates.columns
    ]

    coordinate_key = [
        "GBK",
        "Record_Type",
        "Record_Number",
    ]

    coordinates_small = coordinates[coordinate_columns].copy()

    # The coordinate key must identify exactly one cand_cluster.
    duplicate_coordinates = coordinates_small.duplicated(
        subset=coordinate_key,
        keep=False,
    )

    if duplicate_coordinates.any():
        duplicated = coordinates_small.loc[
            duplicate_coordinates,
            coordinate_key,
        ].drop_duplicates()

        raise ValueError(
            "Duplicate cand_cluster coordinate keys found:\n{}".format(
                duplicated.to_string(index=False)
            )
        )

    records = records.merge(
        coordinates_small,
        on=coordinate_key,
        how="left",
    )

    mix = mix.merge(
        coordinates_small,
        on=coordinate_key,
        how="left",
    )

    # Merge annotations and clustering by BiG-SCAPE Record.
    annotations_small = records[
        [
            "Record",
            "Class",
            "Category",
        ]
    ].drop_duplicates(subset=["Record"])

    bigscape = mix.merge(
        annotations_small,
        on="Record",
        how="left",
    )

    # Keep names used by the final Krill integration.
    bigscape["contig"] = bigscape["contig"].fillna(
        bigscape["GBK"].map(_contig_from_gbk)
    )

    # Coordinates in the final table deliberately use the same names as
    # Krill's BGC table so the final merge is coordinate-based.
    bigscape["Start"] = pd.to_numeric(bigscape["Start"], errors="coerce").astype("Int64")
    bigscape["End"] = pd.to_numeric(bigscape["End"], errors="coerce").astype("Int64")

    # Remove possible duplicate rows before the final merge.
    bigscape = bigscape.drop_duplicates(
        subset=["contig", "Start", "End"],
        keep="first",
    ).reset_index(drop=True)

    output_file = os.path.join(
        bigscape_dir,
        "big_scape_extraction.tsv",
    )

    os.makedirs(bigscape_dir, exist_ok=True)

    bigscape.to_csv(
        output_file,
        sep="\t",
        index=False,
    )

    print(
        "BiG-SCAPE extraction table written to:\n{}".format(
            output_file
        )
    )

    return output_file


def integrate_bigscape_with_krill(
    path,
    bigscape_dir=None,
    cutoff="0.3",
):
    """
    Integrate coordinate-aware BiG-SCAPE results with Krill.

    Input:
        DBsReportOutput/DBs_BGCs_with_Hits.tsv

    Output:
        DBsReportOutput/DBs_BGCs_with_Hits_BiGSCAPE.tsv
        DBsReportOutput/DBs_BGCs_BiGSCAPE_unmatched.tsv
        BiGSCAPE/candidate_clusters_coordinates.tsv
        BiGSCAPE/big_scape_extraction.tsv

    The original Krill table is never overwritten.

    The final join uses:
        contig + Start + End

    'kind' is intentionally NOT part of the key because it is an annotation,
    not genomic identity.
    """
    if bigscape_dir is None:
        bigscape_dir = os.path.join(path, "BiGSCAPE")

    krill_file = os.path.join(
        path,
        "DBsReportOutput",
        "DBs_BGCs_with_Hits.tsv",
    )

    output_file = os.path.join(
        path,
        "DBsReportOutput",
        "DBs_BGCs_with_Hits_BiGSCAPE.tsv",
    )

    unmatched_file = os.path.join(
        path,
        "DBsReportOutput",
        "DBs_BGCs_BiGSCAPE_unmatched.tsv",
    )

    if not os.path.isfile(krill_file):
        raise FileNotFoundError(
            "Krill result table not found:\n{}".format(krill_file)
        )

    print("\n############################################")
    print("# Integrating BiG-SCAPE with Krill")
    print("############################################")

    # Step 1 + Step 2.
    bigscape_file = build_bigscape_extraction_table(
        path=path,
        bigscape_dir=bigscape_dir,
        cutoff=cutoff,
    )

    krill = pd.read_csv(
        krill_file,
        sep="\t",
        dtype=str,
    )

    bigscape = pd.read_csv(
        bigscape_file,
        sep="\t",
        dtype=str,
    )

    merge_keys = [
        "contig",
        "Start",
        "End",
    ]

    for key in merge_keys:
        if key not in krill.columns:
            raise ValueError(
                "Required Krill column '{}' not found in {}".format(
                    key, krill_file
                )
            )

        if key not in bigscape.columns:
            raise ValueError(
                "Required BiG-SCAPE column '{}' not found in {}".format(
                    key, bigscape_file
                )
            )

        krill[key] = krill[key].astype(str).str.strip()
        bigscape[key] = bigscape[key].astype(str).str.strip()

    # Select only information that should be added to Krill.
    bigscape_columns = [
        "Record",
        "Family",
        "CC",
        "Class",
        "Category",
    ]

    for column in bigscape_columns:
        if column not in bigscape.columns:
            raise ValueError(
                "Required BiG-SCAPE column '{}' not found in {}".format(
                    column, bigscape_file
                )
            )

    bigscape_sub = bigscape[
        merge_keys + bigscape_columns
    ].copy()

    bigscape_sub.rename(
        columns={
            "Record": "BiGSCAPE_Record",
            "Family": "BiGSCAPE_GCF",
            "CC": "BiGSCAPE_CC",
            "Class": "BiGSCAPE_Class",
            "Category": "BiGSCAPE_Category",
        },
        inplace=True,
    )

    # A genomic coordinate should identify one candidate cluster.
    duplicate_key = bigscape_sub.duplicated(
        subset=merge_keys,
        keep=False,
    )

    if duplicate_key.any():
        duplicated = bigscape_sub.loc[
            duplicate_key,
            merge_keys,
        ].drop_duplicates()

        raise ValueError(
            "Duplicate BiG-SCAPE records found for the same genomic "
            "coordinate key:\n{}".format(
                duplicated.to_string(index=False)
            )
        )

    # Preserve the original Krill row count.
    merged = krill.merge(
        bigscape_sub,
        on=merge_keys,
        how="left",
        validate="many_to_one",
    )

    added_columns = [
        "BiGSCAPE_Record",
        "BiGSCAPE_GCF",
        "BiGSCAPE_CC",
        "BiGSCAPE_Class",
        "BiGSCAPE_Category",
    ]

    for column in added_columns:
        merged[column] = merged[column].fillna("")

    # Singleton status is based on the complete BiG-SCAPE extraction table.
    valid_gcf = (
        bigscape["Family"]
        .replace("", pd.NA)
        .dropna()
        .astype(str)
        .str.strip()
    )

    gcf_counts = valid_gcf.value_counts()

    def singleton_status(gcf):
        gcf = _clean(gcf)
        if not gcf:
            return ""
        return (
            "True"
            if gcf_counts.get(gcf, 0) == 1
            else "False"
        )

    merged["BiGSCAPE_Singleton"] = (
        merged["BiGSCAPE_GCF"].map(singleton_status)
    )

    # Save integrated result.
    os.makedirs(
        os.path.dirname(output_file),
        exist_ok=True,
    )

    merged.to_csv(
        output_file,
        sep="\t",
        index=False,
    )

    # Save unmatched rows.
    unmatched_mask = merged["BiGSCAPE_Record"].eq("")

    unmatched_columns = [
        c for c in [
            "Database",
            "OriginalName",
            "Original_contig",
            "contig",
            "cluster_number",
            "Start",
            "End",
            "kind",
        ]
        if c in merged.columns
    ]

    merged.loc[
        unmatched_mask,
        unmatched_columns,
    ].to_csv(
        unmatched_file,
        sep="\t",
        index=False,
    )

    matched = int((~unmatched_mask).sum())
    unmatched = int(unmatched_mask.sum())
    total = len(merged)

    print("\n============================================")
    print("BiG-SCAPE / Krill integration summary")
    print("============================================")
    print("Krill BGCs       : {}".format(total))
    print("Matched          : {}".format(matched))
    print("Unmatched        : {}".format(unmatched))

    if total:
        print(
            "Match percentage : {:.2f}%".format(
                matched * 100.0 / total
            )
        )

    print("Integrated table : {}".format(output_file))
    print("Unmatched table  : {}".format(unmatched_file))
    print("============================================\n")

    return output_file


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Integrate BiG-SCAPE results into Krill."
    )

    parser.add_argument(
        "--path",
        required=True,
        help="Krill database/run directory.",
    )

    parser.add_argument(
        "--cutoff",
        default="0.3",
        help="BiG-SCAPE GCF cutoff. Default: 0.3",
    )

    parser.add_argument(
        "--bigscape-dir",
        default=None,
        help="BiG-SCAPE directory. Default: <path>/BiGSCAPE",
    )

    args = parser.parse_args()

    integrate_bigscape_with_krill(
        path=args.path,
        bigscape_dir=args.bigscape_dir,
        cutoff=args.cutoff,
    )
