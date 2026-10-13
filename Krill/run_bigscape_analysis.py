#!/usr/bin/env python3

import os
import glob
import subprocess


def find_bigscape_database(path):
    """
    Find the most recent BiGSCAPE.db generated under the
    Krill BiGSCAPE output directory.

    Expected structure:

        PATH/
        └── BiGSCAPE/
            └── YYYY-MM-DD_HH-MM-SS/
                └── BiGSCAPE.db

    Returns
    -------
    str
        Absolute path to the newest BiGSCAPE.db.
    """

    path = os.path.abspath(path)

    bigscape_root = os.path.join(
        path,
        "BiGSCAPE"
    )

    if not os.path.isdir(bigscape_root):
        raise FileNotFoundError(
            "\nBiG-SCAPE output directory was not found:\n"
            f"{bigscape_root}\n\n"
            "Run Krill with --bigscape first."
        )

    pattern = os.path.join(
        bigscape_root,
        "**",
        "BiGSCAPE.db"
    )

    databases = glob.glob(
        pattern,
        recursive=True
    )

    if not databases:
        raise FileNotFoundError(
            "\nNo BiGSCAPE.db was found under:\n"
            f"{bigscape_root}\n\n"
            "Run Krill with --bigscape first."
        )

    # Use modification time so that the newest BiG-SCAPE
    # result is selected.
    databases.sort(
        key=os.path.getmtime,
        reverse=True
    )

    database = os.path.abspath(
        databases[0]
    )

    print(
        "\n=============================================="
    )
    print(
        "BiG-SCAPE database selected for analysis"
    )
    print(
        "=============================================="
    )
    print(
        f"{database}"
    )
    print(
        "==============================================\n"
    )

    return database


def run(
    path,
    env="bigscape",
    threshold=0.3,
    permutations=999,
    output=None
):
    """
    Run Krill_BiGSCAPE_category_analysis.py on the
    newest BiGSCAPE.db found under PATH/BiGSCAPE.

    Parameters
    ----------
    path : str
        Krill input/project directory.

    env : str
        Conda environment used for the analysis.

    threshold : float
        Singleton distance threshold.

    permutations : int
        Number of permutations for PERMANOVA/PERMDISP.

    output : str or None
        Optional output directory.
    """

    project_dir = os.path.dirname(
        os.path.abspath(__file__)
    )

    analysis_script = os.path.join(
        project_dir,
        "Krill_BiGSCAPE_category_analysis.py"
    )

    if not os.path.isfile(analysis_script):
        raise FileNotFoundError(
            "\nKrill_BiGSCAPE_category_analysis.py "
            "was not found:\n"
            f"{analysis_script}"
        )

    path = os.path.abspath(path)

    db_file = find_bigscape_database(
        path
    )

    command = [
        "conda",
        "run",
        "-n",
        env,
        "python",
        analysis_script,
        db_file,
        path,
        "--threshold",
        str(threshold),
        "--permutations",
        str(permutations)
    ]

    if output is not None:
        command.extend(
            [
                "--output",
                os.path.abspath(output)
            ]
        )

    print(
        "\n=============================================="
    )
    print(
        "Running Krill BiG-SCAPE category analysis"
    )
    print(
        "=============================================="
    )

    print(
        "Conda environment:",
        env
    )

    print(
        "Threshold:",
        threshold
    )

    print(
        "Permutations:",
        permutations
    )

    print(
        "==============================================\n"
    )

    try:

        subprocess.run(
            command,
            check=True
        )

    except subprocess.CalledProcessError as error:

        raise RuntimeError(
            "\nBiG-SCAPE category analysis failed."
            f"\nExit code: {error.returncode}"
        )

    print(
        "\n=============================================="
    )
    print(
        "BiG-SCAPE category analysis completed"
    )
    print(
        "=============================================="
    )

    return db_file


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(
        description=(
            "Run Krill BiG-SCAPE category analysis "
            "on an existing BiGSCAPE.db."
        )
    )

    parser.add_argument(
        "path",
        help=(
            "Krill input/project directory "
            "containing the BiGSCAPE directory."
        )
    )

    parser.add_argument(
        "--env",
        default="bigscape",
        help=(
            "Conda environment used to run "
            "the analysis. Default: bigscape"
        )
    )

    parser.add_argument(
        "--threshold",
        type=float,
        default=0.3,
        help=(
            "Singleton distance threshold. "
            "Default: 0.3"
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

    parser.add_argument(
        "--output",
        default=None,
        help=(
            "Optional output directory."
        )
    )

    args = parser.parse_args()

    run(
        path=args.path,
        env=args.env,
        threshold=args.threshold,
        permutations=args.permutations,
        output=args.output
    )