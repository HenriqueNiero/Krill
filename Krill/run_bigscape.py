import os
import subprocess


def run(path, threads=8, cutoff="0.30", env="bigscape", pfam=None):

    project_dir = os.path.dirname(os.path.abspath(__file__))
    print("\n ################## run_bigscape.py #################")
    print("\n> project_dir: ", project_dir)
    print("\n> db: ", path)
    input_dir = os.path.abspath(path)   
    output_dir = os.path.join(input_dir, "BiGSCAPE")

    os.makedirs(output_dir, exist_ok=True)

    # Skip if BiG-SCAPE has already finished
    flag = os.path.join(output_dir, "completed")

    if os.path.exists(flag):
        print("BiG-SCAPE already completed.")

        return output_dir


    try:
        subprocess.run(
            [
                "bash",
                os.path.join(project_dir, "krill_run_bigscape.sh"),
                env,
                input_dir,
                output_dir,
                str(threads),
                str(cutoff),
                pfam if pfam is not None else ""
            ],
            check=True
        )


    except subprocess.CalledProcessError as error:

        raise RuntimeError(
            f"BiG-SCAPE failed with exit code "
            f"{error.returncode}"
        )

    open(flag, "w").close()

    print("\nBiG-SCAPE completed successfully.")

    return output_dir