import os
import tempfile
import subprocess
import sys

from dnabyte.cluster import Cluster


class CloverClusterer(Cluster):
    """
    DNAbyte wrapper around the actual dna-clover installation.

    Input:
        list[str]
            DNA sequences.

    Output:
        dict[int, list[str]]
            DNAbyte ClusteredDNA-compatible clusters.

    Clover itself performs the clustering. This class only:

        1. converts DNAbyte sequences to Clover input format,
        2. runs the installed Clover program,
        3. parses Clover's output,
        4. converts it back to DNAbyte format.
    """

    def __init__(self, params, logger=None):
        self.params = params
        self.logger = logger


        self.read_len = getattr(params, "clover_read_len", 200)

        # Actual Clover parameters.
        self.horizontal_drift = getattr(
            params, "clover_horizontal_drift", 3
        )
        self.vertical_drift = getattr(
            params, "clover_vertical_drift", 2
        )
        self.tree_threshold = getattr(
            params, "clover_tree_threshold", 10
        )
        self.now_clust_threshold = getattr(
            params, "clover_now_clust_threshold", 8
        )

        self.align = getattr(params, "clover_align", False)

        # Keep temporary files unless explicitly disabled.
        # self.keep_temp = getattr(params, "clover_keep_temp", False)
        
        self.keep_temp = True

    def cluster(self, sequences):
        """
        Cluster DNAbyte sequences using the real installed Clover package.
        """

        if not isinstance(sequences, list):
            raise TypeError("Clover expects a list of DNA sequences.")

        if not sequences:
            return {}, {
                "num_clusters": 0,
                "num_sequences": 0,
            }

        # Validate sequences before sending them to Clover.
        for sequence in sequences:
            if not isinstance(sequence, str):
                raise TypeError(
                    "Every sequence passed to Clover must be a string."
                )

            invalid = set(sequence.upper()) - {"A", "C", "G", "T"}
            if invalid:
                raise ValueError(
                    f"Invalid DNA bases found: {invalid}"
                )

        # Clover expects:
        #
        #     index sequence
        #
        # one read per line.
        input_lines = [
            f"{i} {sequence.upper()}"
            for i, sequence in enumerate(sequences)
        ]

        # temp_dir = tempfile.mkdtemp(prefix="dnabyte_clover_")
        temp_dir = os.path.abspath("clover_debug")
        os.makedirs(temp_dir, exist_ok=True)

        input_file = os.path.join(temp_dir, "input.txt")
        output_file = os.path.join(temp_dir, "clover_output.txt")

        # os.makedirs(output_dir, exist_ok=True)  # No longer needed, output is a file, not a directory

        try:
            with open(
                input_file,
                "w",
                encoding="utf-8"
            ) as f:
                f.write("\n".join(input_lines))
                f.write("\n")

            self._run_clover(
                input_file=input_file,
                output_file=output_file,
                read_len=self.read_len,
            )

            clusters = self._parse_output(
                output_file=output_file,
                sequences=sequences,
            )

            info = {
                "num_clusters": len(clusters),
                "num_sequences": sum(
                    len(reads) for reads in clusters.values()
                ),
                "clover_read_len": self.read_len,
                "clover_horizontal_drift": self.horizontal_drift,
                "clover_vertical_drift": self.vertical_drift,
                "clover_tree_threshold": self.tree_threshold,
                "clover_now_clust_threshold": (
                    self.now_clust_threshold
                ),
                "clover_align": self.align,
            }

            return clusters, info

        finally:
            if not self.keep_temp:
                import shutil
                shutil.rmtree(temp_dir, ignore_errors=True)

    def _run_clover(self, input_file, output_file, read_len):
        """
        Run the actual installed dna-clover package.

        Clover's -O argument is an OUTPUT FILE NAME, not a directory.
        """

        command = [
            sys.executable,
            "-m",
            "clover.main",
            "-I",
            input_file,
            "-O",
            output_file,
            "-L",
            str(read_len),
            "-P",
            "0",
            "--no-tag",
        ]

        if self.align:
            command.append("--align")

        if self.logger:
            self.logger.info(
                "Running Clover: %s",
                " ".join(command)
            )

        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
        )

        if self.logger:
            if result.stdout:
                self.logger.info(
                    "Clover stdout:\n%s",
                    result.stdout
                )

            if result.stderr:
                self.logger.info(
                    "Clover stderr:\n%s",
                    result.stderr
                )

        if result.returncode != 0:
            raise RuntimeError(
                f"Clover failed with exit code {result.returncode}\n\n"
                f"STDOUT:\n{result.stdout}\n\n"
                f"STDERR:\n{result.stderr}"
            )

        return result

    def _parse_output(self, output_file, sequences):
        """
        Convert Clover's output back into:

            {
                cluster_id: [sequence, sequence, ...]
            }

        Clover output formats can differ between versions, so we
        inspect the generated files rather than assuming that one
        hard-coded filename exists.
        """

        files = []

        for root, _, filenames in os.walk(os.path.dirname(output_file)):
            for filename in filenames:
                files.append(
                    os.path.join(root, filename)
                )

        if not files:
            raise RuntimeError(
                "Clover completed successfully but produced "
                "no output files."
            )

        # Map original read index -> sequence.
        sequence_by_index = {
            str(i): sequence
            for i, sequence in enumerate(sequences)
        }

        clusters = {}

        # Look through Clover's generated text files.
        for filepath in files:
            try:
                with open(
                    filepath,
                    "r",
                    encoding="utf-8",
                    errors="ignore",
                ) as f:
                    lines = f.readlines()
            except OSError:
                continue

            for line in lines:
                line = line.strip()

                if not line:
                    continue

                parts = line.split()

                # We expect Clover output to contain read/index
                # information. Do not blindly treat arbitrary
                # output/statistics files as clusters.
                if len(parts) < 2:
                    continue

                cluster_id = parts[0]

                read_indices = []

                for part in parts[1:]:
                    if part in sequence_by_index:
                        read_indices.append(part)

                if not read_indices:
                    continue

                cluster = clusters.setdefault(
                    cluster_id,
                    []
                )

                for index in read_indices:
                    sequence = sequence_by_index[index]

                    if sequence not in cluster:
                        cluster.append(sequence)

        if not clusters:
            raise RuntimeError(
                "Clover produced output, but the wrapper could not "
                "parse any clusters.\n\n"
                "Generated files:\n"
                + "\n".join(files)
            )

        # Convert Clover's IDs to DNAbyte integer IDs.
        result = {}

        for new_id, (_, reads) in enumerate(
            clusters.items()
        ):
            result[new_id] = reads

        return result


def check_parameter(
    parameter,
    default,
    minimum,
    maximum,
    inputparams,
):
    if (
        not hasattr(inputparams, parameter)
        or inputparams.__dict__.get(parameter) is None
    ):
        return default

    value = inputparams.__dict__[parameter]

    if not minimum <= value <= maximum:
        raise ValueError(
            f"{parameter} must be greater than or equal to "
            f"{minimum} and less than or equal to {maximum}, "
            f"got {value}"
        )

    return value


def attributes(params):
    """
    DNAbyte plugin configuration.

    These names are prefixed with clover_ so they don't collide
    with parameters belonging to other clustering methods.
    """

    clover_read_len = check_parameter(
        "clover_read_len",
        200,
        1,
        10000,
        params,
    )

    clover_horizontal_drift = check_parameter(
        "clover_horizontal_drift",
        3,
        0,
        clover_read_len,
        params,
    )

    clover_vertical_drift = check_parameter(
        "clover_vertical_drift",
        2,
        0,
        clover_read_len,
        params,
    )

    clover_tree_threshold = check_parameter(
        "clover_tree_threshold",
        10,
        0,
        clover_read_len,
        params,
    )

    clover_now_clust_threshold = check_parameter(
        "clover_now_clust_threshold",
        8,
        0,
        clover_read_len,
        params,
    )

    clover_align = getattr(
        params,
        "clover_align",
        False,
    )

    clover_keep_temp = getattr(
        params,
        "clover_keep_temp",
        False,
    )

    return {
        "clover_read_len": clover_read_len,
        "clover_horizontal_drift": clover_horizontal_drift,
        "clover_vertical_drift": clover_vertical_drift,
        "clover_tree_threshold": clover_tree_threshold,
        "clover_now_clust_threshold": (
            clover_now_clust_threshold
        ),
        "clover_align": clover_align,
        "clover_keep_temp": clover_keep_temp,
    }