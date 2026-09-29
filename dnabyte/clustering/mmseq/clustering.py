import os
import subprocess
import tempfile
from collections import defaultdict

from dnabyte.cluster import Cluster


class MMseqsClusterer(Cluster):

    def __init__(self, params, logger=None):
        self.identity = getattr(params, "mmseqs_identity", 0.9)
        self.coverage = getattr(params, "mmseqs_coverage", 0.8)
        self.cov_mode = getattr(params, "mmseqs_cov_mode", 0)
        self.threads = getattr(params, "mmseqs_threads", 1)
        self.mmseqs_bin = getattr(params, "mmseqs_bin", "mmseqs")
        self.mmseq_min = getattr(params, "mmseq_min", 2)

        self.logger = logger

    def cluster(self, sequences):

        if not sequences:
            return {}, {}

        with tempfile.TemporaryDirectory() as tmp:

            input_fasta = os.path.join(tmp, "input.fasta")
            input_db = os.path.join(tmp, "input_db")
            cluster_db = os.path.join(tmp, "cluster_db")
            output_tsv = os.path.join(tmp, "clusters.tsv")
            tmp_dir = os.path.join(tmp, "tmp")

            os.makedirs(tmp_dir)

            self._write_fasta(
                sequences,
                input_fasta
            )

            self._run([
                self.mmseqs_bin,
                "createdb",
                input_fasta,
                input_db,
            ])

            self._run([
                self.mmseqs_bin,
                "cluster",
                input_db,
                cluster_db,
                tmp_dir,
                "--min-seq-id",
                str(self.identity),
                "-c",
                str(self.coverage),
                "--cov-mode",
                str(self.cov_mode),
                "--threads",
                str(self.threads),
            ])

            self._run([
                self.mmseqs_bin,
                "createtsv",
                input_db,
                input_db,
                cluster_db,
                output_tsv,
            ])

            clusters = self._read_clusters(
                output_tsv,
                sequences
            )

        clusters = {
            i: members
            for i, members in clusters.items()
            if len(members) >= self.mmseq_min
        }
        info = {
            "method": "mmseqs2",
            "identity": self.identity,
            "coverage": self.coverage,
            "cov_mode": self.cov_mode,
        }

        return clusters, info

    def _write_fasta(self, sequences, filename):

        with open(filename, "w") as f:

            for i, seq in enumerate(sequences):
                f.write(f">seq_{i}\n")
                f.write(f"{seq}\n")

    def _read_clusters(self, filename, sequences):

        clusters = defaultdict(list)

        with open(filename) as f:

            for line in f:

                line = line.strip()

                if not line:
                    continue

                representative, member = line.split("\t")

                rep_id = int(
                    representative.replace("seq_", "")
                )

                member_id = int(
                    member.replace("seq_", "")
                )

                clusters[rep_id].append(
                    sequences[member_id]
                )

        return {
            i: members
            for i, members in enumerate(clusters.values())
        }

    def _run(self, command):

        if self.logger:
            self.logger.info(
                "Running: %s",
                " ".join(command)
            )

        subprocess.run(
            command,
            check=True,
        )


def attributes(params):

    return {
        "mmseqs_identity": getattr(
            params,
            "mmseqs_identity",
            0.95
        ),
        "mmseqs_coverage": getattr(
            params,
            "mmseqs_coverage",
            0.8
        ),
        "mmseqs_cov_mode": getattr(
            params,
            "mmseqs_cov_mode",
            0
        ),
        "mmseqs_threads": getattr(
            params,
            "mmseqs_threads",
            1
        ),
        "mmseqs_bin": getattr(
            params,
            "mmseqs_bin",
            "mmseqs"
        ),
        "mmseq_min": getattr(
            params,
            "mmseq_min",
            2
        ),
    }