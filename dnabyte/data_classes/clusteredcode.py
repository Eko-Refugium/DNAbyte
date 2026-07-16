from .base import Data


class ClusteredDNA(Data):
    """
    Represents clustered DNA sequences after post-sequencing clustering.

    Data format:

    {
        cluster_id: [
            sequence1,
            sequence2,
            ...
        ]
    }

    All clustering methods must output this format.
    """


    def __init__(self, data):

        self._validate_cluster_data(data)

        self.data = data
        self.file_paths = []

        self.num_clusters = len(data)

        self.num_sequences = sum(
            len(sequences)
            for sequences in data.values()
        )

        self.total_length = sum(
            len(seq)
            for sequences in data.values()
            for seq in sequences
        )

        self.average_length = (
            self.total_length / self.num_sequences
            if self.num_sequences > 0
            else 0
        )


    def _validate_cluster_data(self, data):

        if not isinstance(data, dict):
            raise TypeError(
                f"ClusteredDNA data must be dict, got {type(data).__name__}"
            )


        for cluster_id, sequences in data.items():

            if not isinstance(sequences, list):
                raise TypeError(
                    f"Cluster {cluster_id} must contain a list"
                )


            for seq in sequences:

                if not isinstance(seq, str):
                    raise TypeError(
                        "Cluster sequences must be strings"
                    )

                invalid = set(seq) - {"A", "C", "G", "T"}

                if invalid:
                    raise ValueError(
                        f"Invalid DNA bases found: {invalid}"
                    )


    def get_cluster(self, cluster_id):

        return self.data[cluster_id]


    def get_sequences(self):

        output = []

        for cluster in self.data.values():
            output.extend(cluster)

        return output


    def add_cluster(self, cluster_id, sequences):

        self.data[cluster_id] = sequences

        self._update_stats()


    def _update_stats(self):

        self.num_clusters = len(self.data)

        self.num_sequences = sum(
            len(v)
            for v in self.data.values()
        )

        self.total_length = sum(
            len(seq)
            for cluster in self.data.values()
            for seq in cluster
        )

        self.average_length = (
            self.total_length / self.num_sequences
            if self.num_sequences
            else 0
        )


    def __len__(self):

        return self.num_clusters


    def __iter__(self):

        return iter(self.data.items())


    def __getitem__(self, key):

        return self.data[key]


    def __str__(self):

        output = (
            f"Type: ClusteredDNA\n"
            f"Clusters: {self.num_clusters}\n"
            f"Sequences: {self.num_sequences}\n"
            f"Total length: {self.total_length} bp\n"
            f"Average length: {self.average_length:.2f} bp\n"
        )


        output += "Cluster sizes:\n"

        for cluster_id, sequences in list(self.data.items())[:10]:

            output += (
                f"  {cluster_id}: "
                f"{len(sequences)} sequences\n"
            )


        if self.num_clusters > 10:
            output += (
                f"  ... {self.num_clusters-10} more clusters\n"
            )


        return output