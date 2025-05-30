def default_error_function(kmer_amount):
    """
    Default function to calculate the error probabilities based on the amount of kmers.
    :param kmer_amount:
    :return:
    """
    return kmer_amount ** 2 * 0.000002


def kmer_counting(seq, k=20, upper_bound=1, error_function=default_error_function):
    """
    Generates a dictionary with all k-mers for the given k occuring in the sequence and checks if they appear more than
    once. The error probabilities are calculated based on the count of the k-mers.
    :param seq:
    :param k:
    :param upper_bound:
    :param error_function:
    :return:
    """
    kmer_dict = {}
    kmer_pos = {}
    kmers = []
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        if kmer not in kmer_pos:
            kmer_pos[kmer] = list()
        kmer_pos[kmer].append(i)

        kmers.append(kmer)
        if kmer not in kmer_dict.keys():
            kmer_dict[kmer] = 1
        else:
            kmer_dict[kmer] += 1
    high_occ_kmers = [{'kmer': i, 'startpos': elem, 'endpos': elem + k - 1, 'errorprob': error_function(v),
                       'identifier': i} for i, v in kmer_dict.items() if v > upper_bound for elem in
                      kmer_pos[i]]

    high_occ_kmers = [x for x in high_occ_kmers if x['errorprob'] > 0.0]
    # optimally this check should be done during first list comprehension,
    # but the main speedup will be seen at the frontend

    if high_occ_kmers:
        return high_occ_kmers  # sum(list(zip(*high_occ_kmers.values()))[1]), high_occ_kmers
    else:
        return {}

