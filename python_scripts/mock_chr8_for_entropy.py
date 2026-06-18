# Author: Weston Hanson
# Place: Fred Hutch Cancer Center, Seattle, WA
# Date Created: 03/09/2025
# Purpose: Calculates a penalized Shannon entropy score for different mock chr8s to test what the entropy is for different CN patterns.

from python_scripts.imports import *

def mock_chr8_entropy(values, function):
    '''
    Note: Adapted from entropy_by_chromosome.py's find_entropy_per_chromosome function. Did not change the code. 
    Parameters:
    -----------
        values: List
            List of integers (mock CN states).

        function: function
            The entropy function.

    Function:
    ---------
        - Calculates a Shannon entropy score for a list of integers, corrected for ploidy 2.
    
    Returns:
    -------
        entropy_score: float
            Shannon entorpy score from values.
    '''
    clean_values = [v for v in values if v != "NA"]
    len_of_chrom = len(clean_values)
    counts = Counter(clean_values)
    counts_probs = {float(key): count / len_of_chrom for key, count in counts.items()}
    # counts_probs = {1: 0.44, 2: 0.06, 3: 0.33, 4: 0.17, 7: 0.004}

    # New map variable to calculate entropy
    # entropy_score = sum(val * abs(i - 2) * (math.log2((1 / val) + 1)) for i, val in counts_probs.items())
    entropy_score = function(counts_probs)

    return entropy_score

##################################
# COMPARING SINGLE ENTROPY SCORES 
##################################

def base_entropy_hn_normalized(counts_probs):
    entropy = sum(val * math.log2((1/val)) for i, val in counts_probs.items())
    print(f'counts_probs: {counts_probs}')
    for i, val in counts_probs.items():
        print(f'{i}: {val * math.log2((1/val))}')
    print(f'summation: {entropy}')
    print(f'K+1: {len(counts_probs) + 1}')
    print(f'2^sum: {pow(2, entropy)}')
    print(f'denominator: {math.log2(len(counts_probs) + 1)}')
    print(f'entropy: {pow(2, entropy) / math.log2(len(counts_probs) + 1)}')
    return pow(2, entropy) / math.log2(len(counts_probs) + 1)

def base_entropy_hn_cohort_normalized(counts_probs):
    entropy = sum(val * math.log2((1/val)) for i, val in counts_probs.items())
    print(f'counts_probs: {counts_probs}')
    for i, val in counts_probs.items():
        print(f'{i}: {val * math.log2((1/val))}')
    print(f'summation: {entropy}')
    print(f'2^sum: {pow(2, entropy)}')
    return pow(2, entropy)

mock_chr_dict = {
    # "CN neutral": [2] * 50,
    # "CN neutral + one CN gain": [2] * 49 + [3],
    # "CN 8": [8] * 50,
    # "2/3 oscillating": [2, 3] * 25,
    # "> low CN": [1] * 10 + [2] * 11 + [3] * 12 + [4] * 13 + [5] * 14 + [6] * 15 + [7] * 16 + [8] * 17 + [9] * 18 + [10] * 19,
    # "< high CN": sum([[cn] * ((cn % 9) + 2) for cn in range(1, 17)], []),
    # 'low_deep_seg_entropy_mod': [key for key, count in {2: 15, 3:85}.items() for _ in range(count)],
    # "FH0084_ULP": [3] * 73 + [1] * 44 + [4] * 10 + [1] * 31 + [2] * 12 + [4] * 28 + [1] * 19,
    "FH0084_deep_ichor": [3] * 73 + [2] * 44 + [4] * 10 + [2] * 43 + [4] * 28 + [2] * 19,
    # "FH0861_ULP": [2] * 5 + [4] * 28 + [2] * 32 + [1] * 31+ [3] * 16 + [1] * 23+ [3] * 52,
    # "FH0842_ULP": [4] * 10 + [1] * 17 + [3] * 6 + [4] * 10 + [2] * 15 + [1] * 9 + [2] * 47 + [1] * 18 + [2] * 42 + [1] * 6 + [2] * 8,
    # "FH0842_ULP_subclone": [4] * 10 + [1.5] * 17 + [2.5] * 6 + [4] * 10 + [2] * 15 + [1] * 9 + [2] * 47 + [1.5] * 18 + [2] * 42 + [1] * 6 + [2] * 8
}
{'4': 20, 1.5: 35, 2.5: 6, '2': 112, '1': 15}
for name, values in mock_chr_dict.items():
    entropy = mock_chr8_entropy(values, base_entropy_hn_normalized)

    # divide by max k
    entropy = entropy / 7

    print(f"{name}: {entropy}\n")