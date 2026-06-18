# Author: Weston Hanson
# Place: Fred Hutch Cancer Center, Seattle, WA
# Date Created: 06/08/26
# Purpose: List of entropy functions used un main.py and entorpy_by_chromosome.py

from python_scripts.imports import *

# List of entropy functions
def base_entropy(counts_probs):
    if isinstance(counts_probs, dict):
        entropy = sum(val * math.log2((1/val)) for i, val in counts_probs.items())
    if isinstance(counts_probs, np.ndarray):
        counts_probs = counts_probs[counts_probs > 0] # Makes matrix into a list of non-zero values
        entropy = np.sum(counts_probs * np.log2((1/counts_probs)))
    return entropy

def base_entropy_plus_CN(counts_probs):
    CN = sum(counts_probs.keys())
    return sum(val * math.log2((1/val)) for i, val in counts_probs.items()) + CN

def entropy_w1_error_term(counts_probs):
    return sum(val * abs(i - 2) * math.log2((1 / val)) for i, val in counts_probs.items())

def entropy_w2_error_term(counts_probs):
    return sum(val * math.log2(abs(i - 2) + 1) * math.log2((1 / val)) for i, val in counts_probs.items())

def entropy_w3_error_term(counts_probs):
    return sum(val * i * math.log2((1 / val)) for i, val in counts_probs.items())

def entropy_w4_error_term(counts_probs):
    return sum(val * max(0, abs(i-2) - 1) * math.log2((1 / val)) for i, val in counts_probs.items()) 

def entropy_w5_error_term(counts_probs):
    return sum(val * max(0, abs(i-2) - 2) * math.log2((1 / val)) for i, val in counts_probs.items()) 

def entropy_w6_error_term(counts_probs):
    return sum(val * max(0, abs(i-2) - 3) * math.log2((1 / val)) for i, val in counts_probs.items())

def entropy_w7_error_term(counts_probs):
    return sum(val * max(0, abs(i-2) - 0.5) * math.log2((1 / val)) for i, val in counts_probs.items())

def entropy_w8_error_term(counts_probs):
    return sum(val * pow(max(0, abs(i-2) - 0.5), 2) * math.log2((1 / val)) for i, val in counts_probs.items())

def entropy_w9_error_term(counts_probs):
    return sum(val * math.log2(abs(i - 2) + 0.5) * math.log2((1 / val)) for i, val in counts_probs.items())

def entropy_w10_error_term(counts_probs):
    return sum(val * math.log2(abs(i + 8)) * math.log2((1 / val)) for i, val in counts_probs.items())

def entropy_w11_error_term(counts_probs):
    return sum(val * math.log2(0.1 * abs(i - 2) + 2) * math.log2((1 / val)) for i, val in counts_probs.items())

def entropy_w12_error_term(counts_probs):
    return sum(val * math.log2(0.3 * abs(i - 2) + 2) * math.log2((1 / val)) for i, val in counts_probs.items())

def entropy_w13_error_term(counts_probs):
    return sum(val * math.log2(0.5 * abs(i - 2) + 2) * math.log2((1 / val)) for i, val in counts_probs.items())

def entropy_w14_error_term(counts_probs):
    return sum(val * math.log2(0.7 * abs(i - 2) + 2) * math.log2((1 / val)) for i, val in counts_probs.items())

def entropy_w15_error_term(counts_probs):
    return sum(val * math.log2(abs(i-2) + 2) * math.log2((1 / val)) for i, val in counts_probs.items())

def entropy_w16_error_term(counts_probs):
    return sum(val * math.log2(2 * abs(i - 2) + 2) * math.log2((1 / val)) for i, val in counts_probs.items())

def entropy_w17_error_term(counts_probs):
    return sum(val * math.log2(3 * abs(i - 2) + 2) * math.log2((1 / val)) for i, val in counts_probs.items())

def entropy_w18_error_term(counts_probs):
    return sum(val * math.log2(5 * abs(i - 2) + 2) * math.log2((1 / val)) for i, val in counts_probs.items())

def entropy_w19_error_term(counts_probs):
    return sum(val * math.log2(7 * abs(i - 2) + 2) * math.log2((1 / val)) for i, val in counts_probs.items())

def entropy_w20_error_term(counts_probs):
    return sum(val * math.log2(10 * abs(i - 2) + 2) * math.log2((1 / val)) for i, val in counts_probs.items())

def base_entropy_exponentiate(counts_probs):
    if isinstance(counts_probs, dict):
        entropy = sum(val * math.log2((1/val)) for i, val in counts_probs.items())
    if isinstance(counts_probs, np.ndarray):
        counts_probs = counts_probs[counts_probs > 0] # Makes matrix into a list of non-zero values
        entropy = np.sum(counts_probs * np.log2((1/counts_probs)))
    return pow(2, entropy)

def base_entropy_hn_normalized(counts_probs):
    if isinstance(counts_probs, dict):
        entropy = sum(val * math.log2((1/val)) for i, val in counts_probs.items())
    if isinstance(counts_probs, np.ndarray):
        counts_probs = counts_probs[counts_probs > 0] # Makes matrix into a list of non-zero values
        entropy = np.sum(counts_probs * np.log2((1/counts_probs)))
    return pow(2, entropy) / math.log2(len(counts_probs) + 1)

def base_entropy_hn_cohort_normalized(counts_probs):
    if isinstance(counts_probs, dict):
        entropy = sum(val * math.log2((1/val)) for i, val in counts_probs.items())
    if isinstance(counts_probs, np.ndarray):
        counts_probs = counts_probs[counts_probs > 0] # Makes matrix into a list of non-zero values
        entropy = np.sum(counts_probs * np.log2((1/counts_probs)))
    return pow(2, entropy)
    
def base_entropy_hn_cohort_normalized_from_probability_matrix(counts_probs):
    entropy = sum(val * math.log2((1/val)) for i, val in counts_probs.items())
    return pow(2, entropy)

def base_entropy_hn_normalized_plus_CN(counts_probs):
    CN = sum(counts_probs.keys())
    entropy = sum(val * math.log2((1/val)) for i, val in counts_probs.items())
    return (pow(2, entropy) / math.log2(len(counts_probs) + 1)) + CN

def base_entropy_hn_normalized_logged(counts_probs):
    entropy = sum(val * math.log2((1/val)) for i, val in counts_probs.items())
    entropy_score = pow(2, entropy) / math.log2(len(counts_probs) + 1)
    return math.log2(entropy_score)

def entropy_w1_error_term_hn_normalized(counts_probs):
    entropy = sum(val * abs(i - 2) * math.log2((1 / val)) for i, val in counts_probs.items())
    return pow(2, entropy) / math.log2(len(counts_probs) + 1)

def entropy_w2_error_term_log_normalized(counts_probs):
    entropy = sum(val * math.log2(abs(i - 2) + 1) * math.log2((1 / val)) for i, val in counts_probs.items())
    return entropy / math.log2(len(counts_probs) + 1)

def entropy_w2_error_term_hn_normalized(counts_probs):
    entropy = sum(val * math.log2(abs(i - 2) + 1) * math.log2((1 / val)) for i, val in counts_probs.items())
    return pow(2, entropy) / math.log2(len(counts_probs) + 1)

def entropy_w2_error_term_hn_normalized_logged(counts_probs):
    entropy = sum(val * math.log2(abs(i - 2) + 1) * math.log2((1 / val)) for i, val in counts_probs.items())
    entropy_score = pow(2, entropy) / math.log2(len(counts_probs) + 1)
    return math.log2(entropy_score)

def entropy_w3_error_term_hn_normalized(counts_probs):
    entropy = sum(val * i * math.log2((1 / val)) for i, val in counts_probs.items())
    return pow(2, entropy) / math.log2(len(counts_probs) + 1)
