import streamlit as st
import pandas as pd
from collections import deque, defaultdict
import time
import matplotlib.pyplot as plt

# ------------------ SEARCH ALGORITHMS ------------------

def naive_search(text, pattern):
    if not pattern or len(pattern) > len(text):
        return []
    occurrences = []
    n, m = len(text), len(pattern)
    for i in range(n - m + 1):
        if text[i:i + m] == pattern:
            occurrences.append(i)
    return occurrences

def compute_lps(pattern):
    lps = [0] * len(pattern)
    j = 0
    i = 1
    while i < len(pattern):
        if pattern[i] == pattern[j]:
            j += 1
            lps[i] = j
            i += 1
        else:
            if j != 0:
                j = lps[j - 1]
            else:
                lps[i] = 0
