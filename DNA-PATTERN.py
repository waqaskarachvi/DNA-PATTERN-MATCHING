import streamlit as st
import pandas as pd
from collections import deque, defaultdict

# -----------------------------------------------------------------------------------
# SEARCH ALGORITHMS (SINGLE PATTERN)
# -----------------------------------------------------------------------------------

def naive_search(text, pattern):
    if not pattern or len(pattern) > len(text):
        return []
    occurrences = []
    n = len(text)
    m = len(pattern)
    for i in range(n - m + 1):
        if text[i:i + m] == pattern:
            occurrences.append(i)
    return occurrences


def compute_lps(pattern):
    if not pattern:
        return []
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
                i += 1
    return lps


def kmp_search(text, pattern):
    if not pattern or len(pattern) > len(text):
        return []
    lps = compute_lps(pattern)
    occurrences = []
    i = j = 0
    while i < len(text):
        if pattern[j] == text[i]:
            i += 1
            j += 1
        if j == len(pattern):
            occurrences.append(i - j)
            j = lps[j - 1]
        elif i < len(text) and pattern[j] != text[i]:
            if j != 0:
                j = lps[j - 1]
            else:
                i += 1
    return occurrences


def bad_character_table(pattern):
    table = {}
    length = len(pattern)
    for i in range(length):
        table[pattern[i]] = i
    return table


def boyer_moore_search(text, pattern):
    if not pattern or len(pattern) > len(text):
        return []
    occurrences = []
    bad_char = bad_character_table(pattern)
    n = len(text)
    m = len(pattern)
    s = 0
    while s <= n - m:
        j = m - 1
        while j >= 0 and pattern[j] == text[s + j]:
            j -= 1
        if j < 0:
            occurrences.append(s)
            # shift by full pattern length or based on next char
            if s + m < n:
                occurrences.append if False else None  # placeholder (no-op)
                s += (m - bad_char.get(text[s + m], -1)) if text[s + m] in bad_char else 1
            else:
                s += 1
        else:
            s += max(1, j - bad_char.get(text[s + j], -1))
    return occurrences


def rabin_karp_search(text, pattern, d=256, q=101):
    # Robust Rabin-Karp: return [] if pattern longer than text
    if not pattern or len(pattern) > len(text):
        return []
    n = len(text)
    m = len(pattern)
    occurrences = []
    h = pow(d, m - 1) % q
    p = 0
    t = 0
    # initial hash
    for i in range(m):
        p = (d * p + ord(pattern[i])) % q
        t = (d * t + ord(text[i])) % q
    for s in range(n - m + 1):
        if p == t:
            if text[s:s + m] == pattern:
                occurrences.append(s)
        if s < n - m:
            t = (d * (t - ord(text[s]) * h) + ord(text[s + m])) % q
            if t < 0:
                t += q
    return occurrences


# -----------------------------------------------------------------------------------
# AHO–CORASICK (MULTI-PATTERN)
# -----------------------------------------------------------------------------------

class AhoNode:
    def __init__(self):
        self.children = {}
        self.fail = None
        self.output = []


class AhoCorasick:
    def __init__(self, patterns):
        self.root = AhoNode()
        self.build_trie(patterns)
        self.build_failure_links()

    def build_trie(self, patterns):
        for index, pat in enumerate(patterns):
            node = self.root
            for char in pat:
                if char not in node.children:
                    node.children[char] = AhoNode()
                node = node.children[char]
            node.output.append((index, pat))

    def build_failure_links(self):
        queue = deque()
        # Set fail of direct children to root
        for child in self.root.children.values():
            child.fail = self.root
            queue.append(child)

        while queue:
            current = queue.popleft()
            for char, next_node in current.children.items():
                fail_node = current.fail
                while fail_node and char not in fail_node.children:
                    fail_node = fail_node.fail
                next_node.fail = fail_node.children[char] if (fail_node and char in fail_node.children) else self.root
                # Inherit output
                if next_node.fail and next_node.fail.output:
                    next_node.output += next_node.fail.output
                queue.append(next_node)

    def search(self, text):
        node = self.root
        results = defaultdict(list)

        for i, char in enumerate(text):
            while node is not None and char not in node.children:
                node = node.fail
            if node is None:
                node = self.root
                continue
            node = node.children[char]
            if node.output:
                for pattern_index, pattern in node.output:
                    results[pattern].append(i - len(pattern) + 1)

        return results


# -----------------------------------------------------------------------------------
# STREAMLIT UI
# -----------------------------------------------------------------------------------

st.set_page_config(page_title="DNA Pattern Search Tool", layout="wide")
st.title("🔬 DNA Pattern Search Tool")
st.write("You can either upload a text/FASTA file or paste the sequence manually. Multi-pattern uses Aho–Corasick.")

# Input area: file uploader OR manual paste
uploaded_file = st.file_uploader("Upload DNA file (txt, fasta). If multiple records, upload a single file.", type=["txt", "fasta", "fa"])
manual_seq = st.text_area("Or paste DNA sequence here (A/T/C/G only)", height=150)

# choose source priority: manual > file
sequence = ""
if manual_seq and manual_seq.strip():
    # sanitize: keep only A/T/C/G (uppercase)
    sequence = "".join([c for c in manual_seq.upper() if c in "ATCG"])
elif uploaded_file:
    raw = uploaded_file.read().decode("utf-8")
    # If FASTA, remove headers and newlines
    lines = raw.splitlines()
    seq_lines = [ln.strip() for ln in lines if not ln.startswith(">")]
    sequence = "".join(seq_lines).upper()
    sequence = "".join([c for c in sequence if c in "ATCG"])

if not sequence:
    st.info("No sequence provided yet. Paste a sequence or upload a file to start searching.")
# Mode selection
mode = st.radio("Select Search Mode", ["Single Pattern Search", "Multi-Pattern Search"], index=0)

# Single pattern UI
if mode == "Single Pattern Search":
    pattern = st.text_input("Enter pattern (single)", help="Sequence of A/T/C/G")
    algorithm = st.selectbox("Select algorithm", ["Naive", "KMP", "Boyer-Moore", "Rabin-Karp"])
    if st.button("Search Single Pattern"):
        if not sequence:
            st.warning("Please provide a DNA sequence (paste or upload a file).")
        elif not pattern:
            st.warning("Please enter a pattern to search.")
        else:
            # sanitize pattern
            pat = "".join([c for c in pattern.upper() if c in "ATCG"])
            if pat == "":
                st.error("Pattern contains no valid A/T/C/G characters.")
            else:
                if algorithm == "Naive":
                    found = naive_search(sequence, pat)
                elif algorithm == "KMP":
                    found = kmp_search(sequence, pat)
                elif algorithm == "Boyer-Moore":
                    found = boyer_moore_search(sequence, pat)
                else:
                    found = rabin_karp_search(sequence, pat)

                st.success(f"Found {len(found)} matches.")
                if found:
                    st.write(found)
                    # show first 400 chars with highlighted matches (simple)
                    display_snippet = list(sequence[:400])
                    for pos in found:
                        if pos < 400:
                            for j in range(len(pat)):
