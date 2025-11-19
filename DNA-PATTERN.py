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
            s += m - bad_char.get(text[s + m], -1) if s + m < n else 1
        else:
            s += max(1, j - bad_char.get(text[s + j], -1))
    return occurrences


def rabin_karp_search(text, pattern, d=256, q=101):
    if not pattern or len(pattern) > len(text):
        return []
    n = len(text)
    m = len(pattern)
    occurrences = []
    h = pow(d, m - 1) % q
    p = 0
    t = 0
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
            for pattern_index, pattern in node.output:
                results[pattern].append(i - len(pattern) + 1)
        return results

# -----------------------------------------------------------------------------------
# STREAMLIT UI
# -----------------------------------------------------------------------------------

st.set_page_config(page_title="DNA Pattern Search Tool", layout="wide")
st.title("🔬 DNA Pattern Search Tool")
st.write("Single pattern: Naive/KMP/Boyer–Moore/Rabin–Karp. Multi-pattern: Aho–Corasick.")

uploaded_file = st.file_uploader("Upload DNA file (txt/fasta)", type=["txt", "fa", "fasta"])
manual_seq = st.text_area("Or paste DNA sequence (A/T/C/G only)", height=150)

# Prepare sequence
sequence = ""
if manual_seq and manual_seq.strip():
    sequence = "".join([c for c in manual_seq.upper() if c in "ATCG"])
elif uploaded_file:
    raw = uploaded_file.read().decode("utf-8")
    lines = raw.splitlines()
    seq_lines = [ln.strip() for ln in lines if not ln.startswith(">")]
    sequence = "".join(seq_lines).upper()
    sequence = "".join([c for c in sequence if c in "ATCG"])

if not sequence:
    st.info("No DNA sequence provided. Paste it or upload a file.")

mode = st.radio("Select Search Mode", ["Single Pattern Search", "Multi-Pattern Search"], index=0)

# ------------------ Single Pattern ------------------
if mode == "Single Pattern Search":
    pattern = st.text_input("Enter single pattern (A/T/C/G)")
    algorithm = st.selectbox("Select algorithm", ["Naive", "KMP", "Boyer-Moore", "Rabin-Karp"])
    if st.button("Search Single Pattern"):
        if not sequence:
            st.warning("Provide a DNA sequence first!")
        elif not pattern:
            st.warning("Enter a pattern to search.")
        else:
            pat = "".join([c for c in pattern.upper() if c in "ATCG"])
            if not pat:
                st.error("Pattern contains no valid characters (A/T/C/G).")
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
                    snippet = list(sequence[:400])
                    for pos in found:
                        if pos < 400:
                            for j in range(len(pat)):
                                if pos + j < 400:
                                    snippet[pos + j] = f"**{snippet[pos + j]}**"
                    st.markdown("".join(snippet))
                else:
                    st.info("No matches found.")

# ------------------ Multi Pattern ------------------
else:
    patterns_raw = st.text_area("Enter multiple patterns (one per line, min 2)")
    if st.button("Search Multiple Patterns (Aho–Corasick)"):
        if not sequence:
            st.warning("Provide a DNA sequence first!")
        else:
            patterns = [p.strip().upper() for p in patterns_raw.splitlines() if p.strip()]
            patterns = ["".join([c for c in p if c in "ATCG"]) for p in patterns]
            patterns = [p for p in patterns if p]
            if len(patterns) < 2:
                st.error("Enter at least 2 valid patterns (A/T/C/G) for multi-pattern search.")
            else:
                aho = AhoCorasick(patterns)
                results = aho.search(sequence)
                if not results:
                    st.info("No matches found.")
                else:
                    rows = []
                    for pat, pos_list in results.items():
                        for pos in pos_list:
                            rows.append((pat, pos))
                    df = pd.DataFrame(rows, columns=["Pattern", "Position"])
                    st.write(df)

                    csv = df.to_csv(index=False).encode("utf-8")
                    st.download_button("Download CSV", csv, "aho_results.csv", "text/csv")
