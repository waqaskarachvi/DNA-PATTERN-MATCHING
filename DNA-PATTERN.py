# ==========================
# DNA Pattern Matching Analyzer 🧬 (Simplified)
# ==========================
import streamlit as st
import re, time, io
import pandas as pd
import matplotlib.pyplot as plt
from collections import deque, defaultdict

# ==========================
# 1️⃣ Page Config & Styling
# ==========================
st.set_page_config(page_title="DNA Pattern Matching Analyzer", page_icon="🧬", layout="wide")

# Simple styling for dark theme and buttons
st.markdown("""
<style>
body, .stApp { background-color: #0E1117; color: #FAFAFA; }
h1,h2,h3 { color: #00B4D8 !important; }
.stButton>button { background-color:#00B4D8;color:white;font-weight:bold;border-radius:8px;padding:0.6em 1.2em; }
.stButton>button:hover { background-color:#0077B6; }
</style>
""", unsafe_allow_html=True)

# ==========================
# 2️⃣ Header
# ==========================
st.markdown("<h1 style='text-align:center;'>🧬 DNA Pattern Matching Analyzer</h1>", unsafe_allow_html=True)
st.markdown("<p style='text-align:center;'>Compare DNA pattern matching algorithms</p>", unsafe_allow_html=True)
st.markdown("---", unsafe_allow_html=True)

# ==========================
# 3️⃣ Upload / Enter DNA Sequences
# ==========================
uploaded_files = st.file_uploader("📁 Upload FASTA files (multiple allowed)", type=["fasta","fa","txt"], accept_multiple_files=True)

sequences = {}
if uploaded_files:
    for f in uploaded_files:
        # Read file
        text = f.getvalue().decode("utf-8").strip()
        lines = text.split("\n")
        header = lines[0] if lines[0].startswith(">") else f.name
        seq = "".join([line.strip() for line in lines if not line.startswith(">")]).upper()
        seq = re.sub(r'[^ATCG]', '', seq)  # Remove invalid characters
        sequences[header] = seq
        st.success(f"✅ Loaded: {header} ({len(seq)} bp)")
else:
    # Manual input
    seq_input = st.text_area("🧬 Enter DNA Sequence manually", height=120).strip().upper()
    if seq_input:
        sequences["Manual Entry"] = re.sub(r'[^ATCG]', '', seq_input)

# ==========================
# 4️⃣ Enter Patterns & Select Algorithms
# ==========================
pattern_input = st.text_input("🔍 Enter Pattern(s) (comma separated)")
algorithms = ["Naïve Search", "KMP", "Boyer–Moore", "Rabin–Karp", "Aho–Corasick"]
selected_algos = st.multiselect("⚙️ Select Algorithms", algorithms, default=algorithms)

# ==========================
# 5️⃣ Define Algorithms
# ==========================

# ----- 5a. Naïve Search -----
def naive_search(text, pattern):
    matches = []
    for i in range(len(text)-len(pattern)+1):
        if text[i:i+len(pattern)] == pattern:
            matches.append(i)
    return matches

# ----- 5b. KMP (Knuth–Morris–Pratt) -----
def kmp_search(text, pattern):
    # Step 1: Build LPS (Longest Prefix Suffix) array
    lps = [0]*len(pattern)
    j = 0
    for i in range(1,len(pattern)):
        while j>0 and pattern[i] != pattern[j]:
            j = lps[j-1]
        if pattern[i]==pattern[j]:
            j+=1
            lps[i]=j

    # Step 2: Search using LPS to skip unnecessary comparisons
    res = []
    j = 0
    for i in range(len(text)):
        while j>0 and text[i]!=pattern[j]:
            j = lps[j-1]
        if text[i]==pattern[j]:
            j+=1
        if j==len(pattern):
            res.append(i-j+1)
            j=lps[j-1]
    return res

# ----- 5c. Boyer–Moore -----
def boyer_moore_search(text, pattern):
    m, n = len(pattern), len(text)
    # Build Bad Character Table
    bad_char = {pattern[i]:i for i in range(m)}
    res = []
    s = 0
    while s <= n-m:
        j = m-1
        while j>=0 and pattern[j]==text[s+j]:
            j-=1
        if j<0:
            res.append(s)
            s += (m - bad_char.get(text[s+m],-1)) if s+m<n else 1
        else:
            s += max(1, j-bad_char.get(text[s+j],-1))
    return res

# ----- 5d. Rabin–Karp -----
def rabin_karp(text, pattern, d=256, q=101):
    m, n = len(pattern), len(text)
    h = pow(d,m-1) % q
    p = t = 0
    res = []
    # Compute initial hash
    for i in range(m):
        p = (d*p + ord(pattern[i])) % q
        t = (d*t + ord(text[i])) % q
    # Slide over text
    for s in range(n-m+1):
        if p==t and text[s:s+m]==pattern:
            res.append(s)
        if s<n-m:
            t = (d*(t-ord(text[s])*h) + ord(text[s+m])) % q
            if t<0: t+=q
    return res

# ----- 5e. Aho–Corasick (for multiple patterns) -----
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

    # Build trie for all patterns
    def build_trie(self, patterns):
        for pat in patterns:
            node = self.root
            for c in pat:
                if c not in node.children:
                    node.children[c] = AhoNode()
                node = node.children[c]
            node.output.append(pat)

    # Set failure links to allow efficient backtracking
    def build_failure_links(self):
        queue = deque()
        for child in self.root.children.values():
            child.fail = self.root
            queue.append(child)
        while queue:
            node = queue.popleft()
            for c, child in node.children.items():
                f = node.fail
                while f and c not in f.children:
                    f = f.fail
                child.fail = f.children[c] if f and c in f.children else self.root
                child.output += child.fail.output
                queue.append(child)

    # Search all patterns at once
    def search(self, text):
        node = self.root
        res = defaultdict(list)
        for i,c in enumerate(text):
            while node and c not in node.children:
                node = node.fail
            if not node:
                node = self.root
                continue
            node = node.children[c]
            for pat in node.output:
                res[pat].append(i-len(pat)+1)
        return res

# ==========================
# 6️⃣ Run Analysis
# ==========================
if "results_stored" not in st.session_state:
    st.session_state.results_stored = None

if st.button("🔍 Search Pattern"):
    if not sequences or not pattern_input:
        st.warning("⚠️ Enter sequences and patterns first!")
    else:
        patterns_list = [p.strip() for p in pattern_input.split(",") if p.strip()]
        all_results = []

        for header, seq in sequences.items():
            st.markdown(f"## 🧫 Results for **{header}** ({len(seq)} bp)")
            results = []

            for algo in selected_algos:
                start = time.time()
                if algo == "Aho–Corasick":
                    if len(patterns_list)<2:
                        st.warning("⚠️ Aho–Corasick needs multiple patterns"); continue
                    matches = AhoCorasick(patterns_list).search(seq)
                    elapsed = time.time()-start
                    total_matches = sum(len(pos) for pos in matches.values())
                    results.append({"Algorithm":algo,"Matches":total_matches,"Time (s)":round(elapsed,5)})
                else:
                    # Run each pattern sequentially
                    total_matches = 0
                    func_map = {"Naïve Search": naive_search,
                                "KMP": kmp_search,
                                "Boyer–Moore": boyer_moore_search,
                                "Rabin–Karp": rabin_karp}
                    for pat in patterns_list:
                        total_matches += len(func_map[algo](seq, pat))
                    elapsed = time.time()-start
                    results.append({"Algorithm":algo,"Matches":total_matches,"Time (s)":round(elapsed,5)})

            df = pd.DataFrame(results)
            st.dataframe(df, use_container_width=True)
            all_results.append(df)

            # Simple bar chart for algorithm times
            fig, ax = plt.subplots()
            ax.bar(df["Algorithm"], df["Time (s)"], color="#00B4D8")
            ax.set_ylabel("Time (s)")
            ax.set_title(f"Performance for {header}")
            plt.xticks(rotation=45)
            st.pyplot(fig)

        st.session_state.results_stored = pd.concat(all_results, ignore_index=True)

# ==========================
# 7️⃣ Download CSV
# ==========================
if st.session_state.results_stored is not None:
    csv_buffer = io.StringIO()
    st.session_state.results_stored.to_csv(csv_buffer, index=False)
    st.download_button("📥 Download Results as CSV", data=csv_buffer.getvalue(), file_name="dna_results.csv")
