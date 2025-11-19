# ==========================
# DNA Pattern Matching Analyzer 🧬 (Streamlit App)
# ==========================
import streamlit as st
import re, time, io
import pandas as pd
import matplotlib.pyplot as plt

# ==========================
# PAGE CONFIG
# ==========================
st.set_page_config(
    page_title="DNA Pattern Matching Analyzer",
    page_icon="🧬",
    layout="wide"
)

# ==========================
# STYLING
# ==========================
st.markdown("""
    <style>
        body, .stApp { background-color: #0E1117; color: #FAFAFA; }
        h1, h2, h3 { color: #00B4D8 !important; }
        .stButton>button {
            background-color: #00B4D8; color: white; font-weight: bold;
            border-radius: 8px; border: none; padding: 0.6em 1.2em;
        }
        .stButton>button:hover { background-color: #0077B6; }
        .result-box {
            background-color: #1E2636; padding: 15px;
            border-radius: 10px; border: 1px solid #00B4D8;
            font-family: monospace; line-height: 1.6;
            word-wrap: break-word;
        }
    </style>
""", unsafe_allow_html=True)

# ==========================
# HEADER
# ==========================
st.markdown("<h1 style='text-align:center;'>🧬 DNA Pattern Matching Analyzer</h1>", unsafe_allow_html=True)
st.markdown("<p style='text-align:center;'>Compare multiple string matching algorithms for DNA sequences</p>", unsafe_allow_html=True)
st.markdown("---", unsafe_allow_html=True)

# ==========================
# FILE UPLOAD / INPUT
# ==========================
uploaded_files = st.file_uploader(
    "📁 Upload FASTA files (you can select multiple)",
    type=["fasta", "fa", "txt"],
    accept_multiple_files=True
)

sequences = {}

if uploaded_files:
    for uploaded_file in uploaded_files:
        fasta = uploaded_file.getvalue().decode("utf-8")
        lines = fasta.strip().split("\n")
        header = lines[0] if lines[0].startswith(">") else uploaded_file.name
        sequence = "".join([l.strip() for l in lines if not l.startswith(">")]).upper()
        sequence = re.sub(r'[^ATCG]', '', sequence)
        sequences[header] = sequence
        st.success(f"✅ Loaded: {header} ({len(sequence)} bp)")
else:
    seq_input = st.text_area(
        "🧬 Enter DNA Sequence",
        placeholder="ATCGGATCGATCG...",
        height=120
    ).strip().upper()
    if seq_input:
        sequences["Manual Entry"] = re.sub(r'[^ATCG]', '', seq_input)

# MULTI-PATTERN INPUT
pattern_input = st.text_input(
    "🔍 Enter Pattern(s) (comma-separated for multiple patterns)",
    placeholder="ATC, GGA, TTA"
).strip().upper()

patterns = [p.strip() for p in pattern_input.split(",") if p.strip()]

# ==========================
# ALGORITHMS
# ==========================

# --- Naive Search ---
def naive_search(text, pattern):
    return [i for i in range(len(text)-len(pattern)+1) if text[i:i+len(pattern)] == pattern]

# --- KMP ---
def kmp_search(text, pattern):
    lps = [0]*len(pattern)
    j = 0
    for i in range(1, len(pattern)):
        while j > 0 and pattern[i] != pattern[j]:
            j = lps[j-1]
        if pattern[i] == pattern[j]:
            j += 1
            lps[i] = j

    res, j = [], 0
    for i in range(len(text)):
        while j > 0 and text[i] != pattern[j]:
            j = lps[j-1]
        if text[i] == pattern[j]:
            j += 1
        if j == len(pattern):
            res.append(i-j+1)
            j = lps[j-1]
    return res

# --- Boyer Moore ---
def boyer_moore_search(text, pattern):
    m, n = len(pattern), len(text)
    bad_char = {pattern[i]: i for i in range(m)}
    res, s = [], 0

    while s <= n - m:
        j = m - 1
        while j >= 0 and pattern[j] == text[s + j]:
            j -= 1
        if j < 0:
            res.append(s)
            s += (m - bad_char.get(text[s+m], -1)) if s + m < n else 1
        else:
            s += max(1, j - bad_char.get(text[s+j], -1))
    return res

# --- Rabin Karp ---
def rabin_karp(text, pattern, d=256, q=101):
    m, n = len(pattern), len(text)
    p = t = 0
    h = pow(d, m-1) % q
    res = []

    for i in range(m):
        p = (d*p + ord(pattern[i])) % q
        t = (d*t + ord(text[i])) % q

    for s in range(n-m+1):
        if p == t and text[s:s+m] == pattern:
            res.append(s)
        if s < n-m:
            t = (d*(t - ord(text[s])*h) + ord(text[s+m])) % q
            if t < 0:
                t += q

    return res

# --- AHO CORASICK (full implementation) ---
class AhoCorasick:
    def __init__(self, patterns):
        self.num_nodes = 1
        self.edges = [{}]
        self.fail = [0]
        self.output = [[]]

        for pattern in patterns:
            self._add_pattern(pattern)

        self._build_fail_links()

    def _add_pattern(self, pattern):
        node = 0
        for char in pattern:
            if char not in self.edges[node]:
                self.edges[node][char] = self.num_nodes
                self.edges.append({})
                self.fail.append(0)
                self.output.append([])
                self.num_nodes += 1
            node = self.edges[node][char]
        self.output[node].append(pattern)

    def _build_fail_links(self):
        from collections import deque
        q = deque()

        for char, nxt in self.edges[0].items():
            q.append(nxt)

        while q:
            r = q.popleft()
            for char, nxt in self.edges[r].items():
                q.append(nxt)
                f = self.fail[r]
                while f > 0 and char not in self.edges[f]:
                    f = self.fail[f]
                self.fail[nxt] = self.edges[f].get(char, 0)
                self.output[nxt] += self.output[self.fail[nxt]]

    def search(self, text):
        node = 0
        results = []  # (pattern, index)

        for i, char in enumerate(text):
            while node > 0 and char not in self.edges[node]:
                node = self.fail[node]

            node = self.edges[node].get(char, 0)

            for p in self.output[node]:
                results.append((p, i - len(p) + 1))

        return results

def aho_corasick_search(text, patterns):
    ac = AhoCorasick(patterns)
    matches = ac.search(text)
    results = {p: [] for p in patterns}
    for pat, pos in matches:
        results[pat].append(pos)
    return results


# Algorithm dictionary
algo_funcs = {
    "Naïve Search": lambda t, pats: {p: naive_search(t, p) for p in pats},
    "KMP": lambda t, pats: {p: kmp_search(t, p) for p in pats},
    "Boyer–Moore": lambda t, pats: {p: boyer_moore_search(t, p) for p in pats},
    "Rabin–Karp": lambda t, pats: {p: rabin_karp(t, p) for p in pats},
    "Aho–Corasick": lambda t, pats: aho_corasick_search(t, pats)
}

algorithms = list(algo_funcs.keys())

selected_algos = st.multiselect(
    "⚙️ Select Algorithms",
    algorithms,
    default=["KMP", "Boyer–Moore"]
)

# ==========================
# RUN ANALYSIS
# ==========================
if "results_stored" not in st.session_state:
    st.session_state.results_stored = None

if st.button("🔍 Search Pattern"):
    if not sequences or not patterns:
        st.warning("⚠️ Please enter both sequence(s) and at least one pattern.")
    else:
        all_results = []

        for header, dna_sequence in sequences.items():
            st.markdown(f"## 🧫 Results for **{header}** ({len(dna_sequence)} bp)")

            results = []
            for algo in selected_algos:
                start = time.time()

                match_dict = algo_funcs[algo](dna_sequence, patterns)
                match_count = sum(len(v) for v in match_dict.values())

                elapsed = round(time.time() - start, 5)

                results.append({
                    "Sequence": header,
                    "Algorithm": algo,
                    "Patterns": ", ".join(patterns),
                    "Total Matches": match_count,
                    "Time (s)": elapsed
                })

            df = pd.DataFrame(results)
            all_results.append(df)

            st.markdown("### 📊 Algorithm Comparison")
            st.dataframe(df, use_container_width=True)

            # Highlight Results
            st.markdown("### 🎨 Sequence Visualization")

            colors = ["#FFD60A", "#FF6B6B", "#4CC9F0", "#7CFC00", "#FF1493"]

            for algo in selected_algos:
                match_dict = algo_funcs[algo](dna_sequence, patterns)

                highlighted = list(dna_sequence)

                for idx, pat in enumerate(patterns):
                    col = colors[idx % len(colors)]
                    for pos in match_dict[pat]:
                        for j in range(len(pat)):
                            if pos + j < len(highlighted):
                                highlighted[pos + j] = (
                                    f"<span style='color:{col}; font-weight:bold'>{highlighted[pos + j]}</span>"
                                )

                st.markdown(f"**{algo}:**", unsafe_allow_html=True)
                st.markdown(
                    f"<div class='result-box'>{''.join(highlighted[:400])}...</div>",
                    unsafe_allow_html=True
                )

            # Performance Chart
            st.markdown("### 📈 Performance Chart")
            fig, ax = plt.subplots()
            ax.bar(df["Algorithm"], df["Time (s)"])
            ax.set_ylabel("Execution Time (s)")
            ax.set_title("Algorithm Performance Comparison")
            st.pyplot(fig)

        combined_df = pd.concat(all_results, ignore_index=True)
        st.session_state.results_stored = combined_df
        st.success("✅ Analysis complete! CSV ready for download.")

# ==========================
# DOWNLOAD CSV SECTION
# ==========================
if st.session_state.results_stored is not None:
    csv_buffer = io.StringIO()
    st.session_state.results_stored.to_csv(csv_buffer, index=False)
    st.download_button(
        label="📥 Download Results as CSV",
        data=csv_buffer.getvalue(),
        file_name="dna_pattern_results.csv",
        mime="text/csv"
    )
