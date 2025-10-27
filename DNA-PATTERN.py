import streamlit as st
import re, time
import pandas as pd
import matplotlib.pyplot as plt

# ============================
# PAGE CONFIG
# ============================
st.set_page_config(page_title="DNA Pattern Matching Analyzer", page_icon="🧬", layout="wide")

# ============================
# STYLE
# ============================
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
        }
        .highlight { color: #FFD60A; font-weight: bold; }
    </style>
""", unsafe_allow_html=True)

# ============================
# HEADER
# ============================
st.markdown("<h1 style='text-align:center;'>🧬 DNA Pattern Matching Analyzer</h1>", unsafe_allow_html=True)
st.markdown("<p style='text-align:center;'>Compare multiple string matching algorithms for DNA sequences</p>", unsafe_allow_html=True)
st.markdown("---", unsafe_allow_html=True)

# ============================
# INPUT SECTION
# ============================
uploaded_file = st.file_uploader("📁 Upload FASTA file (optional)", type=["fasta", "fa", "txt"])
dna_sequence = ""
if uploaded_file:
    fasta = uploaded_file.getvalue().decode("utf-8")
    lines = fasta.strip().split("\n")
    header = lines[0] if lines[0].startswith(">") else "Uploaded Sequence"
    dna_sequence = "".join([l.strip() for l in lines if not l.startswith(">")]).upper()
    dna_sequence = re.sub(r'[^ATCG]', '', dna_sequence)
    st.success(f"✅ Loaded sequence: {header} ({len(dna_sequence)} bp)")
else:
    dna_sequence = st.text_area("🧬 Enter DNA Sequence", placeholder="ATCGGATCGATCG...", height=120).strip().upper()

pattern = st.text_input("🔍 Enter Pattern to Search", placeholder="CGATCGA").strip().upper()

algorithms = ["Naïve Search", "KMP", "Boyer–Moore", "Rabin–Karp", "Aho–Corasick"]
selected_algos = st.multiselect("⚙️ Select Algorithms", algorithms, default=["KMP", "Boyer–Moore"])

# ============================
# ALGORITHMS
# ============================
def naive_search(text, pattern):
    return [i for i in range(len(text)-len(pattern)+1) if text[i:i+len(pattern)] == pattern]

def kmp_search(text, pattern):
    lps = [0]*len(pattern)
    j = 0
    for i in range(1, len(pattern)):
        while j > 0 and pattern[i] != pattern[j]:
            j = lps[j-1]
        if pattern[i] == pattern[j]:
            j += 1; lps[i] = j
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

def boyer_moore_search(text, pattern):
    m, n = len(pattern), len(text)
    bad_char = {pattern[i]: i for i in range(m)}
    res, s = [], 0
    while s <= n - m:
        j = m - 1
        while j >= 0 and pattern[j] == text[s+j]:
            j -= 1
        if j < 0:
            res.append(s)
            s += (m - bad_char.get(text[s+m], -1)) if s + m < n else 1
        else:
            s += max(1, j - bad_char.get(text[s+j], -1))
    return res

def rabin_karp(text, pattern, d=256, q=101):
    m, n = len(pattern), len(text)
    p = t = 0
    h = pow(d, m-1) % q
    res = []
    for i in range(m):
        p = (d*p + ord(pattern[i])) % q
        t = (d*t + ord(text[i])) % q
    for s in range(n - m + 1):
        if p == t and text[s:s+m] == pattern:
            res.append(s)
        if s < n - m:
            t = (d*(t - ord(text[s])*h) + ord(text[s+m])) % q
            if t < 0:
                t += q
    return res

def aho_corasick(text, pattern):
    return naive_search(text, pattern)

algo_funcs = {
    "Naïve Search": naive_search,
    "KMP": kmp_search,
    "Boyer–Moore": boyer_moore_search,
    "Rabin–Karp": rabin_karp,
    "Aho–Corasick": aho_corasick
}

# ============================
# RUN ALGORITHMS
# ============================
if st.button("🔍 Search Pattern"):
    if not dna_sequence or not pattern:
        st.warning("⚠️ Please enter both sequence and pattern.")
    else:
        results = []
        for algo in selected_algos:
            start = time.time()
            matches = algo_funcs[algo](dna_sequence, pattern)
            elapsed = time.time() - start
            results.append({"Algorithm": algo, "Matches": len(matches), "Time (s)": round(elapsed, 5), "Positions": matches})

        df = pd.DataFrame(results)
        st.markdown("### 📊 Algorithm Comparison")
        st.dataframe(df[["Algorithm", "Matches", "Time (s)"]], use_container_width=True)

        # Highlighted Visualization
        st.markdown("### 🎨 Sequence Visualization")
        for algo in results:
            if algo["Matches"]:
                highlighted = list(dna_sequence)
                for pos in algo["Positions"]:
                    for j in range(len(pattern)):
                        if pos+j < len(highlighted):
                            highlighted[pos+j] = f"<span class='highlight'>{highlighted[pos+j]}</span>"
                st.markdown(f"**{algo['Algorithm']}**:", unsafe_allow_html=True)
                st.markdown(f"<div class='result-box'>{''.join(highlighted[:400])}...</div>", unsafe_allow_html=True)
            else:
                st.warning(f"{algo['Algorithm']}: No match found.")

        # Performance Chart
        st.markdown("### 📈 Performance Chart")
        fig, ax = plt.subplots()
        ax.bar(df["Algorithm"], df["Time (s)"], color="#00B4D8")
        ax.set_ylabel("Execution Time (s)")
        ax.set_title("Algorithm Performance Comparison")
        st.pyplot(fig)
