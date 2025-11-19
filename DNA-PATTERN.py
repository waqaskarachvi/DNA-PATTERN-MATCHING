# ==========================
# DNA Pattern Analyzer 🧬 (Single Pattern, Multiple Algorithms)
# ==========================
import streamlit as st
import re
import time
import pandas as pd
import matplotlib.pyplot as plt

# ==========================
# PAGE CONFIG
# ==========================
st.set_page_config(page_title="DNA Pattern Analyzer", page_icon="🧬", layout="wide")
st.title("🧬 DNA Pattern Analyzer")
st.write("Single Pattern Search with multiple algorithms and comparison chart.")

# ==========================
# FILE UPLOAD / MANUAL INPUT
# ==========================
uploaded_file = st.file_uploader("Upload FASTA file (.fasta, .fa, .txt)", type=["fasta","fa","txt"])
sequences = {}

if uploaded_file:
    content = uploaded_file.getvalue().decode("utf-8").strip()
    lines = content.split("\n")
    header = None
    seq = ""
    for line in lines:
        line = line.strip()
        if line.startswith(">"):
            if header and seq:
                sequences[header] = re.sub(r'[^ATCG]', '', seq.upper())
            header = line[1:]
            seq = ""
        else:
            seq += line
    if header and seq:
        sequences[header] = re.sub(r'[^ATCG]', '', seq.upper())
    for name, s in sequences.items():
        st.success(f"✅ Loaded: {name} ({len(s)} bp)")
else:
    seq_input = st.text_area("Enter DNA Sequence manually", height=120)
    if seq_input:
        sequences["Manual Entry"] = re.sub(r'[^ATCG]', '', seq_input.upper())

if not sequences:
    st.warning("Please upload a FASTA file or enter a DNA sequence manually.")
    st.stop()

# ==========================
# PATTERN INPUT
# ==========================
pattern = st.text_input("Enter DNA pattern to search").strip().upper()
if not pattern:
    st.warning("Please enter a DNA pattern to search.")
    st.stop()

# ==========================
# SELECT ALGORITHMS
# ==========================
algorithms = st.multiselect(
    "Select Algorithm(s) to use",
    ["Naive", "KMP", "Boyer-Moore", "Rabin-Karp"],
    default=["Naive", "KMP"]
)

if not algorithms:
    st.warning("Please select at least one algorithm.")
    st.stop()

# ==========================
# SEARCH ALGORITHMS
# ==========================
def naive_search(text, pattern):
    return [i for i in range(len(text)-len(pattern)+1) if text[i:i+len(pattern)]==pattern]

def kmp_search(text, pattern):
    lps = [0]*len(pattern)
    j=0
    for i in range(1,len(pattern)):
        while j>0 and pattern[i]!=pattern[j]:
            j=lps[j-1]
        if pattern[i]==pattern[j]: j+=1; lps[i]=j
    res=[]; j=0
    for i in range(len(text)):
        while j>0 and text[i]!=pattern[j]:
            j=lps[j-1]
        if text[i]==pattern[j]: j+=1
        if j==len(pattern): res.append(i-j+1); j=lps[j-1]
    return res

def boyer_moore_search(text, pattern):
    bad = {c:i for i,c in enumerate(pattern)}
    res=[]; s=0; m=len(pattern); n=len(text)
    while s<=n-m:
        j=m-1
        while j>=0 and text[s+j]==pattern[j]: j-=1
        if j<0: res.append(s); s+=1
        else: s+=max(1,j-bad.get(text[s+j],-1))
    return res

def rabin_karp_search(text, pattern):
    d,q=256,101; m,n=len(pattern),len(text)
    p=t=0; h=pow(d,m-1)%q; res=[]
    for i in range(m): p=(d*p+ord(pattern[i]))%q; t=(d*t+ord(text[i]))%q
    for s in range(n-m+1):
        if p==t and text[s:s+m]==pattern: res.append(s)
        if s<n-m: t=(d*(t-ord(text[s])*h)+ord(text[s+m]))%q; t+=q if t<0 else 0
    return res

algo_funcs = {
    "Naive": naive_search,
    "KMP": kmp_search,
    "Boyer-Moore": boyer_moore_search,
    "Rabin-Karp": rabin_karp_search
}

# ==========================
# RUN SEARCH
# ==========================
if st.button("🔍 Search Pattern"):
    all_results = []

    for header, dna_sequence in sequences.items():
        st.markdown(f"### 🧫 Results for **{header}** ({len(dna_sequence)} bp)")
        results=[]

        for algo in algorithms:
            start = time.time()
            matches = algo_funcs[algo](dna_sequence, pattern)
            elapsed = round(time.time()-start,5)
            results.append({"Algorithm": algo, "Matches": len(matches), "Time (s)": elapsed})
            st.write(f"**{algo}** found {len(matches)} matches in {elapsed} seconds.")

        # Display comparison table
        df = pd.DataFrame(results)
        st.dataframe(df)

        # Performance chart
        fig, ax = plt.subplots()
        ax.bar(df["Algorithm"], df["Time (s)"], color="#00B4D8")
        ax.set_ylabel("Execution Time (s)")
        ax.set_title("Algorithm Performance Comparison")
        st.pyplot(fig)
