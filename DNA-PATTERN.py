# ==========================
# DNA Pattern Analyzer 🧬 (Single Pattern, Visualization)
# ==========================
import streamlit as st
import re
import time
import pandas as pd

# ==========================
# PAGE CONFIG
# ==========================
st.set_page_config(page_title="DNA Pattern Visualizer", page_icon="🧬", layout="wide")
st.title("🧬 DNA Pattern Analyzer with Visualization")
st.write("Single pattern search: see how each algorithm finds matches in DNA sequence.")

# ==========================
# INPUT
# ==========================
seq_input = st.text_area("Enter DNA Sequence", height=150, placeholder="ATCGGATCGATCG...")
sequence = re.sub(r'[^ATCG]', '', seq_input.upper())

pattern = st.text_input("Enter DNA Pattern to Search", placeholder="ATCG")
pattern = pattern.upper()

algorithms = st.multiselect("Select Algorithm(s)", ["Naive", "KMP", "Boyer-Moore", "Rabin-Karp"], default=["Naive","KMP"])

if not sequence or not pattern or not algorithms:
    st.warning("Please enter sequence, pattern, and select at least one algorithm.")
    st.stop()

# ==========================
# SEARCH ALGORITHMS
# ==========================
def naive_search(text, pattern):
    res=[]
    for i in range(len(text)-len(pattern)+1):
        if text[i:i+len(pattern)]==pattern:
            res.append(i)
    return res

def kmp_search(text, pattern):
    lps=[0]*len(pattern)
    j=0
    for i in range(1,len(pattern)):
        while j>0 and pattern[i]!=pattern[j]: j=lps[j-1]
        if pattern[i]==pattern[j]: j+=1; lps[i]=j
    res=[]; j=0
    for i in range(len(text)):
        while j>0 and text[i]!=pattern[j]: j=lps[j-1]
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
# RUN SEARCH & VISUALIZE
# ==========================
for algo in algorithms:
    st.markdown(f"### 🧪 {algo} Algorithm")
    start = time.time()
    matches = algo_funcs[algo](sequence, pattern)
    elapsed = round(time.time()-start,5)
    st.write(f"Found {len(matches)} matches in {elapsed} seconds.")

    # Highlight matches in sequence (show first 300 bases for readability)
    seq_list = list(sequence[:300])
    for pos in matches:
        for j in range(len(pattern)):
            if pos+j < len(seq_list):
                seq_list[pos+j] = f"<span style='background-color:yellow'>{seq_list[pos+j]}</span>"

    st.markdown("<div style='font-family:monospace; line-height:1.5'>" + "".join(seq_list) + "</div>", unsafe_allow_html=True)
