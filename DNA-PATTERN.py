# ==========================
# DNA Pattern Analyzer with FASTA Upload and Algorithm Visualization 🧬
# ==========================
import streamlit as st
import re
import time

# ==========================
# PAGE CONFIG
# ==========================
st.set_page_config(page_title="DNA Pattern Analyzer", page_icon="🧬", layout="wide")
st.title("🧬 DNA Pattern Analyzer with Algorithm Visualization")
st.write("Single pattern search: Naive, KMP, Boyer-Moore, Rabin-Karp.")

# ==========================
# FILE UPLOAD / MANUAL INPUT
# ==========================
uploaded_file = st.file_uploader("Upload FASTA file (.fasta, .fa, .txt)", type=["fasta","fa","txt"])
sequence = ""

if uploaded_file:
    content = uploaded_file.getvalue().decode("utf-8").strip()
    lines = content.split("\n")
    seq_list = []
    for line in lines:
        if not line.startswith(">"):
            seq_list.append(line.strip())
    sequence = re.sub(r'[^ATCG]', '', "".join(seq_list).upper())
else:
    seq_input = st.text_area("Enter DNA Sequence manually", height=150)
    if seq_input:
        sequence = re.sub(r'[^ATCG]', '', seq_input.upper())

if not sequence:
    st.warning("Please enter sequence manually or upload a FASTA file.")
    st.stop()

# ==========================
# PATTERN INPUT
# ==========================
pattern = st.text_input("Enter DNA Pattern to Search", placeholder="ATCG").upper()
if not pattern:
    st.warning("Please enter a pattern to search.")
    st.stop()

algorithms = st.multiselect("Select algorithm(s)", ["Naive","KMP","Boyer-Moore","Rabin-Karp"], default=["Naive","KMP"])

if not algorithms:
    st.warning("Select at least one algorithm.")
    st.stop()

# ==========================
# SEARCH ALGORITHMS
# ==========================
def naive_search(text, pattern):
    return [i for i in range(len(text)-len(pattern)+1) if text[i:i+len(pattern)]==pattern]

def compute_lps(pattern):
    lps = [0]*len(pattern)
    j = 0
    for i in range(1,len(pattern)):
        while j>0 and pattern[i]!=pattern[j]: j=lps[j-1]
        if pattern[i]==pattern[j]: j+=1; lps[i]=j
    return lps

def kmp_search(text, pattern):
    lps = compute_lps(pattern)
    res=[]
    i=j=0
    while i<len(text):
        if pattern[j]==text[i]:
            i+=1; j+=1
        if j==len(pattern):
            res.append(i-j)
            j=lps[j-1]
        elif i<len(text) and pattern[j]!=text[i]:
            if j!=0: j=lps[j-1]
            else: i+=1
    return res, lps

def boyer_moore_search(text, pattern):
    bad = {c:i for i,c in enumerate(pattern)}
    res=[]; s=0; m=len(pattern); n=len(text)
    while s<=n-m:
        j=m-1
        while j>=0 and text[s+j]==pattern[j]: j-=1
        if j<0: res.append(s); s+=1
        else: s+=max(1,j-bad.get(text[s+j],-1))
    return res, bad

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
    
    if algo=="KMP":
        matches,lps = kmp_search(sequence, pattern)
    elif algo=="Boyer-Moore":
        matches,bad = boyer_moore_search(sequence, pattern)
    else:
        matches = algo_funcs[algo](sequence, pattern)
    
    elapsed = round(time.time()-start,5)
    st.write(f"Found {len(matches)} matches in {elapsed} seconds.")
    
    # Highlight matches in sequence (first 300 bases)
    seq_list = list(sequence[:300])
    for pos in matches:
        for j in range(len(pattern)):
            if pos+j < len(seq_list):
                seq_list[pos+j] = f"<span style='background-color:yellow'>{seq_list[pos+j]}</span>"
    st.markdown("<div style='font-family:monospace; line-height:1.5'>" + "".join(seq_list) + "</div>", unsafe_allow_html=True)

    # Show KMP LPS table
    if algo=="KMP":
        st.markdown("**LPS Table:**")
        st.write(lps)

    # Show Boyer-Moore bad character table
    if algo=="Boyer-Moore":
        st.markdown("**Bad Character Table:**")
        st.write(bad)
