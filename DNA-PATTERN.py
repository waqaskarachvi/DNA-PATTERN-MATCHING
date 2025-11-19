# ==========================
# DNA Pattern Analyzer 🧬 with Step-by-Step Visualization & Rabin-Karp Hash
# ==========================
import streamlit as st
import re
import time

# ==========================
# PAGE CONFIG
# ==========================
st.set_page_config(page_title="DNA Pattern Step Visualizer", page_icon="🧬", layout="wide")
st.title("🧬 DNA Pattern Analyzer with Step-by-Step Visualization")
st.write("See how Naive, KMP, Boyer-Moore, and Rabin-Karp scan DNA sequences with rolling hash visualization.")

# ==========================
# INPUT
# ==========================
uploaded_file = st.file_uploader("Upload FASTA file (.fasta/.fa/.txt)", type=["fasta","fa","txt"])
sequence = ""

if uploaded_file:
    content = uploaded_file.getvalue().decode("utf-8")
    lines = content.split("\n")
    seq_list = [line.strip() for line in lines if not line.startswith(">")]
    sequence = re.sub(r'[^ATCG]', '', "".join(seq_list).upper())
else:
    seq_input = st.text_area("Enter DNA Sequence manually", height=150)
    if seq_input:
        sequence = re.sub(r'[^ATCG]', '', seq_input.upper())

if not sequence:
    st.warning("Please enter sequence or upload a FASTA file.")
    st.stop()

pattern = st.text_input("Enter DNA Pattern to Search").upper()
if not pattern:
    st.warning("Please enter a pattern.")
    st.stop()

algorithms = st.multiselect("Select Algorithm(s)", ["Naive","KMP","Boyer-Moore","Rabin-Karp"], default=["Naive","KMP"])

if not algorithms:
    st.warning("Select at least one algorithm.")
    st.stop()

# ==========================
# SEARCH ALGORITHMS
# ==========================
def naive_search_steps(text, pattern):
    steps=[]
    matches=[]
    n,m=len(text),len(pattern)
    for i in range(n-m+1):
        step_text=list(text[:300])
        for j in range(m):
            step_text[i+j]=f"<span style='background-color:yellow'>{step_text[i+j]}</span>"
        steps.append("".join(step_text))
        if text[i:i+m]==pattern:
            matches.append(i)
    return matches, steps

def compute_lps(pattern):
    lps=[0]*len(pattern)
    j=0
    for i in range(1,len(pattern)):
        while j>0 and pattern[i]!=pattern[j]:
            j=lps[j-1]
        if pattern[i]==pattern[j]:
            j+=1; lps[i]=j
    return lps

def kmp_search_steps(text, pattern):
    lps=compute_lps(pattern)
    steps=[]
    matches=[]
    i=j=0
    while i<len(text):
        step_text=list(text[:300])
        if j<len(pattern):
            step_text[i]=f"<span style='background-color:yellow'>{step_text[i]}</span>"
        steps.append("".join(step_text))
        if j<len(pattern) and pattern[j]==text[i]:
            i+=1; j+=1
        if j==len(pattern):
            matches.append(i-j)
            j=lps[j-1]
        elif i<len(text) and j<len(pattern) and pattern[j]!=text[i]:
            if j!=0: j=lps[j-1]
            else: i+=1
    return matches, steps, lps

def boyer_moore_search_steps(text, pattern):
    bad={c:i for i,c in enumerate(pattern)}
    steps=[]
    matches=[]
    s=0; m=len(pattern); n=len(text)
    while s<=n-m:
        j=m-1
        while j>=0 and text[s+j]==pattern[j]: j-=1
        step_text=list(text[:300])
        for k in range(m):
            if s+k<len(step_text):
                step_text[s+k]=f"<span style='background-color:yellow'>{step_text[s+k]}</span>"
        steps.append("".join(step_text))
        if j<0: matches.append(s); s+=1
        else: s+=max(1,j-bad.get(text[s+j],-1))
    return matches, steps, bad

def rabin_karp_search_steps(text, pattern):
    d,q=256,101; m,n=len(pattern),len(text)
    p=t=0; h=pow(d,m-1)%q; matches=[]; steps=[]; hashes=[]
    for i in range(m): 
        p=(d*p+ord(pattern[i]))%q
        t=(d*t+ord(text[i]))%q
    for s in range(n-m+1):
        step_text=list(text[:300])
        for k in range(m):
            if s+k<len(step_text):
                step_text[s+k]=f"<span style='background-color:yellow'>{step_text[s+k]}</span>"
        steps.append("".join(step_text))
        hashes.append(t)
        if p==t and text[s:s+m]==pattern: matches.append(s)
        if s<n-m: t=(d*(t-ord(text[s])*h)+ord(text[s+m]))%q; t+=q if t<0 else 0
    return matches, steps, hashes

algo_funcs={
    "Naive": naive_search_steps,
    "KMP": kmp_search_steps,
    "Boyer-Moore": boyer_moore_search_steps,
    "Rabin-Karp": rabin_karp_search_steps
}

# ==========================
# RUN & VISUALIZE
# ==========================
for algo in algorithms:
    st.markdown(f"### 🧪 {algo} Algorithm")
    start=time.time()
    if algo=="KMP":
        matches, steps, lps = kmp_search_steps(sequence, pattern)
    elif algo=="Boyer-Moore":
        matches, steps, bad = boyer_moore_search_steps(sequence, pattern)
    elif algo=="Rabin-Karp":
        matches, steps, hashes = rabin_karp_search_steps(sequence, pattern)
    else:
        matches, steps = algo_funcs[algo](sequence, pattern)
    elapsed=round(time.time()-start,5)
    st.write(f"Found {len(matches)} matches in {elapsed} seconds.")

    # Show steps
    step_index=st.slider("Step", 0, len(steps)-1, 0)
    st.markdown("<div style='font-family:monospace; line-height:1.5'>" + steps[step_index] + "</div>", unsafe_allow_html=True)

    # Show LPS / Bad / Hashes table
    if algo=="KMP":
        st.markdown("**LPS Table:**")
        st.write(lps)
    if algo=="Boyer-Moore":
        st.markdown("**Bad Character Table:**")
        st.write(bad)
    if algo=="Rabin-Karp":
        st.markdown("**Rolling Hash Values (q=101)**")
        st.write(hashes)
