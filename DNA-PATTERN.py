# ==========================
# DNA Pattern Matching Analyzer 🧬
# ==========================
import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import time
from collections import deque, defaultdict

# ==========================
# PAGE CONFIG
# ==========================
st.set_page_config(
    page_title="DNA Pattern Matching Analyzer",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 DNA Pattern Matching Analyzer")
st.write("Search DNA sequences using multiple algorithms. Single pattern = multiple algorithms, Multi-pattern = Aho–Corasick.")

# ==========================
# FILE UPLOAD / INPUT
# ==========================
uploaded_file = st.file_uploader("Upload DNA sequence file", type=["txt"])
seq_input = st.text_area("Or enter DNA sequence manually", height=100)

sequence = ""
if uploaded_file:
    sequence = uploaded_file.read().decode("utf-8").replace("\n","").upper()
elif seq_input:
    sequence = seq_input.replace("\n","").upper()

if not sequence:
    st.warning("Please provide a DNA sequence.")
    st.stop()

# ==========================
# PATTERN INPUT
# ==========================
patterns_raw = st.text_area("Enter pattern(s) (one per line)")
patterns = [p.strip().upper() for p in patterns_raw.splitlines() if p.strip()]

if not patterns:
    st.warning("Please enter at least one pattern.")
    st.stop()

mode = "Single" if len(patterns) == 1 else "Multi"

# ==========================
# ALGORITHMS
# ==========================
def naive_search(text, pattern):
    return [i for i in range(len(text)-len(pattern)+1) if text[i:i+len(pattern)] == pattern]

def kmp_search(text, pattern):
    lps = [0]*len(pattern)
    j = 0
    for i in range(1,len(pattern)):
        while j>0 and pattern[i]!=pattern[j]:
            j = lps[j-1]
        if pattern[i]==pattern[j]:
            j+=1; lps[i]=j
    res=[]; j=0
    for i in range(len(text)):
        while j>0 and text[i]!=pattern[j]:
            j=lps[j-1]
        if text[i]==pattern[j]:
            j+=1
        if j==len(pattern):
            res.append(i-j+1)
            j=lps[j-1]
    return res

def boyer_moore_search(text, pattern):
    bad = {c:i for i,c in enumerate(pattern)}
    res=[]; s=0; m=len(pattern); n=len(text)
    while s<=n-m:
        j=m-1
        while j>=0 and text[s+j]==pattern[j]:
            j-=1
        if j<0:
            res.append(s)
            s+=1
        else:
            s+=max(1,j-bad.get(text[s+j],-1))
    return res

def rabin_karp_search(text, pattern):
    d,q=256,101
    m,n=len(pattern),len(text)
    p=t=0; h=pow(d,m-1)%q
    res=[]
    for i in range(m):
        p=(d*p+ord(pattern[i]))%q
        t=(d*t+ord(text[i]))%q
    for s in range(n-m+1):
        if p==t and text[s:s+m]==pattern:
            res.append(s)
        if s<n-m:
            t=(d*(t-ord(text[s])*h)+ord(text[s+m]))%q
            if t<0: t+=q
    return res

# ==========================
# AHO-CORASICK FOR MULTI-PATTERN
# ==========================
class Node:
    def __init__(self):
        self.children = {}
        self.fail = None
        self.output = []

class AhoCorasick:
    def __init__(self, patterns):
        self.root = Node()
        self.build_trie(patterns)
        self.build_failures()

    def build_trie(self, patterns):
        for pat in patterns:
            node=self.root
            for c in pat:
                if c not in node.children:
                    node.children[c]=Node()
                node=node.children[c]
            node.output.append(pat)

    def build_failures(self):
        queue=deque()
        for child in self.root.children.values():
            child.fail=self.root
            queue.append(child)
        while queue:
            current=queue.popleft()
            for c,next_node in current.children.items():
                fail=current.fail
                while fail and c not in fail.children:
                    fail=fail.fail
                next_node.fail=fail.children[c] if fail and c in fail.children else self.root
                next_node.output += next_node.fail.output
                queue.append(next_node)

    def search(self, text):
        node=self.root
        results=defaultdict(list)
        for i,c in enumerate(text):
            while node and c not in node.children:
                node=node.fail
            if not node:
                node=self.root
                continue
            node=node.children[c]
            for pat in node.output:
                results[pat].append(i-len(pat)+1)
        return results

# ==========================
# SELECT ALGORITHMS
# ==========================
if mode=="Single":
    selected_algos = st.multiselect("Select algorithms", ["Naive","KMP","Boyer-Moore","Rabin-Karp"], default=["KMP","Boyer-Moore"])
else:
    selected_algos = ["Aho–Corasick"]

# ==========================
# RUN SEARCH
# ==========================
if st.button("Search"):
    if mode=="Single":
        pattern = patterns[0]
        results=[]
        for algo in selected_algos:
            start=time.time()
            if algo=="Naive": matches=naive_search(sequence,pattern)
            elif algo=="KMP": matches=kmp_search(sequence,pattern)
            elif algo=="Boyer-Moore": matches=boyer_moore_search(sequence,pattern)
            else: matches=rabin_karp_search(sequence,pattern)
            elapsed=time.time()-start
            results.append({"Algorithm":algo,"Matches":len(matches),"Time (s)":round(elapsed,5)})
        df=pd.DataFrame(results)
        st.dataframe(df)

        # Performance chart
        fig,ax=plt.subplots()
        ax.bar(df["Algorithm"], df["Time (s)"], color="#00B4D8")
        ax.set_ylabel("Time (s)")
        ax.set_title("Algorithm Performance")
        st.pyplot(fig)

    else:  # Multi-pattern
        aho = AhoCorasick(patterns)
        results = aho.search(sequence)
        if results:
            df=pd.DataFrame([(pat,pos) for pat, positions in results.items() for pos in positions], columns=["Pattern","Position"])
            st.write(df)
            csv = df.to_csv(index=False).encode("utf-8")
            st.download_button("Download CSV", csv, "aho_results.csv","text/csv")
        else:
            st.info("No matches found.")
