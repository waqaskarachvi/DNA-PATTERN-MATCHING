import streamlit as st
import pandas as pd
from collections import deque, defaultdict
import time
import matplotlib.pyplot as plt

# ---------------- SEARCH ALGORITHMS ----------------

def naive_search(text, pattern):
    return [i for i in range(len(text)-len(pattern)+1) if text[i:i+len(pattern)] == pattern]

def compute_lps(pattern):
    lps = [0]*len(pattern)
    j = 0
    i = 1
    while i < len(pattern):
        if pattern[i]==pattern[j]:
            j+=1
            lps[i]=j
            i+=1
        else:
            if j!=0:
                j=lps[j-1]
            else:
                lps[i]=0
                i+=1
    return lps

def kmp_search(text, pattern):
    lps=compute_lps(pattern)
    occurrences=[]
    i=j=0
    while i<len(text):
        if pattern[j]==text[i]:
            i+=1
            j+=1
        if j==len(pattern):
            occurrences.append(i-j)
            j=lps[j-1]
        elif i<len(text) and pattern[j]!=text[i]:
            if j!=0:
                j=lps[j-1]
            else:
                i+=1
    return occurrences

def bad_character_table(pattern):
    return {c:i for i,c in enumerate(pattern)}

def boyer_moore_search(text, pattern):
    occurrences=[]
    bad_char=bad_character_table(pattern)
    n,m=len(text),len(pattern)
    s=0
    while s<=n-m:
        j=m-1
        while j>=0 and pattern[j]==text[s+j]:
            j-=1
        if j<0:
            occurrences.append(s)
            s += m-bad_char.get(text[s+m],-1) if s+m<n else 1
        else:
            s += max(1,j-bad_char.get(text[s+j],-1))
    return occurrences

def rabin_karp_search(text, pattern,d=256,q=101):
    n,m=len(text),len(pattern)
    occurrences=[]
    h=pow(d,m-1)%q
    p=t=0
    for i in range(m):
        p=(d*p+ord(pattern[i]))%q
        t=(d*t+ord(text[i]))%q
    for s in range(n-m+1):
        if p==t and text[s:s+m]==pattern:
            occurrences.append(s)
        if s<n-m:
            t=(d*(t-ord(text[s])*h)+ord(text[s+m]))%q
            if t<0:
                t+=q
    return occurrences

# ---------------- AHO-CORASICK ----------------

class AhoNode:
    def __init__(self):
        self.children={}
        self.fail=None
        self.output=[]

class AhoCorasick:
    def __init__(self,patterns):
        self.root=AhoNode()
        self.build_trie(patterns)
        self.build_failure_links()
    def build_trie(self,patterns):
        for idx,pat in enumerate(patterns):
            node=self.root
            for c in pat:
                if c not in node.children:
                    node.children[c]=AhoNode()
                node=node.children[c]
            node.output.append((idx,pat))
    def build_failure_links(self):
        queue=deque()
        for child in self.root.children.values():
            child.fail=self.root
            queue.append(child)
        while queue:
            current=queue.popleft()
            for c,next_node in current.children.items():
                fail_node=current.fail
                while fail_node and c not in fail_node.children:
                    fail_node=fail_node.fail
                next_node.fail=fail_node.children[c] if fail_node and c in fail_node.children else self.root
                if next_node.fail and next_node.fail.output:
                    next_node.output+=next_node.fail.output
                queue.append(next_node)
    def search(self,text):
        node=self.root
        results=defaultdict(list)
        for i,c in enumerate(text):
            while node and c not in node.children:
                node=node.fail
            if not node:
                node=self.root
                continue
            node=node.children[c]
            for _,pat in node.output:
                results[pat].append(i-len(pat)+1)
        return results

# ---------------- STREAMLIT UI ----------------

st.set_page_config(page_title="DNA Pattern Search Tool",layout="wide")
st.title("🔬 DNA Pattern Search Tool")
st.write("Single Pattern: select multiple algorithms to compare. Multi-Pattern: Aho–Corasick.")

uploaded_file = st.file_uploader("Upload DNA file (.txt/.fa/.fasta)")
manual_seq = st.text_area("Or paste DNA sequence (A/T/C/G only)", height=150)

sequence = ""
if manual_seq.strip():
    sequence = "".join([c for c in manual_seq.upper() if c in "ATCG"])
elif uploaded_file:
    raw = uploaded_file.read().decode("utf-8")
    lines = [l.strip() for l in raw.splitlines() if not l.startswith(">")]
    sequence = "".join(lines).upper()
    sequence = "".join([c for c in sequence if c in "ATCG"])

if not sequence:
    st.info("No DNA sequence provided yet. Paste or upload to start searching.")

mode = st.radio("Select Search Mode", ["Single Pattern Search", "Multi-Pattern Search"], index=0)

# ---------------- SINGLE PATTERN ----------------
if mode == "Single Pattern Search":
    pattern_raw = st.text_area("Enter single pattern (A/T/C/G)", height=50)
    algo_options = st.multiselect("Select algorithms", ["Naive", "KMP", "Boyer-Moore", "Rabin-Karp"], default=["KMP","Boyer-Moore"])

    if st.button("Search Single Pattern"):
        if not sequence:
            st.warning("Provide a DNA sequence first!")
        elif not pattern_raw.strip():
            st.warning("Enter a pattern to search.")
        else:
            patterns = [p.strip() for p in pattern_raw.strip().splitlines() if p.strip()]
            if len(patterns) != 1:
                st.error("Single Pattern Search mode accepts only ONE pattern. Please enter a single line pattern.")
            else:
                pat = "".join([c for c in patterns[0].upper() if c in "ATCG"])
                if not pat:
                    st.error("Pattern contains no valid A/T/C/G characters.")
                else:
                    algo_funcs = {
                        "Naive": naive_search,
                        "KMP": kmp_search,
                        "Boyer-Moore": boyer_moore_search,
                        "Rabin-Karp": rabin_karp_search
                    }
                    results_list = []
                    for algo in algo_options:
                        start = time.time()
                        found = algo_funcs[algo](sequence, pat)
                        elapsed = time.time() - start
                        results_list.append({"Algorithm": algo, "Matches": len(found), "Time(s)": round(elapsed,5), "Positions": found})

                    df = pd.DataFrame(results_list)
                    st.markdown("### 📊 Algorithm Comparison")
                    st.dataframe(df[["Algorithm","Matches","Time(s)"]])

                    fig, ax = plt.subplots()
                    ax.bar(df["Algorithm"], df["Time(s)"], color="#00B4D8")
                    ax.set_ylabel("Time (s)")
                    ax.set_title("Execution Time Comparison")
                    st.pyplot(fig)

                    st.markdown("### 🎨 Sequence Snippet Highlight")
                    for row in results_list:
                        positions = row["Positions"]
                        if positions:
                            snippet = list(sequence[:400])
                            for pos in positions:
                                if pos < 400:
                                    for j in range(len(pat)):
                                        if pos+j < 400:
                                            snippet[pos+j] = f"**{snippet[pos+j]}**"
                            st.markdown(f"**{row['Algorithm']}**: "+"".join(snippet))
                        else:
                            st.info(f"{row['Algorithm']}: No matches found.")

# ---------------- MULTI-PATTERN ----------------
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
                st.error("Enter at least 2 valid patterns for multi-pattern search.")
            else:
                aho = AhoCorasick(patterns)
                results = aho.search(sequence)
                if not results:
                    st.info("No matches found.")
                else:
                    rows = []
                    for pat,pos_list in results.items():
                        for pos in pos_list:
                            rows.append((pat,pos))
                    df = pd.DataFrame(rows,columns=["Pattern","Position"])
                    st.write(df)
                    csv = df.to_csv(index=False).encode("utf-8")
                    st.download_button("Download CSV",csv,"aho_results.csv","text/csv")
