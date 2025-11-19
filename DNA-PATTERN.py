# ==========================

# DNA Pattern Matching Analyzer 🧬 (Streamlit App)

# Multi-Pattern Aho–Corasick Supported

# ==========================

import streamlit as st
import re, time, io
import pandas as pd
import matplotlib.pyplot as plt
from collections import deque

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

st.markdown(""" <style>
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
.highlight { color: #FFD60A; font-weight: bold; } </style>
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
if sequence:
sequences[header] = sequence
st.success(f"✅ Loaded: {header} ({len(sequence)} bp)")
else:
st.warning(f"⚠️ No valid DNA sequence found in {header}.")
else:
seq_input = st.text_area("🧬 Enter DNA Sequence", placeholder="ATCGGATCGATCG...", height=120).strip().upper()
if seq_input:
cleaned_seq = re.sub(r'[^ATCG]', '', seq_input)
if cleaned_seq:
sequences["Manual Entry"] = cleaned_seq
else:
st.warning("⚠️ Enter a valid DNA sequence containing only A, T, C, G.")

# Allow multiple patterns separated by comma or newline

pattern_input = st.text_area("🔍 Enter Pattern(s) to Search (comma or newline separated)", placeholder="CGATCGA\nATCG")
patterns = [re.sub(r'[^ATCG]', '', p.strip().upper()) for p in re.split('[,\n]', pattern_input) if p.strip()]

algorithms = ["Naïve Search", "KMP", "Boyer–Moore", "Rabin–Karp", "Aho–Corasick"]
selected_algos = st.multiselect(
"⚙️ Select Algorithms",
algorithms,
default=["KMP", "Boyer–Moore"]
)

# ==========================

# ALGORITHMS WITH SAFEGUARDS

# ==========================

def naive_search(text, patterns):
if not text or not patterns:
return {}
res = {}
for pat in patterns:
res[pat] = [i for i in range(len(text)-len(pat)+1) if text[i:i+len(pat)] == pat] if len(text) >= len(pat) else []
return res

def kmp_search(text, patterns):
if not text or not patterns:
return {}
def kmp_single(text, pattern):
if len(text) < len(pattern): return []
lps = [0]*len(pattern); j = 0
for i in range(1,len(pattern)):
while j>0 and pattern[i]!=pattern[j]: j=lps[j-1]
if pattern[i]==pattern[j]: j+=1; lps[i]=j
res,j=[],0
for i in range(len(text)):
while j>0 and text[i]!=pattern[j]: j=lps[j-1]
if text[i]==pattern[j]: j+=1
if j==len(pattern): res.append(i-len(pattern)+1); j=lps[j-1]
return res
return {pat: kmp_single(text, pat) for pat in patterns}

def boyer_moore_search(text, patterns):
if not text or not patterns:
return {}
def bm_single(text, pat):
if len(text)<len(pat): return []
m,n=len(pat),len(text)
bad_char={pat[i]:i for i in range(m)}
res,s=[],0
while s<=n-m:
j=m-1
while j>=0 and pat[j]==text[s+j]: j-=1
if j<0: res.append(s); s += (m - bad_char.get(text[s+m],-1)) if s+m<n else 1
else: s += max(1,j - bad_char.get(text[s+j],-1))
return res
return {pat: bm_single(text, pat) for pat in patterns}

def rabin_karp(text, patterns, d=256,q=101):
if not text or not patterns:
return {}
def rk_single(text, pat):
m,n=len(pat),len(text)
if n<m: return []
h=pow(d,m-1)%q; p=t=0
for i in range(m): p=(d*p+ord(pat[i]))%q; t=(d*t+ord(text[i]))%q
res=[]
for s in range(n-m+1):
if p==t and text[s:s+m]==pat: res.append(s)
if s<n-m: t=(d*(t - ord(text[s])*h)+ord(text[s+m]))%q; t+=q if t<0 else 0
return res
return {pat: rk_single(text, pat) for pat in patterns}

# ==========================

# Multi-Pattern Aho–Corasick

# ==========================

class TrieNode:
def **init**(self):
self.children = {}
self.fail = None
self.output = []

def build_aho_trie(patterns):
root = TrieNode()
for pat in patterns:
node = root
for char in pat:
if char not in node.children: node.children[char] = TrieNode()
node = node.children[char]
node.output.append(pat)
queue=deque()
for child in root.children.values(): child.fail=root; queue.append(child)
while queue:
node=queue.popleft()
for char, child in node.children.items():
fail_node=node.fail
while fail_node and char not in fail_node.children: fail_node=fail_node.fail
child.fail=fail_node.children[char] if fail_node and char in fail_node.children else root
child.output += child.fail.output
queue.append(child)
return root

def aho_corasick(text, patterns):
if not text or not patterns:
return {}
root = build_aho_trie(patterns)
node = root
res = {pat: [] for pat in patterns}
for i,char in enumerate(text):
while node and char not in node.children: node=node.fail
node=node.children[char] if node and char in node.children else root
for pat in node.output: res[pat].append(i-len(pat)+1)
return res

algo_funcs = {
"Naïve Search": naive_search,
"KMP": kmp_search,
"Boyer–Moore": boyer_moore_search,
"Rabin–Karp": rabin_karp,
"Aho–Corasick": aho_corasick
}

# ==========================

# RUN ANALYSIS

# ==========================

if "results_stored" not in st.session_state:
st.session_state.results_stored = None

if st.button("🔍 Search Pattern(s)"):
if not sequences or not patterns:
st.warning("⚠️ Please enter both sequence(s) and pattern(s).")
else:
all_results = []
for header, dna_sequence in sequences.items():
st.markdown(f"## 🧫 Results for **{header}** ({len(dna_sequence)} bp)")

```
        results_list = []
        for algo in selected_algos:
            start = time.time()
            matches_dict = algo_funcs[algo](dna_sequence, patterns)
            elapsed = time.time() - start
            for pat, matches in matches_dict.items():
                results_list.append({
                    "Sequence Name": header,
                    "Pattern": pat,
                    "Algorithm": algo,
                    "Matches": len(matches),
                    "Time (s)": round(elapsed,5)
                })

        df = pd.DataFrame(results_list)
        all_results.append(df)
        st.markdown("### 📊 Algorithm Comparison")
        st.dataframe(df, use_container_width=True)

        # Visualization (first 400 bp)
        st.markdown("### 🎨 Sequence Visualization")
        for algo in selected_algos:
            for pat in patterns:
                matches = algo_funcs[algo](dna_sequence[:400],[pat])[pat]
                if matches:
                    highlighted=list(dna_sequence[:400])
                    for pos in matches:
                        for j in range(len(pat)):
                            if pos+j < len(highlighted): highlighted[pos+j]=f"<span class='highlight'>{highlighted[pos+j]}</span>"
                    st.markdown(f"**{algo}** pattern **{pat}**:", unsafe_allow_html=True)
                    st.markdown(f"<div class='result-box'>{''.join(highlighted)}...</div>", unsafe_allow_html=True)
                else:
                    st.warning(f"{algo} pattern {pat}: No match found.")

        # Performance chart
        st.markdown("### 📈 Performance Chart")
        chart_df=df.groupby("Algorithm")["Time (s)"].max().reset_index()
        fig,ax=plt.subplots()
        ax.bar(chart_df["Algorithm"], chart_df["Time (s)"], color="#00B4D8")
        ax.set_ylabel("Execution Time (s)")
        ax.set_title("Algorithm Performance Comparison")
        st.pyplot(fig)

    # Store combined results
    combined_df=pd.concat(all_results, ignore_index=True)
    st.session_state.results_stored=combined_df
    st.success("✅ Analysis complete! CSV file ready for download!")
```

# ==========================

# DOWNLOAD CSV

# ==========================

if st.session_state.results_stored is not None:
csv_buffer=io.StringIO()
st.session_state.results_stored.to_csv(csv_buffer,index=False)
st.download_button(
label="📥 Download Results as CSV",
data=csv_buffer.getvalue(),
file_name="dna_pattern_results.csv",
mime="text/csv"
)
