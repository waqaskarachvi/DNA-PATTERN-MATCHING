# ==========================
# DNA Pattern Matching Analyzer 🧬 (Streamlit App)
# ==========================
import streamlit as st
import re, time, io
import pandas as pd
import matplotlib.pyplot as plt
from collections import deque, defaultdict

# ==========================
# PAGE CONFIG
# ==========================
st.set_page_config(page_title="DNA Pattern Matching Analyzer", page_icon="🧬", layout="wide")

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
        .result-box { background-color: #1E2636; padding: 15px; border-radius: 10px; border: 1px solid #00B4D8;
                      font-family: monospace; line-height: 1.6; word-wrap: break-word; }
        .highlight { color: #FFD60A; font-weight: bold; }
    </style>
""", unsafe_allow_html=True)

# ==========================
# HEADER
# ==========================
st.markdown("<h1 style='text-align:center;'>🧬 DNA Pattern Matching Analyzer</h1>", unsafe_allow_html=True)
st.markdown("<p style='text-align:center;'>Compare string matching algorithms for DNA sequences</p>", unsafe_allow_html=True)
st.markdown("---", unsafe_allow_html=True)

# ==========================
# FILE UPLOAD / INPUT
# ==========================
uploaded_files = st.file_uploader(
    "📁 Upload FASTA files (you can select multiple)",
    type=["fasta","fa","txt"],
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
    seq_input = st.text_area("🧬 Enter DNA Sequence", placeholder="ATCGGATCGATCG...", height=120).strip().upper()
    if seq_input:
        sequences["Manual Entry"] = re.sub(r'[^ATCG]', '', seq_input)

pattern_input = st.text_input(
    "🔍 Enter Pattern(s) to Search (comma separated for multiple)",
    placeholder="CGATCGA,ATGCGT"
).strip().upper()

algorithms = ["Naïve Search", "KMP", "Boyer–Moore", "Rabin–Karp", "Aho–Corasick"]
selected_algos = st.multiselect("⚙️ Select Algorithms", algorithms, default=["KMP","Boyer–Moore"])

# ==========================
# ALGORITHMS
# ==========================
def naive_search(text, pattern):
    return [i for i in range(len(text)-len(pattern)+1) if text[i:i+len(pattern)] == pattern]

def kmp_search(text, pattern):
    lps=[0]*len(pattern); j=0
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
    m,n=len(pattern),len(text)
    bad_char={pattern[i]:i for i in range(m)}
    res=[]; s=0
    while s<=n-m:
        j=m-1
        while j>=0 and pattern[j]==text[s+j]: j-=1
        if j<0: res.append(s); s+=(m-bad_char.get(text[s+m],-1)) if s+m<n else 1
        else: s+=max(1,j-bad_char.get(text[s+j],-1))
    return res

def rabin_karp(text, pattern, d=256, q=101):
    m,n=len(pattern),len(text); p=t=0; h=pow(d,m-1)%q; res=[]
    for i in range(m): p=(d*p+ord(pattern[i]))%q; t=(d*t+ord(text[i]))%q
    for s in range(n-m+1):
        if p==t and text[s:s+m]==pattern: res.append(s)
        if s<n-m: t=(d*(t-ord(text[s])*h)+ord(text[s+m]))%q; t+=q if t<0 else 0
    return res

# ==========================
# AHO-CORASICK
# ==========================
class AhoNode:
    def __init__(self): self.children={}; self.fail=None; self.output=[]
class AhoCorasick:
    def __init__(self,patterns):
        self.root=AhoNode(); self.build_trie(patterns); self.build_failure_links()
    def build_trie(self,patterns):
        for pat in patterns:
            node=self.root
            for c in pat:
                if c not in node.children: node.children[c]=AhoNode()
                node=node.children[c]
            node.output.append(pat)
    def build_failure_links(self):
        queue=deque()
        for child in self.root.children.values(): child.fail=self.root; queue.append(child)
        while queue:
            current=queue.popleft()
            for c,child in current.children.items():
                f=current.fail
                while f and c not in f.children: f=f.fail
                child.fail=f.children[c] if f and c in f.children else self.root
                child.output+=child.fail.output
                queue.append(child)
    def search(self,text):
        node=self.root
        res=defaultdict(list)
        for i,c in enumerate(text):
            while node and c not in node.children: node=node.fail
            if not node: node=self.root; continue
            node=node.children[c]
            for pat in node.output: res[pat].append(i-len(pat)+1)
        return res

# ==========================
# RUN ANALYSIS
# ==========================
if "results_stored" not in st.session_state: st.session_state.results_stored=None

if st.button("🔍 Search Pattern"):
    if not sequences or not pattern_input:
        st.warning("⚠️ Please enter sequence(s) and pattern(s).")
    else:
        patterns_list = [p.strip() for p in pattern_input.split(",") if p.strip()]
        all_results=[]
        for header,dna_sequence in sequences.items():
            st.markdown(f"## 🧫 Results for **{header}** ({len(dna_sequence)} bp)")
            results=[]
            for algo in selected_algos:
                if algo=="Aho–Corasick":
                    if len(patterns_list)<2: st.warning("⚠️ Aho–Corasick requires multiple patterns."); continue
                    start=time.time()
                    matches=AhoCorasick(patterns_list).search(dna_sequence)
                    elapsed=time.time()-start
                    for pat,pos_list in matches.items():
                        results.append({"Sequence Name":header,"Algorithm":f"Aho–Corasick ({pat})","Matches":len(pos_list),"Time (s)":round(elapsed,5)})
                else:
                    start=time.time()
                    # Run single-pattern algorithms sequentially for multiple patterns
                    total_matches=0
                    for pat in patterns_list:
                        total_matches+=len({"Naïve Search": naive_search,
                                            "KMP": kmp_search,
                                            "Boyer–Moore": boyer_moore_search,
                                            "Rabin–Karp": rabin_karp}[algo](dna_sequence, pat))
                    elapsed=time.time()-start
                    results.append({"Sequence Name":header,"Algorithm":algo,"Matches":total_matches,"Time (s)":round(elapsed,5)})
            df=pd.DataFrame(results)
            all_results.append(df)
            st.markdown("### 📊 Algorithm Comparison")
            st.dataframe(df,use_container_width=True)
            st.markdown("### 📈 Performance Chart")
            fig,ax=plt.subplots(); ax.bar(df["Algorithm"],df["Time (s)"],color="#00B4D8"); ax.set_ylabel("Time (s)"); ax.set_title("Algorithm Performance")
            st.pyplot(fig)

        combined_df=pd.concat(all_results,ignore_index=True)
        st.session_state.results_stored=combined_df

# ==========================
# DOWNLOAD CSV
# ==========================
if st.session_state.results_stored is not None:
    csv_buffer=io.StringIO()
    st.session_state.results_stored.to_csv(csv_buffer,index=False)
    st.download_button(label="📥 Download Results as CSV",data=csv_buffer.getvalue(),file_name="dna_pattern_results.csv",mime="text/csv")
