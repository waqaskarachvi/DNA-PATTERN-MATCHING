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
selected_algos = st.multiselect("⚙️ Select Algorithms", algorithms, default=algorithms)

# ==========================
# ALGORITHMS (with comparison counting)
# ==========================
def naive_search(text, pattern):
    comparisons = 0
    results = []
    for i in range(len(text)-len(pattern)+1):
        comparisons += 1
        match = True
        for j in range(len(pattern)):
            comparisons += 1
            if text[i+j] != pattern[j]:
                match = False
                break
        if match:
            results.append(i)
    return results, comparisons

def kmp_search(text, pattern):
    comparisons = 0
    lps=[0]*len(pattern); j=0
    for i in range(1,len(pattern)):
        while j>0 and pattern[i]!=pattern[j]: 
            j=lps[j-1]
            comparisons += 1
        comparisons += 1
        if pattern[i]==pattern[j]: j+=1; lps[i]=j
    res=[]; j=0
    for i in range(len(text)):
        while j>0 and text[i]!=pattern[j]: 
            j=lps[j-1]
            comparisons += 1
        comparisons += 1
        if text[i]==pattern[j]: j+=1
        if j==len(pattern): res.append(i-j+1); j=lps[j-1]
    return res, comparisons

def boyer_moore_search(text, pattern):
    comparisons = 0
    m,n=len(pattern),len(text)
    bad_char={pattern[i]:i for i in range(m)}
    res=[]; s=0
    while s<=n-m:
        j=m-1
        while j>=0:
            comparisons += 1
            if pattern[j]!=text[s+j]:
                break
            j-=1
        if j<0: 
            res.append(s)
            s+=(m-bad_char.get(text[s+m],-1)) if s+m<n else 1
        else: 
            s+=max(1,j-bad_char.get(text[s+j],-1))
    return res, comparisons

def rabin_karp(text, pattern, d=256, q=101):
    comparisons = 0
    m,n=len(pattern),len(text); p=t=0; h=pow(d,m-1)%q; res=[]
    for i in range(m): p=(d*p+ord(pattern[i]))%q; t=(d*t+ord(text[i]))%q
    for s in range(n-m+1):
        comparisons += 1
        if p==t:
            if text[s:s+m]==pattern: 
                res.append(s)
                comparisons += m
        if s<n-m: t=(d*(t-ord(text[s])*h)+ord(text[s+m]))%q; t+=q if t<0 else 0
    return res, comparisons

# ==========================
# AHO-CORASICK (with comparison counting)
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
        comparisons = 0
        for i,c in enumerate(text):
            comparisons += 1
            while node and c not in node.children: 
                node=node.fail
                comparisons += 1
            if not node: node=self.root; continue
            node=node.children[c]
            for pat in node.output: res[pat].append(i-len(pat)+1)
        return res, comparisons

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
                    matches, comparisons = AhoCorasick(patterns_list).search(dna_sequence)
                    elapsed=time.time()-start
                    
                    total_matches = sum(len(pos_list) for pos_list in matches.values())
                    
                    # Single entry for Aho-Corasick (finds all patterns simultaneously)
                    results.append({
                        "Algorithm": algo,
                        "Pattern": "ALL PATTERNS",
                        "Matches": total_matches,
                        "Comparisons": comparisons,
                        "Time (s)": round(elapsed, 5)
                    })
                    # Add empty row for spacing
                    results.append({
                        "Algorithm": "",
                        "Pattern": "",
                        "Matches": "",
                        "Comparisons": "",
                        "Time (s)": ""
                    })
                else:
                    # Run single-pattern algorithms sequentially for multiple patterns
                    pattern_results = []
                    total_comparisons = 0
                    
                    for pat in patterns_list:
                        start_pat = time.time()
                        result = {"Naïve Search": naive_search,
                                 "KMP": kmp_search,
                                 "Boyer–Moore": boyer_moore_search,
                                 "Rabin–Karp": rabin_karp}[algo](dna_sequence, pat)
                        elapsed_pat = time.time() - start_pat
                        
                        matches_list, comp = result
                        match_count = len(matches_list)
                        total_comparisons += comp
                        
                        results.append({
                            "Algorithm": algo,
                            "Pattern": pat,
                            "Matches": match_count,
                            "Comparisons": comp,
                            "Time (s)": round(elapsed_pat, 5)
                        })
                        pattern_results.append(match_count)
                    
                    # Add TOTAL row for sequential algorithms
                    total_time = sum([r["Time (s)"] for r in results if r["Algorithm"]==algo and r["Pattern"]!="TOTAL" and r["Time (s)"]!=""])
                    results.append({
                        "Algorithm": algo,
                        "Pattern": "TOTAL",
                        "Matches": sum(pattern_results),
                        "Comparisons": total_comparisons,
                        "Time (s)": round(total_time, 5)
                    })
                    # Add empty row for spacing
                    results.append({
                        "Algorithm": "",
                        "Pattern": "",
                        "Matches": "",
                        "Comparisons": "",
                        "Time (s)": ""
                    })
            
            df=pd.DataFrame(results)
            all_results.append(df)
            
            # Display results table
            st.markdown("### 📊 Performance Results")
            st.dataframe(df, use_container_width=True)
            
            # Visualization of pattern matches - separate for each pattern
            st.markdown("### 🎨 Pattern Visualizations (500 bp window for each pattern)")
            
            pattern_colors = ["#FF6B6B", "#4ECDC4", "#45B7D1", "#FFA07A", "#98D8C8", "#F7DC6F", "#BB8FCE", "#85C1E2"]
            
            for idx, pat in enumerate(patterns_list):
                color = pattern_colors[idx % len(pattern_colors)]
                
                # Find all positions of this pattern in entire sequence
                positions = []
                for i in range(len(dna_sequence) - len(pat) + 1):
                    if dna_sequence[i:i+len(pat)] == pat:
                        positions.append(i)
                
                # Always show pattern header
                st.markdown(f"#### 🔍 Pattern: **{pat}**")
                
                if positions:
                    # Get first occurrence position
                    first_match_pos = positions[0]
                    
                    # Extract 500 bp window starting from first match
                    start_pos = first_match_pos
                    end_pos = min(start_pos + 500, len(dna_sequence))
                    preview_seq = dna_sequence[start_pos:end_pos]
                    
                    st.success(f"✅ {len(positions)} match(es) found in entire sequence")
                    st.info(f"📍 Showing sequence from position **{start_pos}** to **{end_pos}** (length: {len(preview_seq)} bp)")
                    
                    # Find pattern positions within this window
                    positions_in_window = [pos - start_pos for pos in positions if start_pos <= pos < end_pos]
                    
                    # Create highlighted HTML
                    highlighted = [False] * len(preview_seq)
                    
                    # Mark all positions that should be highlighted
                    for pos in positions_in_window:
                        for i in range(pos, min(pos + len(pat), len(preview_seq))):
                            if i < len(highlighted):
                                highlighted[i] = True
                    
                    # Build HTML with highlighted sections
                    html_parts = ['<div style="font-family: monospace; line-height: 2; word-wrap: break-word; background-color: #1E2636; padding: 15px; border-radius: 10px; border: 1px solid #00B4D8;">']
                    
                    i = 0
                    while i < len(preview_seq):
                        if highlighted[i]:
                            start_hl = i
                            while i < len(preview_seq) and highlighted[i]:
                                i += 1
                            html_parts.append(f'<span style="background-color: {color}; padding: 2px 4px; border-radius: 3px; font-weight: bold;">{preview_seq[start_hl:i]}</span>')
                        else:
                            html_parts.append(preview_seq[i])
                            i += 1
                    
                    html_parts.append('</div>')
                    
                    # Legend for this pattern
                    legend_html = '<div style="margin-top: 10px; padding: 10px; background-color: #1E2636; border-radius: 8px;">'
                    legend_html += f'<strong>Pattern:</strong> <span style="background-color: {color}; padding: 3px 8px; border-radius: 3px; margin: 0 5px; font-weight: bold;">{pat}</span>'
                    legend_html += f' <span style="margin-left: 15px;">{len(positions_in_window)} occurrence(s) highlighted in this window</span>'
                    legend_html += '</div>'
                    
                    st.markdown(''.join(html_parts) + legend_html, unsafe_allow_html=True)
                    
                else:
                    st.warning(f"⚠️ No matches found for pattern '{pat}' in the sequence.")
                
                st.markdown("---")
            
            # Chart - Performance comparison
            st.markdown("### 📈 Performance Chart")
            # Create chart data using TOTAL rows only
            df_chart = df[(df["Pattern"] == "TOTAL") | (df["Pattern"] == "ALL PATTERNS")].copy()
            
            fig, ax = plt.subplots(figsize=(10,5))
            # Filter out empty rows for chart
            df_chart_filtered = df_chart[df_chart["Algorithm"] != ""]
            ax.bar(df_chart_filtered["Algorithm"], df_chart_filtered["Time (s)"], color="#00B4D8")
            ax.set_ylabel("Time (s)")
            ax.set_title("Algorithm Performance (Total Time)")
            ax.set_xticklabels(df_chart_filtered["Algorithm"], rotation=45, ha='right', fontsize=10)
            plt.tight_layout()
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
