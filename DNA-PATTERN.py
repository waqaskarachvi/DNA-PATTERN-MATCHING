# ================================
# DNA Pattern Matching Analyzer🧬 
# ================================
# This application compares different string matching algorithms for DNA sequence analysis
# It allows users to upload FASTA files or input DNA sequences manually and search for patterns

import streamlit as st
import re, time, io
import pandas as pd
import matplotlib.pyplot as plt
from collections import deque, defaultdict

# ==========================
# PAGE CONFIG
# ==========================
# Configure the Streamlit page with title, icon, and wide layout for better visualization
st.set_page_config(page_title="DNA Pattern Matching Analyzer", page_icon="🧬", layout="wide")

# ==========================
# STYLING
# ==========================
# Custom CSS styling for dark theme with blue accent colors
# Makes the app visually appealing with proper color schemes for DNA analysis
st.markdown("""
    <style>
        /* Dark background for main app body */
        body, .stApp { background-color: #0E1117; color: #FAFAFA; }
        
        /* Blue headers for visual hierarchy */
        h1, h2, h3 { color: #00B4D8 !important; }
        
        /* Styled buttons with hover effects */
        .stButton>button {
            background-color: #00B4D8; color: white; font-weight: bold;
            border-radius: 8px; border: none; padding: 0.6em 1.2em;
        }
        .stButton>button:hover { background-color: #0077B6; }
        
        /* Result boxes with monospace font for DNA sequences */
        .result-box { background-color: #1E2636; padding: 15px; border-radius: 10px; border: 1px solid #00B4D8;
                      font-family: monospace; line-height: 1.6; word-wrap: break-word; }
        
        /* Highlight color for matched patterns */
        .highlight { color: #FFD60A; font-weight: bold; }
    </style>
""", unsafe_allow_html=True)

# ==========================
# HEADER
# ==========================
# Display main title and subtitle centered on the page
st.markdown("<h1 style='text-align:center;'>🧬 DNA Pattern Matching Analyzer</h1>", unsafe_allow_html=True)
st.markdown("<p style='text-align:center;'>Compare string matching algorithms for DNA sequences</p>", unsafe_allow_html=True)
st.markdown("---", unsafe_allow_html=True)

# ==========================
# FILE UPLOAD / INPUT
# ==========================
# Allow users to upload multiple FASTA files containing DNA sequences
uploaded_files = st.file_uploader(
    "📁 Upload FASTA files (you can select multiple)",
    type=["fasta","fa","txt"],
    accept_multiple_files=True
)

# Dictionary to store all loaded sequences with their headers as keys
sequences = {}

# Process uploaded FASTA files
if uploaded_files:
    for uploaded_file in uploaded_files:
        # Decode file content from bytes to string
        fasta = uploaded_file.getvalue().decode("utf-8")
        lines = fasta.strip().split("\n")
        
        # Extract header (first line starting with >) or use filename as header
        header = lines[0] if lines[0].startswith(">") else uploaded_file.name
        
        # Join all non-header lines into single sequence and convert to uppercase
        sequence = "".join([l.strip() for l in lines if not l.startswith(">")]).upper()
        
        # Remove any non-DNA characters (keep only A, T, C, G)
        sequence = re.sub(r'[^ATCG]', '', sequence)
        
        # Store sequence in dictionary and show success message
        sequences[header] = sequence
        st.success(f"✅ Loaded: {header} ({len(sequence)} bp)")
else:
    # If no files uploaded, provide text area for manual DNA sequence entry
    seq_input = st.text_area("🧬 Enter DNA Sequence", placeholder="ATCGGATCGATCG...", height=120).strip().upper()
    if seq_input:
        # Clean manual input to contain only valid DNA bases
        sequences["Manual Entry"] = re.sub(r'[^ATCG]', '', seq_input)

# Input field for patterns to search (supports multiple comma-separated patterns)
pattern_input = st.text_input(
    "🔍 Enter Pattern(s) to Search (comma separated for multiple)",
    placeholder="CGATCGA,ATGCGT"
).strip().upper()

# List of available algorithms for pattern matching
algorithms = ["Naïve Search", "KMP", "Boyer–Moore", "Rabin–Karp", "Aho–Corasick"]

# Multi-select widget allowing users to choose which algorithms to run
selected_algos = st.multiselect("⚙️ Select Algorithms", algorithms, default=algorithms)

# ============
# ALGORITHMS 
# ============

# NAÏVE SEARCH ALGORITHM
# Simplest approach: checks every possible position in the text
# Time Complexity: O(n*m) where n=text length, m=pattern length
def naive_search(text, pattern):
    comparisons = 0  # Counter to track number of character comparisons
    results = []      # List to store positions where pattern is found
    
    # Try every possible starting position in the text
    for i in range(len(text)-len(pattern)+1):
        comparisons += 1
        match = True
        
        # Check if pattern matches at current position
        for j in range(len(pattern)):
            comparisons += 1
            if text[i+j] != pattern[j]:
                match = False
                break
        
        # If all characters matched, record this position
        if match:
            results.append(i)
    
    return results, comparisons

# KMP (KNUTH-MORRIS-PRATT) ALGORITHM
# Uses preprocessing to avoid re-checking characters
# Time Complexity: O(n+m) - much more efficient for longer patterns
def kmp_search(text, pattern):
    comparisons = 0
    
    # Build LPS (Longest Proper Prefix which is also Suffix) array
    # This array helps skip characters we've already matched
    lps=[0]*len(pattern)
    j=0
    
    # Preprocess pattern to build LPS array
    for i in range(1,len(pattern)):
        while j>0 and pattern[i]!=pattern[j]: 
            j=lps[j-1]
            comparisons += 1
        comparisons += 1
        if pattern[i]==pattern[j]: 
            j+=1
            lps[i]=j
    
    # Search for pattern in text using LPS array
    res=[]
    j=0
    for i in range(len(text)):
        # Use LPS to skip already matched characters
        while j>0 and text[i]!=pattern[j]: 
            j=lps[j-1]
            comparisons += 1
        comparisons += 1
        
        # Character matches, move to next pattern character
        if text[i]==pattern[j]: 
            j+=1
        
        # Full pattern found
        if j==len(pattern): 
            res.append(i-j+1)
            j=lps[j-1]  # Continue searching for more matches
    
    return res, comparisons

# BOYER-MOORE ALGORITHM
# Searches from right to left, uses bad character heuristic
# Time Complexity: O(n/m) best case, O(n*m) worst case
def boyer_moore_search(text, pattern):
    comparisons = 0
    m,n=len(pattern),len(text)
    
    # Build bad character table - stores last occurrence of each character
    bad_char={pattern[i]:i for i in range(m)}
    
    res=[]
    s=0  # Current shift position
    
    # Slide pattern over text
    while s<=n-m:
        j=m-1  # Start comparing from rightmost character of pattern
        
        # Compare pattern from right to left
        while j>=0:
            comparisons += 1
            if pattern[j]!=text[s+j]:
                break
            j-=1
        
        if j<0:  # Pattern found
            res.append(s)
            # Shift pattern to align with next occurrence of current character
            s+=(m-bad_char.get(text[s+m],-1)) if s+m<n else 1
        else:  
            # Shift pattern based on bad character heuristic
            s+=max(1,j-bad_char.get(text[s+j],-1))
    
    return res, comparisons

# RABIN-KARP ALGORITHM
# Uses rolling hash function for fast pattern matching
# Time Complexity: O(n+m) average case, O(n*m) worst case
def rabin_karp(text, pattern, d=256, q=101):
    comparisons = 0
    m,n=len(pattern),len(text)
    p=t=0  # Pattern hash and text window hash
    h=pow(d,m-1)%q  # Hash multiplier for removing leading digit
    res=[]
    
    # Calculate hash for pattern and first text window
    for i in range(m): 
        p=(d*p+ord(pattern[i]))%q
        t=(d*t+ord(text[i]))%q
    
    # Slide pattern over text one character at a time
    for s in range(n-m+1):
        comparisons += 1
        
        # If hash values match, verify actual string match
        if p==t:
            if text[s:s+m]==pattern: 
                res.append(s)
                comparisons += m  # Count character comparisons during verification
        
        # Calculate hash for next window using rolling hash
        if s<n-m: 
            t=(d*(t-ord(text[s])*h)+ord(text[s+m]))%q
            t+=q if t<0 else 0  # Ensure hash is positive
    
    return res, comparisons

# ==============
# AHO-CORASICK 
# ==============
# AHO-CORASICK ALGORITHM
# Efficient for searching multiple patterns simultaneously
# Uses trie structure with failure links
# Time Complexity: O(n+m+z) where z=number of matches

# Node structure for Aho-Corasick trie
class AhoNode:
    def __init__(self): 
        self.children={}     # Dictionary of child nodes
        self.fail=None       # Failure link for mismatches
        self.output=[]       # Patterns that end at this node

# Main Aho-Corasick automaton class
class AhoCorasick:
    def __init__(self,patterns):
        self.root=AhoNode()
        self.build_trie(patterns)        # Build trie from patterns
        self.build_failure_links()       # Add failure links for efficient traversal
    
    # Build trie structure from all patterns
    def build_trie(self,patterns):
        for pat in patterns:
            node=self.root
            # Insert pattern character by character
            for c in pat:
                if c not in node.children: 
                    node.children[c]=AhoNode()
                node=node.children[c]
            # Mark end of pattern
            node.output.append(pat)
    
    # Build failure links for efficient pattern matching
    # Failure links point to longest proper suffix that's also a prefix
    def build_failure_links(self):
        queue=deque()
        
        # Root's children fail back to root
        for child in self.root.children.values(): 
            child.fail=self.root
            queue.append(child)
        
        # BFS to build failure links for all nodes
        while queue:
            current=queue.popleft()
            for c,child in current.children.items():
                f=current.fail
                
                # Find deepest node that has a child with character c
                while f and c not in f.children: 
                    f=f.fail
                
                # Set failure link
                child.fail=f.children[c] if f and c in f.children else self.root
                
                # Inherit output from failure link
                child.output+=child.fail.output
                queue.append(child)
    
    # Search for all patterns in text simultaneously
    def search(self,text):
        node=self.root
        res=defaultdict(list)  # Dictionary mapping patterns to their positions
        comparisons = 0
        
        # Process each character in text
        for i,c in enumerate(text):
            comparisons += 1
            
            # Follow failure links until we find a match or reach root
            while node and c not in node.children: 
                node=node.fail
                comparisons += 1
            
            if not node: 
                node=self.root
                continue
            
            # Move to next node
            node=node.children[c]
            
            # Record all patterns that end at current position
            for pat in node.output: 
                res[pat].append(i-len(pat)+1)
        
        return res, comparisons

# ==========================
# RUN ANALYSIS
# ==========================
# Initialize session state to store results between reruns
if "results_stored" not in st.session_state: 
    st.session_state.results_stored=None

# Main search button - triggers pattern matching analysis
if st.button("🔍 Search Pattern"):
    # Validate that we have both sequences and patterns
    if not sequences or not pattern_input:
        st.warning("⚠️ Please enter sequence(s) and pattern(s).")
    else:
        # Parse comma-separated patterns
        patterns_list = [p.strip() for p in pattern_input.split(",") if p.strip()]
        all_results=[]
        
        # Process each loaded DNA sequence
        for header,dna_sequence in sequences.items():
            st.markdown(f"## 🧫 Results for **{header}** ({len(dna_sequence)} bp)")
            results=[]
            
            # Run each selected algorithm
            for algo in selected_algos:
                # AHO-CORASICK: processes all patterns at once
                if algo=="Aho–Corasick":
                    # Requires at least 2 patterns
                    if len(patterns_list)<2: 
                        st.warning("⚠️ Aho–Corasick requires multiple patterns.")
                        continue
                    
                    # Time the algorithm execution
                    start=time.time()
                    matches, comparisons = AhoCorasick(patterns_list).search(dna_sequence)
                    elapsed=time.time()-start
                    
                    # Count total matches across all patterns
                    total_matches = sum(len(pos_list) for pos_list in matches.values())
                    
                    # Record results for Aho-Corasick (one entry for all patterns)
                    results.append({
                        "Algorithm": algo,
                        "Pattern": "ALL PATTERNS",
                        "Matches": total_matches,
                        "Comparisons": comparisons,
                        "Time (s)": round(elapsed, 5)
                    })
                    # Add empty row for visual spacing
                    results.append({
                        "Algorithm": "",
                        "Pattern": "",
                        "Matches": "",
                        "Comparisons": "",
                        "Time (s)": ""
                    })
                else:
                    # SINGLE-PATTERN ALGORITHMS: run separately for each pattern
                    pattern_results = []
                    total_comparisons = 0
                    
                    # Test each pattern individually
                    for pat in patterns_list:
                        start_pat = time.time()
                        
                        # Call appropriate algorithm function
                        result = {"Naïve Search": naive_search,
                                 "KMP": kmp_search,
                                 "Boyer–Moore": boyer_moore_search,
                                 "Rabin–Karp": rabin_karp}[algo](dna_sequence, pat)
                        
                        elapsed_pat = time.time() - start_pat
                        
                        # Extract matches and comparison count
                        matches_list, comp = result
                        match_count = len(matches_list)
                        total_comparisons += comp
                        
                        # Record results for this pattern
                        results.append({
                            "Algorithm": algo,
                            "Pattern": pat,
                            "Matches": match_count,
                            "Comparisons": comp,
                            "Time (s)": round(elapsed_pat, 5)
                        })
                        pattern_results.append(match_count)
                    
                    # Add TOTAL row summarizing all patterns for this algorithm
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
            
            # Convert results to DataFrame for display
            df=pd.DataFrame(results)
            all_results.append(df)
            
            # ==========================
            # DISPLAY RESULTS TABLE
            # ==========================
            st.markdown("### 📊 Performance Results")
            st.dataframe(df, use_container_width=True)
            
            # ==========================
            # PATTERN VISUALIZATIONS
            # ==========================
            st.markdown("### 🎨 Pattern Visualizations (500 bp window for each pattern)")
            
            # Color palette for highlighting different patterns
            pattern_colors = ["#FF6B6B", "#4ECDC4", "#45B7D1", "#FFA07A", "#98D8C8", "#F7DC6F", "#BB8FCE", "#85C1E2"]
            
            # Visualize each pattern separately
            for idx, pat in enumerate(patterns_list):
                color = pattern_colors[idx % len(pattern_colors)]
                
                # Find all occurrences of this pattern in the entire sequence
                positions = []
                for i in range(len(dna_sequence) - len(pat) + 1):
                    if dna_sequence[i:i+len(pat)] == pat:
                        positions.append(i)
                
                # Display pattern header
                st.markdown(f"#### 🔍 Pattern: **{pat}**")
                
                if positions:
                    # Extract 500 bp window starting from first match
                    first_match_pos = positions[0]
                    start_pos = first_match_pos
                    end_pos = min(start_pos + 500, len(dna_sequence))
                    preview_seq = dna_sequence[start_pos:end_pos]
                    
                    # Display match statistics
                    st.success(f"✅ {len(positions)} match(es) found in entire sequence")
                    st.info(f"📍 Showing sequence from position **{start_pos}** to **{end_pos}** (length: {len(preview_seq)} bp)")
                    
                    # Find pattern positions within the 500 bp window
                    positions_in_window = [pos - start_pos for pos in positions if start_pos <= pos < end_pos]
                    
                    # Create boolean array marking which positions to highlight
                    highlighted = [False] * len(preview_seq)
                    
                    # Mark all positions that should be highlighted
                    for pos in positions_in_window:
                        for i in range(pos, min(pos + len(pat), len(preview_seq))):
                            if i < len(highlighted):
                                highlighted[i] = True
                    
                    # Build HTML with syntax highlighting for matched patterns
                    html_parts = ['<div style="font-family: monospace; line-height: 2; word-wrap: break-word; background-color: #1E2636; padding: 15px; border-radius: 10px; border: 1px solid #00B4D8;">']
                    
                    # Process sequence and create highlighted spans
                    i = 0
                    while i < len(preview_seq):
                        if highlighted[i]:
                            # Found start of highlighted region
                            start_hl = i
                            # Find end of highlighted region
                            while i < len(preview_seq) and highlighted[i]:
                                i += 1
                            # Add highlighted span with pattern-specific color
                            html_parts.append(f'<span style="background-color: {color}; padding: 2px 4px; border-radius: 3px; font-weight: bold;">{preview_seq[start_hl:i]}</span>')
                        else:
                            # Regular unhighlighted character
                            html_parts.append(preview_seq[i])
                            i += 1
                    
                    html_parts.append('</div>')
                    
                    # Create legend showing pattern color and match count
                    legend_html = '<div style="margin-top: 10px; padding: 10px; background-color: #1E2636; border-radius: 8px;">'
                    legend_html += f'<strong>Pattern:</strong> <span style="background-color: {color}; padding: 3px 8px; border-radius: 3px; margin: 0 5px; font-weight: bold;">{pat}</span>'
                    legend_html += f' <span style="margin-left: 15px;">{len(positions_in_window)} occurrence(s) highlighted in this window</span>'
                    legend_html += '</div>'
                    
                    # Display highlighted sequence
                    st.markdown(''.join(html_parts) + legend_html, unsafe_allow_html=True)
                    
                else:
                    # No matches found for this pattern
                    st.warning(f"⚠️ No matches found for pattern '{pat}' in the sequence.")
                
                st.markdown("---")
            
            # ==========================
            # PERFORMANCE CHART
            # ==========================
            st.markdown("### 📈 Performance Chart")
            
            # Create chart showing only total/summary times for each algorithm
            df_chart = df[(df["Pattern"] == "TOTAL") | (df["Pattern"] == "ALL PATTERNS")].copy()
            
            # Create bar chart
            fig, ax = plt.subplots(figsize=(10,5))
            df_chart_filtered = df_chart[df_chart["Algorithm"] != ""]
            ax.bar(df_chart_filtered["Algorithm"], df_chart_filtered["Time (s)"], color="#00B4D8")
            ax.set_ylabel("Time (s)")
            ax.set_title("Algorithm Performance (Total Time)")
            ax.set_xticklabels(df_chart_filtered["Algorithm"], rotation=45, ha='right', fontsize=10)
            plt.tight_layout()
            st.pyplot(fig)

        # Combine all results from multiple sequences
        combined_df=pd.concat(all_results,ignore_index=True)
        
        # Store results in session state for download
        st.session_state.results_stored=combined_df

# ==========================
# DOWNLOAD CSV
# ==========================
# Provide download button if results are available
if st.session_state.results_stored is not None:
    # Create CSV buffer
    csv_buffer=io.StringIO()
    st.session_state.results_stored.to_csv(csv_buffer,index=False)
    
    # Display download button
    st.download_button(
        label="📥 Download Results as CSV",
        data=csv_buffer.getvalue(),
        file_name="dna_pattern_results.csv",
        mime="text/csv"
    )
