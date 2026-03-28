# %%
from Bio import SeqIO # BioPython for sequence input/output
from Bio.Seq import Seq # BioPython for sequence object
from Bio.Restriction import RestrictionBatch, Analysis, CommOnly, EcoRI, BamHI, HindIII
# CommOnly = BioPython built-in curated list of 200 commercially available enzymes (the one you buy and use in lab)

import matplotlib.pyplot as plt
import matplotlib_inline
import matplotlib.patches as mpatches
import seaborn as sns

import pandas as pd
import numpy as np
import streamlit as st

# %% [markdown]
# # Restriction enzyme sites/specific seq motifs.

# %% [markdown]
# Restriction enzymes (also called restriction endonucleases) are proteins that cut DNA at specific short sequences called recognition sites. 

# %% [markdown]
# Objective: read FASTA file, print the sequence. Calculate GC content. Compare GC content across different genes/organisms (identify pattern). Translate DNA seq to protein. Search for restriction enzyme sites/specific seq motifs. Create visualizations showing where interesting features occur along the sequence.

# %% [markdown]
# This is fundamental in molecular biology for cloning, gel electrophoresis, genome mapping, and diagnostics.

# %% [markdown]
# FASTA Parser → Find Cut Sites → Wrap & Table → Visualize Map → Density Plot → Streamlit UI

# %%
# 1 FASTA file parser

def read_fasta(file):
    sequences = {}
    if isinstance(file, str): #isinstance() function returns True if the specified object is of the specified type, otherwise False.
        # so if the file is just string (like file path, so open with r (read) mode)
        f = open(file, "r")
    else:
        f = file  # file-like object from uploader, just use it directly

    name = ""
    seq = ""
    for line in f:
        line = line.decode("utf-8").strip() if isinstance(line, bytes) else line.strip()
        # decode("utf-8") -> character encode if it is bytes
        if line.startswith(">"):
            if name:
                sequences[name] = seq
            name = line[1:] # exclude the >, start from idx 1
            seq = ""
        else:
            seq += line # combine all line
    if name:
        sequences[name] = seq

    if isinstance(file, str):
        f.close()
    return sequences

# 2. To find restriction enzyme cut site

def find_restriction_enzymes(seq: str, enzymes = None):
    bio_seq = Seq(seq.upper()) # convert string into BioPython object and make it UPPERCASE
    rb = RestrictionBatch(CommOnly) if enzymes is None else RestrictionBatch(enzymes)
    # RB(RestrictionBatch represents collection of restriction enzymes (commOnly = commercially available restriction enzyme usually used in lab))
    ana = Analysis(rb, bio_seq, linear=True) # to treat the sequence as linear
    results = ana.full() #Returns all enzymes and their cut positions

    filtered_results = {}

    # Loop through all enzymes and positions
    for enzymes, positions in results.items():
        if positions:
            filtered_results[str(enzymes)] = positions
    return filtered_results

# 3. Wrap repetitive task and to put it into Pandas Table

def analyze_sequence(seq, enzymes=None):
    """Runs all sequence-analysis steps automatically."""
    seq = seq.replace("\n", "").strip()
    results = find_restriction_enzymes(seq, enzymes)
    df = restriction_table({"Sequence": seq}, enzymes)
    return seq, results, df

def restriction_table(data: dict, enzymes= None):
    
    record = {}

    for name, seq in data.items():
        site = find_restriction_enzymes(seq, enzymes)
        record[name] = {enz: len(pos) for enz, pos in site.items()}
    

    df = pd.DataFrame(record).T.fillna(0).astype(int)
    df = df.loc[:, (df>0).any(axis = 0)]
    return df

# 4. Visualize the restriction map

def plot_restriction_map(seq_len, enzyme_dict, seq_name="Sequence"):
    # enzyme_dict = keys = enzyme names, values = list of cut positions

    fig, ax = plt.subplots(figsize=(20, 2.5)) #size of the figure

    ax.hlines(0.5, 0, seq_len, colors='black', linewidth=2)
    
    y_top, y_bottom = 0.9, 0.1
    toggle = True #alternate the placement of enzyme labels so they don’t overlap

    for enz_name, positions in enzyme_dict.items():
        for pos in positions:
            ax.vlines(pos, 0.15, 0.85, color='red', linewidth=1)
            y_label = y_top if toggle else y_bottom
            toggle = not toggle
            # the enzyme names zig-zagging above and below the DNA line, making the plot readable & not overlap
            ax.text(pos, 0.9, enz_name, rotation=90, ha='center', va='bottom', fontsize=8)

    ax.set_xlim(0, seq_len)
    ax.set_ylim(0, 1)
    ax.set_xlabel("Position (bp)")
    ax.set_yticks([])
    ax.set_title(f"Restriction map of {seq_name}")
    plt.tight_layout()
    st.pyplot(fig)


# 5. Plot feature density
def plot_feature_density(enzyme_dict):
    all_positions = []
    for pos_list in enzyme_dict.values():
        all_positions.extend(pos_list)
    if len(all_positions) == 0:
        st.warning("No restriction sites found.")
        return

    fig, ax = plt.subplots(figsize=(10, 3))
    sns.kdeplot(all_positions, fill=True, ax=ax)
    ax.set_xlabel("Position (bp)")
    ax.set_ylabel("Density")
    ax.set_title("Restriction Site Density")
    st.pyplot(fig)


# 6. Streamlit
st.title("Sequence Analyzer Web App")
st.write("Paste your DNA sequence to get automated analysis.")

sequence_input = st.text_area("Enter DNA Sequence", height=200)

uploaded_file = st.file_uploader("Or upload a FASTA file", type=["fasta","fa","fna"])

# --- Widgets FIRST, before any results ---
enzyme_map = {"EcoRI": EcoRI, "BamHI": BamHI, "HindIII": HindIII}
enzyme_choice = st.multiselect("Choose enzymes (leave empty for all)", list(enzyme_map.keys()))
selected = [enzyme_map[e] for e in enzyme_choice] or None

if st.button("Analyze"):
    if uploaded_file:
        sequences = read_fasta(uploaded_file)
    elif sequence_input.strip():
        sequences = {"Sequence": sequence_input}
    else:
        st.error("Please enter a sequence or upload a file.")
        st.stop()

    for name, seq in sequences.items():
        st.subheader(f"Results for: {name}")
        seq, results, df = analyze_sequence(seq, enzymes=selected)  # pass selected here

        st.subheader("Restriction Enzyme Summary Table")
        st.dataframe(df)
        csv = df.to_csv().encode("utf-8")
        st.download_button("Download CSV", csv, f"{name}_results.csv", "text/csv")

        st.subheader("Restriction Map")
        plot_restriction_map(len(seq), results, seq_name=name)

        st.subheader("Restriction Site Density")
        plot_feature_density(results)

    st.success("Analysis complete!")


