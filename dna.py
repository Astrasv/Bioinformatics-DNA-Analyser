import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import streamlit as st
import altair as alt
from PIL import Image
import itertools
def display_about():
    
    
    image = Image.open('icon.png')
    st.image(image, width=400)

    st.title("About DNA Nucleotide Analysis Web App")  

    with st.container(): 
        st.info(
            "This web app empowers you to analyze DNA sequences. Perform operations like counting nucleotides, calculating GC content, transcribing DNA to RNA, translating RNA to protein, and analyzing amino acid frequencies. It also provides interactive visualizations for enhanced understanding."
        )

        st.subheader("Features:") 
        st.write(
            """
            - **DNA Analysis**: Analyze DNA sequences, count nucleotides, and calculate GC content.
            - **Operations**: Complement DNA, reverse sequences, and get reverse complements.
            - **Pair Frequency**: View nucleotide pair frequency using interactive heatmaps and tabular representations.
            - **RNA Transcription**: Transcribe DNA sequences into RNA.
            - **Protein Analysis**: Translate RNA sequences into proteins and analyze amino acid frequencies.
            """
        )

    with st.container(): 
        st.subheader("About the Developer:")
        st.write(
            """
            This web app is developed by Sudharsan Vanamali. It's built using Streamlit, a powerful tool for creating data apps with Python.
            If you have any feedback or suggestions, feel free to reach out!
            """
        )
        st.write("**Sudharsan Vanamali**")
        
    with st.container(): 
        st.subheader("Technologies Used:")
        st.write(
            """
            - **Python**: The backend of this web app is written in Python.
            - **Streamlit**: Streamlit facilitates the interactive user interface.
            - **Pandas**: Pandas aids in data manipulation and analysis.
            - **Matplotlib and Seaborn**: These libraries enhance data visualization.
            - **Altair**: Altair creates interactive visualizations.
            - **PIL (Python Imaging Library)**: PIL assists with image handling.
            """
        )

    with st.container():  
        st.subheader("How to Use:")
        st.write(
            """
            1. **Enter DNA Sequence**: Input your DNA sequence in the provided text area.
            2. **Explore Analysis Options**: Use the navigation bar to access various analysis options like DNA Analysis, Operations, Pair Frequency, etc.
            3. **Visualize Results**: View interactive visualizations like pie charts and heatmaps for comprehensive analysis.
            """
        )

    with st.container():  
        st.subheader("Additional Resources:")
        st.write(
            """
            - Streamlit Documentation: https://docs.streamlit.io/
            - Pandas Documentation: https://pandas.pydata.org/docs/
            - Matplotlib Documentation: https://matplotlib.org/stable/contents.html
            - Seaborn Documentation: https://seaborn.pydata.org/
            - Altair Documentation: https://altair-viz.github.io/
            """
        )

    with st.container(): 
        st.subheader("Get in Touch:")
        st.write(
            """
            If you have any questions, suggestions, or feedback, feel free to reach out:
            - Email: astrasv247@gmail.com
            - GitHub: [Your GitHub Profile](https://github.com/Astrasv)
            """
        )

def display_main():


    image = Image.open('icon.png')
    st.image(image, width=400)

    st.write("""# DNA Nucleotide Analysis Web App""")


    st.header('Enter DNA sequence')

    sequence_input = ">DNA Query \n Enter Nucleotide and press Ctrl + Enter"

    st.subheader("Sequence Input")
    sequence = st.text_area("", sequence_input, height=170)
    sequence = sequence.splitlines()
    sequence = ''.join(sequence) # Concatenates list to string
    sequence = sequence.replace(" ","")
    
    flag = 0
    for i in sequence.upper():
        if i not in "AGTC":
            st.error("Invalid nucleotide")
            flag = 1
            break
        else:
            flag = 0
    
    if flag == 0:
        
        st.header('Input DNA ')
        st.success(sequence)
        

        st.header('Input DNA Nucleotide Count')

        df = DNA_nucleotide_count(sequence)

        col1, col2, col3, col4 = st.columns(4)
        A_count = int(df[df["nucleotide"] == 'A']["count"])
        G_count = int(df[df["nucleotide"] == 'G']["count"])
        T_count = int(df[df["nucleotide"] == 'T']["count"])
        C_count = int(df[df["nucleotide"] == 'C']["count"])

        col1.metric("A nucleotide", A_count, "0%")
        col2.metric("G nucleotide", G_count, "0%")
        col3.metric("T nucleotide", T_count, "0%")
        col4.metric("C nucleotide", C_count, "0%")

        nucleotides = ['A Nucleotides', 'G Nucleotides', 'T Nucleotides' , 'C Nucleotides']
        counts = [A_count, G_count, T_count, C_count]


        fig, ax = plt.subplots()
        ax.pie(counts, labels=nucleotides, autopct='%1.1f%%', startangle=0)
        ax.axis('equal')  
        

        st.pyplot(fig)
        
        gc_content = calculate_GC_content(sequence)     
        display_GC_content(gc_content)
    
        st.header('Operations')

        operation = st.selectbox("Select Operation", ["Complement", "Reverse", "Reverse-Complement"])
        if operation == "Complement":
            complement_sequence = complement(sequence)
            st.header('Complement DNA')
            st.error(complement_sequence)

            col1, col2, col3, col4 = st.columns(4)

            col1.metric("A nucleotide", T_count, str((T_count - A_count)/100)+'%')
            col2.metric("G nucleotide", C_count, str((C_count - G_count)/100)+'%')
            col3.metric("T nucleotide", A_count, str((A_count - T_count)/100)+'%')
            col4.metric("C nucleotide", G_count, str((G_count - C_count)/100)+'%')
            
        elif operation == "Reverse":
            reverse_sequence = reverse_sequence_func(sequence)
            st.header('Reverse DNA')
            st.error(reverse_sequence)
            col1, col2, col3, col4 = st.columns(4)

            col1.metric("A nucleotide", T_count, "0%")
            col2.metric("G nucleotide", C_count, "0%")
            col3.metric("T nucleotide", A_count, "0%")
            col4.metric("C nucleotide", G_count, "0%")
        elif operation == "Reverse-Complement":
            rev_comp_sequence = reverse_complement_func(sequence)
            st.header('Reverse Complement DNA')
            st.error(rev_comp_sequence)
            col1, col2, col3, col4 = st.columns(4)

            col1.metric("A nucleotide", T_count, str((T_count - A_count)/100)+'%')
            col2.metric("G nucleotide", C_count, str((C_count - G_count)/100)+'%')
            col3.metric("T nucleotide", A_count, str((A_count - T_count)/100)+'%')
            col4.metric("C nucleotide", G_count, str((G_count - C_count)/100)+'%')



        pairs_frequency = nucleotide_pair_frequency(sequence)
        pairs_df = pd.DataFrame(pairs_frequency.items(), columns=['pair', 'frequency'])
        display_nucleotide_pair_frequency(pairs_df)

        rna_sequence = transcribe_DNA_to_RNA(sequence)
        display_RNA_transcription(rna_sequence)

        protein_sequence = translate_RNA_to_protein(rna_sequence)
        display_protein_names(protein_sequence)
        display_translated_protein_sequence(protein_sequence)
        
        
        

        amino_acids_freq = amino_acid_frequency(protein_sequence)
        amino_acids_df = pd.DataFrame(amino_acids_freq.items(), columns=['amino_acid', 'frequency'])
        display_amino_acid_frequency(amino_acids_df)
        
        st.subheader("Percentage composition of Amino Acids")
        amino_acids = list(amino_acids_freq.keys())
        counts = list(amino_acids_freq.values())


        fig, ax = plt.subplots()
        ax.pie(counts, labels=amino_acids, autopct='%1.1f%%', startangle=0 ,pctdistance=1.2, labeldistance=0.9)
        ax.axis('equal')  
        st.pyplot(fig)
        

def DNA_nucleotide_count(seq):
    d = dict([
                ('A',seq.count('A')),
                ('T',seq.count('T')),
                ('G',seq.count('G')),
                ('C',seq.count('C'))
                ])
    df = pd.DataFrame.from_dict(d, orient='index').astype(float)
    df = df.rename({0: 'count'}, axis='columns')
    df.reset_index(inplace=True)
    df = df.rename(columns = {'index':'nucleotide'})
    return df

def complement(seq):
    complement_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    complement_sequence = ''.join(complement_dict[base] for base in seq)
    
    return complement_sequence

def reverse_sequence_func(seq):
    return seq[::-1]

def reverse_complement_func(seq):
    complement_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    complement_sequence = ''.join(complement_dict[base] for base in reversed(seq))
    return complement_sequence

def calculate_GC_content(seq):
    gc_count = seq.count('G') + seq.count('C')
    total_bases = len(seq)
    gc_content = gc_count / total_bases * 100
    return gc_content

def display_GC_content(gc_content):
    if gc_content > 60:
        st.warning(f"🔥 High-GC content {round(gc_content,2)} - Higher Melting termperature")
    elif gc_content < 40:
        st.error(f"🧊 Low-GC content {round(gc_content,2)} - Low Melting termperature")
        
    else:
        st.info(f"🌡️ Neutral -GC content {round(gc_content,2)} - Normal Melting termperature")

def nucleotide_pair_frequency(seq):
    pairs = {}
    nucleotides = ['A', 'T', 'G', 'C']

    # Generate all possible pairs of nucleotides
    all_pairs = [x + y for x in nucleotides for y in nucleotides]

    for pair in all_pairs:
        pairs[pair] = seq.count(pair)

    return pairs

def display_nucleotide_pair_frequency(pairs_df):
    st.header('Nucleotide Pair Frequency')
    
    new_pair_df = {
        'A':{
            'A': pairs_df[pairs_df['pair'] == "AA"]["frequency"].values[0],
            'T': pairs_df[pairs_df['pair'] == "AT"]["frequency"].values[0],
            'G': pairs_df[pairs_df['pair'] == "AG"]["frequency"].values[0],
            'C': pairs_df[pairs_df['pair'] == "AC"]["frequency"].values[0],
        },
        'T':{
            'A': pairs_df[pairs_df['pair'] == "TA"]["frequency"].values[0],
            'T': pairs_df[pairs_df['pair'] == "TT"]["frequency"].values[0],
            'G': pairs_df[pairs_df['pair'] == "TG"]["frequency"].values[0],
            'C': pairs_df[pairs_df['pair'] == "TC"]["frequency"].values[0],
        },
        'G':{
            'A': pairs_df[pairs_df['pair'] == "GA"]["frequency"].values[0],
            'T': pairs_df[pairs_df['pair'] == "GT"]["frequency"].values[0],
            'G': pairs_df[pairs_df['pair'] == "GG"]["frequency"].values[0],
            'C': pairs_df[pairs_df['pair'] == "GC"]["frequency"].values[0],
        },
        'C':{
            'A':pairs_df[pairs_df['pair'] == "CA"]["frequency"].values[0],
            'T':pairs_df[pairs_df['pair'] == "CT"]["frequency"].values[0],
            'G':pairs_df[pairs_df['pair'] == "CG"]["frequency"].values[0],
            'C':pairs_df[pairs_df['pair'] == "CC"]["frequency"].values[0],
        },
    }
    
    new_pair_df = pd.DataFrame(new_pair_df)
    
    st.table(new_pair_df)
    
    fig = plt.figure(figsize=(8, 6))
    sns.heatmap(new_pair_df, annot=True, cmap='coolwarm', fmt='g')
    st.pyplot(fig)

def transcribe_DNA_to_RNA(seq):
    rna_seq = seq.replace('T', 'U')
    return rna_seq

def display_RNA_transcription(rna_sequence):
    st.header('Transcription to RNA')
    st.warning(rna_sequence)
    st.info("Converted Thymine (T) to Uracil (U)")
    
def get_codon():
    codon_table = {
        'UUU': ('F', 'Phenylalanine'), 'UUC': ('F', 'Phenylalanine'), 'UUA': ('L', 'Leucine'), 'UUG': ('L', 'Leucine'),
        'CUU': ('L', 'Leucine'), 'CUC': ('L', 'Leucine'), 'CUA': ('L', 'Leucine'), 'CUG': ('L', 'Leucine'),
        'AUU': ('I', 'Isoleucine'), 'AUC': ('I', 'Isoleucine'), 'AUA': ('I', 'Isoleucine'), 'AUG': ('M', 'Methionine'),
        'GUU': ('V', 'Valine'), 'GUC': ('V', 'Valine'), 'GUA': ('V', 'Valine'), 'GUG': ('V', 'Valine'),
        'UCU': ('S', 'Serine'), 'UCC': ('S', 'Serine'), 'UCA': ('S', 'Serine'), 'UCG': ('S', 'Serine'),
        'CCU': ('P', 'Proline'), 'CCC': ('P', 'Proline'), 'CCA': ('P', 'Proline'), 'CCG': ('P', 'Proline'),
        'ACU': ('T', 'Threonine'), 'ACC': ('T', 'Threonine'), 'ACA': ('T', 'Threonine'), 'ACG': ('T', 'Threonine'),
        'GCU': ('A', 'Alanine'), 'GCC': ('A', 'Alanine'), 'GCA': ('A', 'Alanine'), 'GCG': ('A', 'Alanine'),
        'UAU': ('Y', 'Tyrosine'), 'UAC': ('Y', 'Tyrosine'), 'UAA': ('*', 'Stop'), 'UAG': ('*', 'Stop'),
        'CAU': ('H', 'Histidine'), 'CAC': ('H', 'Histidine'), 'CAA': ('Q', 'Glutamine'), 'CAG': ('Q', 'Glutamine'),
        'AAU': ('N', 'Asparagine'), 'AAC': ('N', 'Asparagine'), 'AAA': ('K', 'Lysine'), 'AAG': ('K', 'Lysine'),
        'GAU': ('D', 'Aspartic Acid'), 'GAC': ('D', 'Aspartic Acid'), 'GAA': ('E', 'Glutamic Acid'), 'GAG': ('E', 'Glutamic Acid'),
        'UGU': ('C', 'Cysteine'), 'UGC': ('C', 'Cysteine'), 'UGA': ('*', 'Stop'), 'UGG': ('W', 'Tryptophan'),
        'CGU': ('R', 'Arginine'), 'CGC': ('R', 'Arginine'), 'CGA': ('R', 'Arginine'), 'CGG': ('R', 'Arginine'),
        'AGU': ('S', 'Serine'), 'AGC': ('S', 'Serine'), 'AGA': ('R', 'Arginine'), 'AGG': ('R', 'Arginine'),
        'GGU': ('G', 'Glycine'), 'GGC': ('G', 'Glycine'), 'GGA': ('G', 'Glycine'), 'GGG': ('G', 'Glycine')
    }
    return codon_table

def display_protein_names(codon):
    protein_names = {
        "Code" : [],
        "Protein Name": []
    }
    codon_table = get_codon()
    for codon, (amino_acid, protein_name) in codon_table.items():
        if amino_acid not in protein_names['Code'] and amino_acid != '*'  :
            protein_names["Code"].append(amino_acid)
            protein_names["Protein Name"].append(protein_name)
        
    
    print(protein_names)

    st.subheader('Protein Names in Codon Table')
    protein_names_df = pd.DataFrame(protein_names)
    protein_names_df.index += 1
    st.table(protein_names_df)
    
def translate_codon(codon):
    codon_table = get_codon()
    return codon_table.get(codon, ('X', 'Unknown'))


def translate_RNA_to_protein(rna_seq):
    
    protein_seq = ""
    for i in range(0, len(rna_seq), 3):
        codon = rna_seq[i:i+3]
        amino_acid = translate_codon(codon)[0]
        if amino_acid == '*':
            break
        protein_seq += amino_acid
    
    return protein_seq

def display_translated_protein_sequence(protein_sequence):
    st.subheader('Translation to Protein and Amino Acid Analysis')
    st.warning(protein_sequence)
    

def amino_acid_frequency(protein_seq):
    amino_acids = {}
    for amino_acid in protein_seq:
        amino_acids[amino_acid] = amino_acids.get(amino_acid, 0) + 1
    return amino_acids

def display_amino_acid_frequency(amino_acids_df):
    st.subheader('Amino Acid Frequency Analysis')
    amino_acids_chart = alt.Chart(amino_acids_df).mark_bar().encode(
        x='amino_acid',
        y='frequency'
    ).properties(
        width=alt.Step(40)
    )
    st.altair_chart(amino_acids_chart)

# Add a navigation bar
nav_choice = st.sidebar.radio("Navigation", ["Analysis", "About"])

# Display different pages based on navigation choice
if nav_choice == "Analysis":
    display_main()
elif nav_choice == "About":
    display_about()
