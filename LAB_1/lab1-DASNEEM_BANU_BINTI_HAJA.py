import streamlit as st
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

Entrez.email = "dasneembanu05@gmail.com"
handle = Entrez.efetch(db='protein', id='P04637', rettype='fasta', retmode='text')

# you may put all your functions at HERE
#Function 2A - retrieve the protein sequence data based on the given id
def retrieve_data(id):
    handle = Entrez.efetch(db='protein', id=id, rettype='fasta', retmode='text')
    record = SeqIO.read(handle, 'fasta')
    handle.close()
    return record

#Function 2B - analyse the given sequence and return all the basic analysis outcome (sequence length, amino acid composition, molecular weight and isoelectric point)
def get_basic_analysis(sequence):
    #calc length using len()
    seq_length = len(sequence)

    seq_analysis = ProteinAnalysis(str(sequence))
    aa_composition = seq_analysis.count_amino_acids()
    molecular_weight = seq_analysis.molecular_weight()
    isoelectric_point = seq_analysis.isoelectric_point()

    return seq_length, aa_composition, molecular_weight, isoelectric_point


st.title('Lab 1 - DASNEEM BANU BINTI HAJA')

protein_id = st.text_input('Enter Uniprot ID')
retrieve = st.button('Retrieve')

if retrieve:
    if protein_id!="":
    	# your function calls and execution of processes at HERE
        record = retrieve_data(protein_id)

        #display basic info
        st.subheader("Protein Basic Information")
        st.write(f"Protein ID: \n{record.id}")
        st.write(f"Name: \n{record.name}")
        st.write(f"Description: \n{record.description}")
        st.write(f"Sequence: \n{record.seq}")

        seq_length, aa_composition, molecular_weight, isoelectric_point = get_basic_analysis(record.seq)

        #display analysis results
        st.subheader("Protein Basic Analysis")
        st.write(f"Sequence length: {seq_length}")
        st.write(f"Molecular weight: {molecular_weight}")
        st.write(f"Isoelectric point: {isoelectric_point}")
        st.write(f"Amino acid composition: {aa_composition}")
    	
    else:
        st.warning('Please enter Uniprot ID')