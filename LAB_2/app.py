
# BIOINFORMATICS LAB 1 + LAB 2
# Protein Analysis and Pairwise Sequence Alignment

import streamlit as st
import pandas as pd

from Bio import Entrez
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

from alignment import global_alignment, local_alignment


# PAGE CONFIGURATION

st.set_page_config(
    page_title="Bioinformatics Lab 1 and Lab 2",
    page_icon="🧬",
    layout="wide"
)


# LAB 1 - PROTEIN ANALYSIS FUNCTIONS

# Function 2A - Retrieve protein sequence data based on ID
def retrieve_data(protein_id):

    Entrez.email = "dasneembanu05@gmail.com"

    handle = Entrez.efetch(
        db="protein",
        id=protein_id,
        rettype="fasta",
        retmode="text"
    )

    record = SeqIO.read(handle, "fasta")
    handle.close()

    return record


# Function 2B - Perform basic protein analysis
def get_basic_analysis(sequence):

    # Calculate sequence length
    seq_length = len(sequence)

    # Analyse protein sequence
    seq_analysis = ProteinAnalysis(str(sequence))

    # Calculate amino acid composition
    aa_composition = seq_analysis.count_amino_acids()

    # Calculate molecular weight
    molecular_weight = seq_analysis.molecular_weight()

    # Calculate isoelectric point
    isoelectric_point = seq_analysis.isoelectric_point()

    return (
        seq_length,
        aa_composition,
        molecular_weight,
        isoelectric_point
    )


# LAB 2 - SEQUENCE VALIDATION

def clean_sequence(sequence):

    # Remove spaces and convert to uppercase
    sequence = "".join(sequence.split()).upper()

    return sequence


def validate_sequence(sequence):

    if sequence == "":
        return False, "Please enter a sequence."

    if not sequence.isalpha():
        return False, "Sequence must contain letters only."

    if not sequence.isascii():
        return False, "Please use standard English alphabet characters."

    return True, ""


# LAB 2 - DISPLAY ALIGNMENT RESULTS

def display_alignment_result(result, seq1, seq2):

    st.subheader("Alignment Results")

    # Display alignment score

    st.markdown("### Alignment Score")

    st.metric(
        label="Score",
        value=result["score"]
    )


    # Display aligned sequences

    st.markdown("### Aligned Sequences")

    aligned_seq1 = result["aligned_seq1"]
    aligned_seq2 = result["aligned_seq2"]

    st.code(
        f"Sequence 1: {aligned_seq1}\n"
        f"Sequence 2: {aligned_seq2}",
        language=None
    )

    # Display alignment statistics
    matches = 0
    mismatches = 0
    gaps = 0

    for char1, char2 in zip(aligned_seq1, aligned_seq2):

        if char1 == "-" or char2 == "-":
            gaps += 1

        elif char1 == char2:
            matches += 1

        else:
            mismatches += 1

    st.markdown("### Alignment Statistics")

    col1, col2, col3 = st.columns(3)

    col1.metric("Matches", matches)
    col2.metric("Mismatches", mismatches)
    col3.metric("Gaps", gaps)

    # Display scoring matrix
    st.markdown("### Scoring Matrix")

    st.caption(
        "Green cells indicate the traceback path."
    )

    matrix = result["matrix"]

    # Create unique row labels
    row_labels = ["0 (-)"] + [
        f"{i} ({char})"
        for i, char in enumerate(seq1, start=1)
    ]

    # Create unique column labels
    col_labels = ["0 (-)"] + [
        f"{j} ({char})"
        for j, char in enumerate(seq2, start=1)
    ]

    # Create DataFrame with unique labels
    matrix_df = pd.DataFrame(
        matrix,
        index=row_labels,
        columns=col_labels
    )

    # Highlight traceback path
    traceback_path = result["traceback_path"]

    def highlight_path(data):

        # Create an empty style DataFrame
        styles = pd.DataFrame(
            "",
            index=data.index,
            columns=data.columns
        )

        # Highlight traceback cells using row/column positions
        for i, j in traceback_path:

            styles.iloc[i, j] = (
                "background-color: lightgreen; "
                "color: black; "
                "font-weight: bold;"
            )

        return styles

    # Apply styling
    styled_matrix = matrix_df.style.apply(
        highlight_path,
        axis=None
    )

    st.dataframe(
        styled_matrix,
        use_container_width=True
    )

    # Display traceback coordinates
    with st.expander("View Traceback Path"):

        st.write(
            "The traceback path shows the cells visited "
            "when reconstructing the alignment."
        )

        st.write(traceback_path)


# MAIN APPLICATION
st.title("🧬 Bioinformatics Laboratory Application")

st.write(
    "This application combines Lab 1 Protein Analysis "
    "and Lab 2 Pairwise Sequence Alignment."
)

st.divider()


# Store retrieved protein in session state
if "protein_record" not in st.session_state:

    st.session_state.protein_record = None



# APPLICATION TABS
tab1, tab2 = st.tabs([
    "Lab 1 - Protein Analysis",
    "Lab 2 - Pairwise Alignment"
])



# TAB 1 - LAB 1 PROTEIN ANALYSIS
with tab1:

    st.header("Lab 1 - Protein Analysis")

    st.write(
        "Retrieve a protein sequence using its protein ID "
        "and perform basic protein analysis."
    )

    protein_id = st.text_input(
        "Enter Protein ID",
        placeholder="Example: P04637"
    )

    retrieve = st.button(
        "Retrieve Protein",
        key="retrieve_protein"
    )

    if retrieve:

        if protein_id.strip() != "":

            try:

                with st.spinner("Retrieving protein data..."):

                    record = retrieve_data(protein_id.strip())

                    st.session_state.protein_record = record

                st.success("Protein retrieved successfully!")

            except Exception as e:

                st.error(
                    f"Unable to retrieve protein data: {e}"
                )

        else:

            st.warning("Please enter a Protein ID.")

    # Display retrieved protein information
    if st.session_state.protein_record is not None:

        record = st.session_state.protein_record

        st.divider()

        st.subheader("Protein Basic Information")

        st.write(f"**Protein ID:** {record.id}")

        st.write(f"**Name:** {record.name}")

        st.write(f"**Description:** {record.description}")

        st.write("**Sequence:**")

        st.code(str(record.seq), language=None)

        # Perform protein analysis
        try:

            (
                seq_length,
                aa_composition,
                molecular_weight,
                isoelectric_point
            ) = get_basic_analysis(record.seq)

            st.divider()

            st.subheader("Protein Basic Analysis")

            col1, col2 = st.columns(2)

            col1.metric(
                "Sequence Length",
                seq_length
            )

            col2.metric(
                "Molecular Weight",
                f"{molecular_weight:.2f} Da"
            )

            st.metric(
                "Isoelectric Point",
                f"{isoelectric_point:.2f}"
            )

            st.subheader("Amino Acid Composition")

            aa_df = pd.DataFrame(
                list(aa_composition.items()),
                columns=["Amino Acid", "Count"]
            )

            st.dataframe(
                aa_df,
                use_container_width=True
            )

        except Exception as e:

            st.error(
                f"Unable to analyse this protein sequence: {e}"
            )


# TAB 2 - LAB 2 PAIRWISE SEQUENCE ALIGNMENT
with tab2:

    st.header("Lab 2 - Pairwise Sequence Alignment")

    st.write(
        "Compare two biological sequences using "
        "Global Alignment or Local Alignment."
    )

    st.divider()

    # Sequence 1 input
    st.subheader("Sequence 1")

    source1 = st.selectbox(
        "Choose Sequence 1 source",
        [
            "Enter sequence manually",
            "Use retrieved Lab 1 protein"
        ],
        key="source1"
    )

    if source1 == "Enter sequence manually":

        input_seq1 = st.text_area(
            "Enter Sequence 1",
            placeholder="Example: ATCG",
            key="input_seq1"
        )

    else:

        if st.session_state.protein_record is not None:

            input_seq1 = str(
                st.session_state.protein_record.seq
            )

            st.success(
                "Using the retrieved Lab 1 protein sequence."
            )

            st.text_area(
                "Sequence 1",
                value=input_seq1,
                disabled=True,
                key="retrieved_seq1"
            )

        else:

            input_seq1 = ""

            st.warning(
                "Please retrieve a protein in Lab 1 first."
            )

    st.divider()


    # Sequence 2 input
    st.subheader("Sequence 2")

    source2 = st.selectbox(
        "Choose Sequence 2 source",
        [
            "Enter sequence manually",
            "Use retrieved Lab 1 protein"
        ],
        key="source2"
    )

    if source2 == "Enter sequence manually":

        input_seq2 = st.text_area(
            "Enter Sequence 2",
            placeholder="Example: ATGG",
            key="input_seq2"
        )

    else:

        if st.session_state.protein_record is not None:

            input_seq2 = str(
                st.session_state.protein_record.seq
            )

            st.success(
                "Using the retrieved Lab 1 protein sequence."
            )

            st.text_area(
                "Sequence 2",
                value=input_seq2,
                disabled=True,
                key="retrieved_seq2"
            )

        else:

            input_seq2 = ""

            st.warning(
                "Please retrieve a protein in Lab 1 first."
            )

    st.divider()


    # Alignment method
    st.subheader("Alignment Method")

    alignment_method = st.radio(
        "Choose alignment method",
        [
            "Global Alignment (Needleman-Wunsch)",
            "Local Alignment (Smith-Waterman)"
        ]
    )

    st.divider()

    # Scoring parameters
    st.subheader("Scoring Parameters")

    col1, col2, col3 = st.columns(3)

    with col1:

        match_score = st.number_input(
            "Match Score",
            value=2,
            step=1
        )

    with col2:

        mismatch_score = st.number_input(
            "Mismatch Score",
            value=-1,
            step=1
        )

    with col3:

        gap_penalty = st.number_input(
            "Gap Penalty",
            value=-2,
            step=1
        )

    st.divider()

    # Perform alignment
    perform_alignment = st.button(
        "Perform Alignment",
        type="primary",
        use_container_width=True
    )

    if perform_alignment:

        # Clean sequences
        seq1 = clean_sequence(input_seq1)

        seq2 = clean_sequence(input_seq2)

        # Validate both sequences
        valid1, message1 = validate_sequence(seq1)

        valid2, message2 = validate_sequence(seq2)

        if not valid1:

            st.error(f"Sequence 1: {message1}")

        elif not valid2:

            st.error(f"Sequence 2: {message2}")

        else:

            try:

                with st.spinner("Performing sequence alignment..."):

                    # Global Alignment
                    if alignment_method.startswith("Global"):

                        result = global_alignment(
                            seq1,
                            seq2,
                            match_score,
                            mismatch_score,
                            gap_penalty
                        )

                    # Local Alignment
                    else:

                        result = local_alignment(
                            seq1,
                            seq2,
                            match_score,
                            mismatch_score,
                            gap_penalty
                        )

                st.success(
                    "Alignment completed successfully!"
                )

                # Display results
                display_alignment_result(
                    result,
                    seq1,
                    seq2
                )

            except Exception as e:

                st.error(
                    f"An error occurred during alignment: {e}"
                )
