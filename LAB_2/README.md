
@"
# Bioinformatics Lab 2 - Pairwise Sequence Alignment

## Overview

This project combines Bioinformatics Lab 1 and Lab 2 into one Streamlit application.

The application provides protein retrieval and basic protein analysis, together with pairwise sequence alignment using Global Alignment and Local Alignment.

## Features

### Lab 1 - Protein Analysis
- Retrieve protein sequences using NCBI Entrez.
- Display protein ID, name, description and sequence.
- Calculate sequence length.
- Calculate amino acid composition.
- Calculate molecular weight.
- Calculate isoelectric point.

### Lab 2 - Pairwise Sequence Alignment
- Accept two input sequences.
- Perform Global Alignment (Needleman-Wunsch).
- Perform Local Alignment (Smith-Waterman).
- Allow user-defined match, mismatch and gap scores.
- Display alignment score and aligned sequences.
- Display the scoring matrix.
- Highlight the traceback path.
- Display alignment statistics.

## Requirements

- Python
- Streamlit
- Biopython
- Pandas

## Installation

1. Open the project folder in VS Code.
2. Create a virtual environment:

   ```bash
   python -m venv .venv
   ```

3. Activate the virtual environment.

   Windows PowerShell:

   ```powershell
   .venv\Scripts\Activate.ps1
   ```

4. Install the required packages:

   ```bash
   pip install -r requirements.txt
   ```

## Running the Application

Run the following command:

```bash
streamlit run app.py
```

The application will open in your web browser.

## Project Structure

```text
LAB_2/
├── app.py
├── alignment.py
├── requirements.txt
├── README.md
└── .gitignore
```

## Implementation

The pairwise alignment algorithms are implemented manually using dynamic programming.

- Global Alignment: Needleman-Wunsch algorithm.
- Local Alignment: Smith-Waterman algorithm.

No external sequence alignment library is used for the alignment calculations.
"@ | Set-Content README.md
