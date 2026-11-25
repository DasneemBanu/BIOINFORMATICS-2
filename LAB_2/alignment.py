
# ============================================================
# LAB 2 - PAIRWISE SEQUENCE ALIGNMENT
# ============================================================


#Global Alignment (Needleman-Wunsch)
def global_alignment(seq1, seq2, match, mismatch, gap):

    rows = len(seq1)
    cols = len(seq2)

    # Create scoring matrix
    matrix = [
        [0 for j in range(cols + 1)]
        for i in range(rows + 1)
    ]

    # Create traceback matrix
    traceback = [
        [None for j in range(cols + 1)]
        for i in range(rows + 1)
    ]

    # Initialize first column
    for i in range(1, rows + 1):
        matrix[i][0] = i * gap
        traceback[i][0] = "U"

    # Initialize first row
    for j in range(1, cols + 1):
        matrix[0][j] = j * gap
        traceback[0][j] = "L"

    # Fill the scoring matrix
    for i in range(1, rows + 1):

        for j in range(1, cols + 1):

            # Calculate diagonal score
            if seq1[i - 1] == seq2[j - 1]:
                diagonal = matrix[i - 1][j - 1] + match
            else:
                diagonal = matrix[i - 1][j - 1] + mismatch

            # Calculate gap scores
            up = matrix[i - 1][j] + gap
            left = matrix[i][j - 1] + gap

            # Select maximum score
            best_score = max(diagonal, up, left)

            matrix[i][j] = best_score

            # Store traceback direction
            if best_score == diagonal:
                traceback[i][j] = "D"

            elif best_score == up:
                traceback[i][j] = "U"

            else:
                traceback[i][j] = "L"

    # Traceback starts from bottom-right
    i = rows
    j = cols

    aligned_seq1 = ""
    aligned_seq2 = ""

    traceback_path = []

    while i > 0 or j > 0:

        traceback_path.append((i, j))

        direction = traceback[i][j]

        if direction == "D":

            aligned_seq1 += seq1[i - 1]
            aligned_seq2 += seq2[j - 1]

            i -= 1
            j -= 1

        elif direction == "U":

            aligned_seq1 += seq1[i - 1]
            aligned_seq2 += "-"

            i -= 1

        elif direction == "L":

            aligned_seq1 += "-"
            aligned_seq2 += seq2[j - 1]

            j -= 1

    traceback_path.append((0, 0))

    # Reverse aligned sequences
    aligned_seq1 = aligned_seq1[::-1]
    aligned_seq2 = aligned_seq2[::-1]

    return {
        "score": matrix[rows][cols],
        "aligned_seq1": aligned_seq1,
        "aligned_seq2": aligned_seq2,
        "matrix": matrix,
        "traceback_path": traceback_path
    }



# Function 2D - Local Alignment (Smith-Waterman)
def local_alignment(seq1, seq2, match, mismatch, gap):

    rows = len(seq1)
    cols = len(seq2)

    # Create scoring matrix
    matrix = [
        [0 for j in range(cols + 1)]
        for i in range(rows + 1)
    ]

    # Create traceback matrix
    traceback = [
        [None for j in range(cols + 1)]
        for i in range(rows + 1)
    ]

    # Store highest score and its position
    max_score = 0
    max_position = (0, 0)

    # Fill the scoring matrix
    for i in range(1, rows + 1):

        for j in range(1, cols + 1):

            # Calculate diagonal score
            if seq1[i - 1] == seq2[j - 1]:
                diagonal = matrix[i - 1][j - 1] + match
            else:
                diagonal = matrix[i - 1][j - 1] + mismatch

            # Calculate gap scores
            up = matrix[i - 1][j] + gap
            left = matrix[i][j - 1] + gap

            # Select maximum score, including zero
            best_score = max(0, diagonal, up, left)

            matrix[i][j] = best_score

            # Store traceback direction
            if best_score == 0:
                traceback[i][j] = None

            elif best_score == diagonal:
                traceback[i][j] = "D"

            elif best_score == up:
                traceback[i][j] = "U"

            else:
                traceback[i][j] = "L"

            # Update maximum score
            if best_score > max_score:
                max_score = best_score
                max_position = (i, j)

    # Traceback starts from highest-scoring cell
    i, j = max_position

    aligned_seq1 = ""
    aligned_seq2 = ""

    traceback_path = []

    # Stop when a zero-score cell is reached
    while i > 0 and j > 0 and matrix[i][j] > 0:

        traceback_path.append((i, j))

        direction = traceback[i][j]

        if direction == "D":

            aligned_seq1 += seq1[i - 1]
            aligned_seq2 += seq2[j - 1]

            i -= 1
            j -= 1

        elif direction == "U":

            aligned_seq1 += seq1[i - 1]
            aligned_seq2 += "-"

            i -= 1

        elif direction == "L":

            aligned_seq1 += "-"
            aligned_seq2 += seq2[j - 1]

            j -= 1

    # Include zero-score stopping cell
    traceback_path.append((i, j))

    # Reverse aligned sequences
    aligned_seq1 = aligned_seq1[::-1]
    aligned_seq2 = aligned_seq2[::-1]

    return {
        "score": max_score,
        "aligned_seq1": aligned_seq1,
        "aligned_seq2": aligned_seq2,
        "matrix": matrix,
        "traceback_path": traceback_path
    }
