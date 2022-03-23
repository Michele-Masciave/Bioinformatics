##
# LOCAL ALIGNMENT
#
import sys
import numpy as np

# Global variables wait to be initialized
match = 0
mismatch = 0
gap = 0


def fill_matrix(reference, read):
    score_matrix = np.zeros((len(read) + 1, len(reference) + 1))
    for col in range(1, len(subject) + 1):
        score_matrix[0][col] = 0
    for row in range(1, len(query) + 1):
        score_matrix[row][0] = 0
    return score_matrix


def get_similarity(diagonal_score, up_score, left_score, there_is_match):
    up_score += gap
    left_score += gap
    if there_is_match:
        diagonal_score += match
    else:
        diagonal_score += mismatch
    return max([max([up_score, left_score, diagonal_score]),0])


def print_alignment(_subject_alignment, _query_alignment):
    print('\nFinal alignment:')
    for b in _subject_alignment:
        print(b, end="")
    print()
    for _ in range(len(_subject_alignment)):
        print("|", end="")
    print()
    for b in _query_alignment:
        print(b, end="")
    print("\n\n")


if __name__ == '__main__':
    DNA_bases = ['A', 'T', 'C', 'G', 'N']

    subject = sys.argv[1].upper()
    for base in subject:
        if base not in DNA_bases:
            raise TypeError("DNA sequence not consistent.")
    query = sys.argv[2].upper()
    for base in subject:
        if base not in DNA_bases:
            raise TypeError("DNA sequence not consistent.")
    match = int(sys.argv[3])
    if type(match) is not int:
        raise TypeError("Only integers are allowed for match cost.")
    mismatch = int(sys.argv[4])
    if type(mismatch) is not int:
        raise TypeError("Only integers are allowed for mismatch cost.")
    gap = int(sys.argv[5])
    if type(gap) is not int:
        raise TypeError("Only integers are allowed for gap cost.")

    # Matrix headers definition
    matrix = fill_matrix(subject, query)

    # Matrix body filling
    for i in range(1, len(query) + 1):
        for j in range(1, len(subject) + 1):
            diagonal = matrix[i-1][j-1]
            up = matrix[i-1][j]
            left = matrix[i][j-1]
            wow_a_match = subject[j-1] == query[i-1]
            matrix[i][j] = get_similarity(diagonal, up, left, wow_a_match)

    # Get score on the right-bottom matrix value
    score = matrix[len(query)][len(subject)]

    # (1) Print score matrix
    print('\nLocal Alignment score:', score)
    print(matrix.T)

    # Find optimal path
    # -- get maximum from which  starting the traceback
    max_i = 0
    max_j = 0
    current_max = -1
    for row in range(len(query)):
        for col in range(len(subject)):
            if matrix[row][col] > current_max:
                current_max = matrix[row][col]
                max_i = row
                max_j = col
    # -- initialize some variables
    subject_alignment = []
    query_alignment = []
    up = 0
    left = 0
    diagonal = 0
    i = max_i
    j = max_j
    # -- start traceback!
    while matrix[i][j] != 0:
        if i > 0:
            up = matrix[i-1][j]
        if j > 0:
            left = matrix[i][j-1]
        if i > 0 and j > 0:
            diagonal = matrix[i-1][j-1]

        if (i > 0 and j > 0 and subject[j-1] == query[i-1]) or (diagonal >= up and diagonal >= left):
            # move diagonally
            # -- match!
            # -- mismatch convenience
            subject_alignment.append(subject[j-1])
            query_alignment.append(query[i-1])
            i -= 1
            j -= 1
        elif up >= diagonal and up >= left:
            # move up
            # -- subject (reference) gap convenience
            subject_alignment.append('-')
            query_alignment.append(query[i-1])
            i -= 1
        else:
            # move left
            # -- query(read) gap  convenience
            subject_alignment.append(subject[j-1])
            query_alignment.append('-')
            j -= 1

    # Reverse subject and query
    subject_alignment.reverse()
    query_alignment.reverse()

    # (2) Print alignment
    print_alignment(subject_alignment, query_alignment)
