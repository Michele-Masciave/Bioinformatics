##
# GLOBAL ALIGNMENT
#
import numpy as np

# Global variables wait to be initialized
match = 0
mismatch = 0
gap = 0

# Just for pretty visualization
START_GREEN = "\033[0;32m"
END_GREEN = "\033[0m"
START_RED = "\033[0;31m"
END_RED = "\033[0m"


def fill_matrix(reference, read):
    score_matrix = np.zeros((len(read) + 1, len(reference) + 1))
    for col in range(1, len(subject) + 1):
        score_matrix[0][col] = gap * col
    for row in range(1, len(query) + 1):
        score_matrix[row][0] = gap * row
    return score_matrix


def get_similarity(diagonal_score, up_score, left_score, there_is_match):
    up_score += gap
    left_score += gap
    if there_is_match:
        diagonal_score += match
    else:
        diagonal_score += mismatch
    if diagonal_score >= up_score and diagonal_score >= left_score:
        return diagonal_score
    elif up_score >= diagonal_score and up_score >= left_score:
        return up_score
    return left_score


def print_alignment(_subject_alignment, _query_alignment):
    print('\n▬ SUBJECT: ', end="\t")
    for base in _subject_alignment:
        if base != '-':
            print(base, end="\t")
        else:
            print(START_RED + base + END_RED, end="\t")
    print("\n\t\t\t", end="")
    for index in range(len(_subject_alignment)):
        if _subject_alignment[index] == _query_alignment[index]:
            print(START_GREEN + "|" + END_GREEN, end="\t")
        else:
            print(START_RED + "|" + END_RED, end="\t")
    print('\n▬ QUERY:   ', end="\t")
    for base in _query_alignment:
        if base != '-':
            print(base, end="\t")
        else:
            print(START_RED + base + END_RED, end="\t")
    print()


if __name__ == '__main__':
    # User inserts subject and query
    subject = input("\n(1) Insert subject (reference): ")
    query = input("(2) Insert query (read): ")
    # User inserts scores
    match = int(input("(3) Insert match score: "))
    mismatch = int(input("(4) Insert mismatch score: "))
    gap = int(input("(5) Insert gap score: "))
    # subject = "AAAAA"
    # query = "AAAAAA"
    # match = 1
    # mismatch = 0
    # gap = -1
    n_matches = 0
    n_mismatches = 0
    n_gaps = 0

    # Matrix headers definition
    matrix = fill_matrix(subject, query)

    # Matrix body filling
    for i in range(1, len(query) + 1):
        for j in range(1, len(subject) + 1):
            diagonal = matrix[i - 1][j - 1]
            up = matrix[i - 1][j]
            left = matrix[i][j - 1]
            wow_a_match = subject[j - 1] == query[i - 1]
            matrix[i][j] = get_similarity(diagonal, up, left, wow_a_match)

    # (1) Print score matrix
    print()
    print(matrix)

    # Get score on the right-bottom matrix value
    score = int(matrix[len(query)][len(subject)])

    # Find optimal path
    subject_alignment = []
    query_alignment = []
    up = 0
    left = 0
    diagonal = 0
    i = len(query)
    j = len(subject)

    while i != 0 or j != 0:
        if i > 0:
            up = matrix[i-1][j]
        else:
            up = -2147483648        # just the minimum int
        if j > 0:
            left = matrix[i][j-1]
        else:
            left = -2147483648      # just the minimum int
        if i > 0 and j > 0:
            diagonal = matrix[i-1][j-1]
        else:
            diagonal = -2147483648  # just the minimum int

        if (i > 0 and j > 0 and subject[j-1] == query[i-1]) or (diagonal >= up and diagonal >= left):
            # move diagonally
            # -- match!
            # -- mismatch convenience
            if subject[j-1] == query[i-1]:
                n_matches += 1
            else:
                n_mismatches += 1
            subject_alignment.append(subject[j-1])
            query_alignment.append(query[i-1])
            i -= 1
            j -= 1
        elif up >= diagonal and up >= left:
            # move up
            # -- subject (reference) gap convenience
            n_gaps += 1
            subject_alignment.append('-')
            query_alignment.append(query[i-1])
            i -= 1
        else:
            # move left
            # -- query(read) gap  convenience
            n_gaps += 1
            subject_alignment.append(subject[j-1])
            query_alignment.append('-')
            j -= 1

    # Reverse subject and query
    subject_alignment.reverse()
    query_alignment.reverse()

    # (2) Print score and statistics
    print('\n▬ SCORE: %s (matches: +%d, mismatch: %d, gap: %d)' % (score, match, mismatch, gap))
    print("▬▬ Match percentage: {:.2f}%".format(n_matches / len(subject_alignment) * 100))
    print("▬▬ Mismatch percentage: {:.2f}%".format(n_mismatches / len(subject_alignment) * 100))
    print("▬▬ Gap percentage: {:.2f}%".format(n_gaps / len(subject_alignment) * 100))

    # (3) Print pretty alignment
    print_alignment(subject_alignment, query_alignment)