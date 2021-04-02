# Bioinformatics
## Needleman–Wunsch global alignment basic algorithm

### input

Program needs a subject (reference) string and a query (read) string as inputs. They can be provided via standard output (option INPUT), via text file (option INPUT_GENOME_FILE) or just running the program without any modality for debugging. The first one seems to be stable, text file read needs to be improved in terms of max length of string (need to edit code for building matrix on heap instead of on the stack). Gaps, matches and mismatches penalities can be modified.

### output

* Calculated score and alignment performed following the path built by the algorithm.
* option MATRIX enables to print matrix on console
* option STATISTICS enables matches, mismatches and gaps count.

### matrix construction

A matrix is initialized for evaluationg the path, matching query (vertical) and subject (horizontal) bases. The first row represents query gaps, first coloumn represents subject gaps. Element (0,0) has zero as value. Step-by-step is calcluated local similarity evaluating each (2x2) sub-matrix of the whole matrix, taking the minimum between:
- the "diagonal" value with respect to the current one +match is bases are equal or -mismatch if they are not equal;
- the right value ("up" with respect to the value to obtain) -gap
- the bottom value ("left" with respect to the value to obtain) -gap

### path construction

Starting with the lower-right element, until the element (0,0) is reached, this evaluation is performed:
- are bases related to this position equal? If yes, path moves diagonally, else:
- upper-value, left-value and diagonal-value are compared with the current one and the greater of them establishes the path to be taken (up, left, up-left)
