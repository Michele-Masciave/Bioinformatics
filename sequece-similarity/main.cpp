/*
 * Masciavè Michele, global alignment algorithm
 * v1.0
 */

#include <iostream>
#include <fstream>

//just to customize console
#define CHECK_MARK "\033[0;32m\xE2\x9C\x94\033[0m "
#define UNCHECK_MARK "\033[0;31m\u2717\033[0m "
#define START_GREEN "\033[0;32m"
#define END_GREEN "\033[0m"
#define START_RED "\033[0;31m"
#define END_RED "\033[0m"

//! *** define your scores ***
#define gap -1
#define mismatch 0
#define match +1

//enable/disable options
#define INPUT 0
#define INPUT_GENOME_FILE 0
#define MATRIX 0
#define STATISTICS 1

/**
 * Reverse a given string
 * @param target
 * @param n
 * @param i
 */
void reverse_string(std::string& target, unsigned long n,int i){
    if(n<=i) return;
    std::swap(target[i],target[n]);
    reverse_string(target,n-1,i+1);
}

/**
 * Calculate similarity at a given step on the path
 * @param diagonal
 * @param up
 * @param left
 * @param there_is_match
 * @return
 */
int max(int diagonal, int up, int left, bool there_is_match){
    up = up + gap;
    left = left + gap;
    diagonal = there_is_match ? diagonal + match : diagonal + mismatch;
    if (diagonal >= up && diagonal >= left) return diagonal;
    if (up >= diagonal && up >= left) return up;
    return left;
}

int main() {
    //declarations and initialization
    unsigned int i,j;
    int score = 0;
    unsigned int n_matches = 0;
    unsigned int n_mismatches = 0;
    unsigned int n_gaps = 0;
    std::string subject;
    std::string query;
#if INPUT
    //user input via console
    std::cout << "1. Insert or paste reference/subject: ";
    std::cin >> std::uppercase >> subject;
    std::cout << "2. Insert or paste read/query: ";
    std::cin >> std::uppercase >> query;
#elif INPUT_GENOME_FILE
    //file input
    std::ifstream genome;
    genome.open("../subject.txt", std::ios::beg);
    if (genome.is_open()) {
        std::string line;
        while (getline(genome,line)) {
            subject.append(line);
        }
    } else {
        genome.close();
        return -1;
    }
    genome.close();

    std::ifstream read;
    read.open("../query.txt", std::ios::beg); //modified genome
    if (read.is_open()) {
        std::string line;
        while (getline(read,line)) {
            query.append(line);
        }
    } else {
        read.close();
        return -1;
    }
    read.close();
#else
    //direct input for debug
    subject = "GGGCGGCGACCTCGCGGGTTTTCGCTATTTATGAAAATTTTC"; //"GACTAC";
    query = "CGGCGGCGACCTCGCGGGTTTTCGTATTTATGAA"; //ACGC;
#endif

    const unsigned int N = subject.size()+1;
    const unsigned int M = query.size()+1;
    int matrix[N][M];

    //fill horizontal line
    for(j=0; j<N; j++) matrix[0][j] = gap*j;

    //fill vertical line
    for(i=1; i<M; i++) matrix[i][0] = gap*i;

    //fill matrix
    for(i=1; i<M; i++)
        for(j=1; j<N; j++){
            int similarity = max(
                    matrix[i-1][j-1],
                    matrix[i-1][j],
                    matrix[i][j-1],
                    subject.at(j-1) == query.at(i-1));
            matrix[i][j] = similarity;
        }
    score = matrix[i-1][j-1];

#if MATRIX
    //print matrix
    std::cout << "\n\t\t"; for(char c : subject) std::cout <<"\t\t" << c << " "; std::cout << std::endl;
    for(i=0; i<M; i++){
        if(i>0) std::cout << query.at(i-1) << "\t";
        else std::cout << "\t";
        for(j=0; j<N; j++)
            std::cout << "\t" <<  matrix[i][j] << "\t";
        std::cout << std::endl;
    }
#endif

    //find optimal path
    std::string subject_alignment;
    std::string query_alignment;
    int up=0, diagonal = 0, left=0;
    i=M-1, j=N-1;

    while(i!=0 || j!=0) {
       if(i>0) up=matrix[i-1][j]; else up=INT_MIN; //-2147483648
       if(j>0)left=matrix[i][j-1]; else left=INT_MIN;
       if(i>0 && j>0) diagonal=matrix[i-1][j-1]; else diagonal=INT_MIN;

       if((i>0 && j>0 && subject.at(j-1) == query.at(i-1)) || (diagonal >= up && diagonal >= left)) {
           //move diagonally
           // - match! :-)
           // - mismatch convenience
#if STATISTICS
           if(subject.at(j-1) == query.at(i-1))
               n_matches++;
           else n_mismatches++;
#endif
           subject_alignment.push_back(subject.at(j-1));
           query_alignment.push_back(query.at(i-1));
           i--;
           j--;
       }  else if(up >= diagonal && up >= left) {
           //move up
           // - subject (reference) gap convenience
#if STATISTICS
            n_gaps++;
#endif
           subject_alignment.push_back('-');
           query_alignment.push_back(query.at(i-1));
           i--;
       } else {
           //move left
           // - query (read) gap convenience
#if STATISTICS
           n_gaps++;
#endif
           subject_alignment.push_back(subject.at(j-1));
           query_alignment.push_back('-');
           j--;
       }
    }

    //output result
    reverse_string(subject_alignment, subject_alignment.size()-1,0);
    reverse_string(query_alignment, query_alignment.size()-1,0);

    std::cout << "\nSCORE: " << score <<
    " (matches: +" << match << ", mismatches: " << mismatch << ", gaps: " << gap << ")" <<  std::endl;

#if STATISTICS
    std::cout << CHECK_MARK " matches: " << n_matches << "/" << subject_alignment.size() << std::endl;
    std::cout << UNCHECK_MARK " mismatches: " << n_mismatches << "/" << subject_alignment.size() << std::endl;
    std::cout << UNCHECK_MARK " gaps: " << n_gaps << "/" << subject_alignment.size() << std::endl << std::endl;
#endif

    std::cout << "S:\t";
    for(char base : subject_alignment){
        if(base != '-')
            std::cout << base << "\t";
        else
            std::cout << START_RED << base << END_RED << "\t";
    }

    std::cout<< std::endl << "\t";
    for(i=0; i<subject.size(); i++) {
        if(subject_alignment.at(i) == query_alignment.at(i))
            std::cout << START_GREEN "|" END_GREEN "\t";
        else
            std::cout << START_RED "|" END_RED "\t";
    }

    std::cout << "\nQ:\t";
    for(char base : query_alignment){
        if(base != '-')
            std::cout << base << "\t";
        else
            std::cout << START_RED << base << END_RED << "\t";
    }

    std::cout << std::endl;
    return 0;
}