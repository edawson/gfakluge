#include <string>
#include <iostream>

#include "pliib.hpp"
#include "tinyfa.hpp"

using namespace std;
using namespace TFA;

int main(int argc, char** argv){


    tiny_faidx_t tf;
    if (argc < 2){
        cerr << "./getseq <fastaFile> <seq> [<start> <end>]" << endl;
        exit(1);
    }


        
        if (!checkFAIndexFileExists(argv[1])){
            createFAIndex(argv[1], tf);
        }
        else{
            parseFAIndex(argv[1], tf);
        }
    if (argc == 3){
        char* test = argv[2];
        char* seq;
        getSequence(tf, test, seq);
        cout << seq << endl;
        delete [] seq;
    }
    else if (argc == 5){
        char* test = argv[2];
        int start = atoi(argv[3]);
        int end = atoi(argv[4]);
        char* seq;
        getSequence(tf, test, seq, start, end);
        cout << seq << endl;
        delete [] seq;
    }
    else{
        cerr << "./getseq <fastaFile> <seq> [<start> <end>]" << endl;
        exit(1);
    }

    

    return 0;
}