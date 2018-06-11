#include <string>
#include <iostream>

#include "pliib.hpp"
#include "tinyfa.hpp"

using namespace std;
using namespace TFA;

int main(int argc, char** argv){

    if (argc < 2){
        cerr << "./index <fastaFile>" << endl;
    }

    tiny_faidx_t tf;
    createFAIndex(argv[1], tf);
    writeFAIndex(argv[1], tf);

    return 0;
}