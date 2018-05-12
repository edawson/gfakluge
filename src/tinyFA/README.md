tinyFA
-----------------
Parse, index and get random access to FASTA files
with no extra dependencies.

### Overview
tinyFA provides a (relatively) fast and highly minimalist header-only library
for reading FASTA files, especially in the case that you'd like random access
via their FAI indices. It requires only a modern (GCC4.8 or newer) C++ compiler.

### Build / install
Make sure both headers are in a path accessible to your code. You can either
copy them to your include directory, or pass the directory to your compiler
with ` -I/path/to/tinyFA`.

Then just include both headers in your C++ code and use the TFA namespace:

```
    #include "tinyFA/tinyfa.hpp"  
    #include "tinyFA/pliib.hpp"

    using namespace TFA;

```

### Usage

```

#include "tinyfa.hpp"
#include "pliib.hpp"


int main (int argc, char** argv){
    // usage: ./getseq <fastaFile> <seqName> <start> <end>
    // Parse a FASTA file and extract a subsequence.
    // If a FASTA index exists, use it, otherwise,
    // build one.
   
    // The tinyFA faidx struct
    tiny_faidx_t tf;
   
   // Check if an index exists, and create one if not.
    if (!checkFAIndexFileExists(argv[1])){
        createFAIndex(argv[1], tf);
    }
    else{
        // Parses an FAI file when passed a FASTA file name.
        parseFAIndex(argv[1], tf);
    }

    char* contigName = argv[2];
    int start = atoi(argv[3]);
    int end = atoi(argv[4]);
    char* seq;
    
    // Not passing start/end will return the whole contig.
    getSequence(tf, contigName, seq, start, end);
    cout << seq << endl;

    // seq gets allocated in getSequence, so you'll want to delete that.
    delete [] seq;

    return 0;
}
```

### Other tools
tinyFA takes a lot of code and inspiration from the
excellent [fastahack](https://github.com/ekg/fastahack).
Fastahack has been extensively used and we recommend it for production environments.
It differs from fastahack in that:

1. There's no default library setup for Fastahack. This could easily be remedied
   with some small makefile tweaks.  
2. tinyFA mostly uses structs and primitive types (rather than STL containers).


There's also [htslib](https://github.com/samtools/htslib), but if you want to
parse a FASTA file you have to build utilities for SAM/BAM/CRAM/VCF, it requires
zlib, lzma, and lz4 (which are often not up to date / installed by default).

However, htslib can parse files compressed with BGZF; this may be an important
addition when dealing with large files.

SeqLib and SeqAn are both excellent tools that add to htslib, with the same cons and
many extra pros.



