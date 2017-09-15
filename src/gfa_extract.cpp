#include <getopt.h>
#include <string>
#include <iostream>
#include "gfakluge.hpp"


using namespace std;
using namespace gfak;

void extract_help(char** argv){
    cout << "gfa_extract: extract a FASTA file from GFA" << endl
        << "Usage: gfa_extract <GFA_FILE> > file.fa " << endl
        << endl;
}

int main(int argc, char** argv){
    string gfa_file = "";

    if (argc < 2){
        extract_help(argv);
        exit(0);
    }

    int c;
    optind = 1;
    while (true){
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {0,0,0,0}
        };
    
        int option_index = 0;
        c = getopt_long(argc, argv, "h", long_options, &option_index);
        if (c == -1){
            break;
        }

        switch (c){
            case 'i':
                gfa_file = optarg;
                break;

            case '?':
            case 'h':
                cerr << "gfa_sort [-b --block-order ] -i <GFA_File> >> my_sorted_gfa.gfa" << endl;
                exit(0);
            default:
                abort();
        }
    }
    gfa_file = argv[optind];

    GFAKluge gg;
    gg.parse_gfa_file(gfa_file);

    map<string, sequence_elem, custom_key> seqs = gg.get_name_to_seq();

    for (auto name_seq : seqs){
        cout << ">" << name_seq.second.name << endl
            << name_seq.second.sequence << endl
            << endl;
    }
	return 0;
}
