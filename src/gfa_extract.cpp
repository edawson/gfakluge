#include <getopt.h>
#include <string>
#include <iostream>
#include "gfakluge.hpp"


using namespace std;
using namespace gfak;

void extract_help(char** argv){
    cout << "gfa_extract: extract a FASTA file from GFA" << endl
        << "Usage: gfa_extract [ -p ] <GFA_FILE> > file.fa " << endl
            << "-p / --include-paths  include paths in output (might be contigs?? We nevr know)" << endl
        << endl;
}

int main(int argc, char** argv){
    string gfa_file = "";
    bool include_paths = false;

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
            {"inlude-paths", no_argument, 0, 'p'},
            {0,0,0,0}
        };
    
        int option_index = 0;
        c = getopt_long(argc, argv, "hp", long_options, &option_index);
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
            case 'p':
                include_paths= true;
                break;
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
    // Do the same for groups/walks/paths, as the concated
    // sequences of the oriented seq_elems create a path (e.g. a chromosome or contig)
    if (include_paths){
        for (auto p : gg.get_groups()){
            stringstream pstr;
            if (p.second.ordered){
                
                pstr << ">" << p.second.id << endl;
                for (int i = 0; i < p.second.items.size(); i++){
                    string seqname = p.second.items[i];
                    bool ori = p.second.orientations[i];
                    pstr << ( ori ? seqs[seqname].sequence : string((seqs[seqname]).sequence.rbegin(), (seqs[seqname]).sequence.rend()));
                }
                pstr << endl;
                cout << pstr.str() << endl;
                pstr.str("");
            } 
        }
    }
    
	return 0;
}
