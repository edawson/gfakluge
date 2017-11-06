#include "gfakluge.hpp"
#include <getopt.h>
#include <string>
#include <iostream>


using namespace std;
using namespace gfak;
int main(int argc, char** argv){
    string gfa_file = "";
    bool block_order = false;

    if (argc == 1){
        cerr << "gfa_diff <a> <b> " << endl;
        exit(0);
    }

    int c;
    while (true){
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {0,0,0,0}
        };
    
        int option_index = 0;
        c = getopt_long(argc, argv, "hbi:", long_options, &option_index);
        if (c == -1){
            break;
        }

        switch (c){
            case '?':
            case 'h':
                cerr << "gfa_diff <gfa1> <gfa2> >> diff.gfa" << endl;
                exit(0);

            case 'b':
                block_order = true;
                break;

            default:
                abort();
        }
    }

    GFAKluge ff;
    GFAKluge gg;
    ff.parse_gfa_file(argv[optind]);
    optind++;
    gg.parse_gfa_file(argv[optind]);
    optind++;
    map<std::string, sequence_elem, custom_key> seq_1 = ff.get_name_to_seq();
    map<std::string, sequence_elem, custom_key> seq_2 = gg.get_name_to_seq();
    map<std::string, vector<edge_elem>> e_1 = ff.get_seq_to_edges();
    map<std::string, vector<edge_elem>> e_2 = ff.get_seq_to_edges();
    if (seq_1.size() != seq_2.size()){
        cerr << "Graphs have different numbers of sequences." << endl;
    }
    else{

    }
    
    return 0;
}
