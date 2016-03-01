#include "gfakluge.hpp"
#include <getopt.h>
#include <string>
#include <iostream>


using namespace std;
using namespace gfak;
int main(int argc, char** argv){
    string gfa_file = "";
    bool show_nodes = false;
    bool show_edges = false;
    bool show_containments = false;
    bool show_alignments = false;
    bool show_length = false;

    if (argc == 1){
        cerr << "gfa_stats [ -n -e -s -l ] -i <GFA_File> " << endl;
        exit(0);
    }

    int c;
    while (true){
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"num-nodes", no_argument, 0, 'n'},
            {"num-edges", no_argument, 0, 'e'},
            {"length", no_argument, 0, 'l'},
            {"all", no_argument, 0, 's'},
            {"gfa-file", required_argument, 0, 'i'},
            {0,0,0,0}
        };
    
        int option_index = 0;
        c = getopt_long(argc, argv, "hnesli:", long_options, &option_index);
        if (c == -1){
            break;
        }

        switch (c){
            case 'i':
                gfa_file = optarg;
                break;

            case '?':
            case 'h':
                cerr << "gfa_stats [ -n -e -s -l ] -i <GFA_File> " << endl;
                exit(0);

            case 'e':
                show_edges = true;
                break;
            case 'n':
                show_nodes = true;
                break;
            case 's':
                show_nodes = true;
                show_length = true;
                show_edges = true;

            
            default:
                abort();
        }
    }

    GFAKluge gg;
    gg.parse_gfa_file(gfa_file);
    if (show_nodes){
        int num_nodes = gg.get_name_to_seq().size();
        cout << "Number of nodes: " << num_nodes << endl;
    }
    if (show_edges){
        
        int num_edges = 0;
        map<string, vector<link_ele> >::iterator it;
        for (gg.get_seq_to_link().size();
        cout << "Number of edges: " << num_edges << endl;
    }
    if (show_length){
        //This one's exciting. Let's iterate over the sequence elements and sum
        //the length of their sequence.
        map<string, sequence_elem, custom_key> my_seqs = gg.get_name_to_seq();
        map<string, sequence_elem, custom_key>::iterator it;
        int64_t total_len = 0;
        for (it = my_seqs.begin(); it != my_seqs.end(); it++){
            total_len += (it->second).sequence.size();
        }
        cout << "Total graph length in basepairs: " << total_len << endl;
    }
    

    
    return 0;
}
