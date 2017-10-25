#include "gfakluge.hpp"
#include <getopt.h>
#include <string>
#include <iostream>

using namespace std;
using namespace gfak;
int main(int argc, char** argv){
    vector<string> g_files;
    bool block_order = false;
    uint64_t start_id = 0;
    uint64_t end_id = UINT64_MAX;

    if (argc == 1){
        cerr << "gfa_subset -S start_id -E end_id <gfa.gfa> > sub.gfa" << endl;
        exit(0);
    }

    int c;
    while (true){
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"block-order", no_argument, 0, 'b'},
            {"spec-version", required_argument, 0, 's'},
            {0,0,0,0}
        };
    
        int option_index = 0;
        c = getopt_long(argc, argv, "S:E:hs:b", long_options, &option_index);
        if (c == -1){
            break;
        }

        switch (c){

            case '?':
            case 'h':
                cerr << "gfa_subset <gfa> > sub.gfa" << endl;
                exit(0);
            case 'S':
                start_id = stoul(optarg);
                break;
            case 'E':
                end_id = stoul(optarg);
                break;
            case 'b':
                block_order = true;
                break;

            default:
                abort();
        }
    }

    vector<string> gfiles;
    while (optind < argc){
        gfiles.push_back(argv[optind]);
        optind++;
    }

    for (auto i : gfiles){
        GFAKluge gg;
        GFAKluge outg;
        gg.parse_gfa_file(i);
        gg.gfa_2_ize();

        map<string, sequence_elem, custom_key> seqs = gg.get_name_to_seq();
        map<string, vector<edge_elem>> edges = gg.get_seq_to_edges();
        map<string, vector<gap_elem>> gaps = gg.get_seq_to_gaps();
        map<string, vector<fragment_elem>> fragments = gg.get_seq_to_fragments();
        
        for (auto s : seqs){
            if (stoul(s.first) <= end_id && stoul(s.first) >= start_id){
                outg.add_sequence(s.second);
                for (auto e : edges[s.first]){
                    if (stoul(e.sink_name) <= end_id && stoul(e.sink_name) >= start_id){
                        outg.add_edge(e.source_name, e);
                    }
                }
                for (auto gap : gaps[s.first]){
                    outg.add_gap(gap);
                }
                for (auto frag : fragments[s.first]){
                    outg.add_fragment(s.first, frag);
                }
            }
            
        }
        // TODO Just copy over all paths, but preferably trim the orientations and items vectors
        for (auto p : gg.get_name_to_path())
        {
            outg.add_path(p.first, p.second);
        }

        cout << outg.to_string_2() << endl;
    }


}
