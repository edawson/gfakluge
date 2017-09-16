/**
 * Coordinates the ID spaces of two GFA files (assumed to be independent graphs
 * Essentially, it either increments or renames each subsequent graph to prevent
 * name overlap in any element.
 *
 */

#include "gfakluge.hpp"
#include <getopt.h>
#include <string>
#include <iostream>
#include <vector>

using namespace std;
using namespace gfak;
int main(int argc, char** argv){
    vector<string> g_files;
    bool block_order = false;

    if (argc == 1){
        cerr << "gfa_ids [-s S_ID:F_ID:GA_ID:GR_ID ] <gfa1> <gfa2> ... [gfaN]" << endl;
        exit(0);
    }

    int c;
    while (true){
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"start-ids", no_argument, 0, 's'},
            {0,0,0,0}
        };

        int option_index = 0;
        c = getopt_long(argc, argv, "hs:", long_options, &option_index);
        if (c == -1){
            break;
        }

        switch (c){
            case 'i':
                g_files.push_back( optarg );
                break;

            case '?':
            case 'h':
                cerr << "gfa_merge [-b --block-order ] -i <GFA_File> -i <OTHER_GFA_FILE> >> my_sorted_gfa.gfa" << endl;
                exit(0);

            case 'b':
                block_order = true;
                break;

            default:
                abort();
        }
    }
    
    int processed = 0;
    tuple<uint64_t, uint64_t, uint64_t, uint64_t> max_ids = std::make_tuple(0,0,0,0);
    for (auto gfi : g_files){
        // get previous ID
        // if it is greater than the minimum ID in this gg,
        // increment all IDs in gg by the prev_id.
        // output this updated gg
        GFAKluge gg;
        gg.parse_gfa_file(gfi);
        gg.re_id(max_ids);
        max_ids = gg.max_ids();
        cout << gg.to_string();
        ++processed;
        cerr << "Processed " << processed << " graphs..." << endl;
    }
    cerr << "Done." << endl;

}

