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


int main(int argc, char** argv){
    vector<string> g_files;
    bool block_order = false;

    if (argc == 1){
        cerr << "gfa_sort [-b ] -i <GFA_File> >> my_sorted_gfa.gfa" << endl;
        exit(0);
    }

    int c;
    while (true){
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"block-order", no_argument, 0, 'b'},
            {"gfa-file", required_argument, 0, 'i'},
            {0,0,0,0}
        };

        int option_index = 0;
        c = getopt_long(argc, argv, "hbi:", long_options, &option_index);
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
    
    int last_prev_id = 0;
    for (auto gfi : g_files){
        GFAKluge gg;
        gg.parse_gfa_file(gfi);

        // get previous ID
        // if it is greater than the minimum ID in this gg,
        // increment all IDs in gg by the prev_id.
        // output this updated gg
        cout << gg.to_string();
    }

}

