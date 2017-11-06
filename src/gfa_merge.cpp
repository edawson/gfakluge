#include "gfakluge.hpp"
#include <getopt.h>
#include <string>
#include <iostream>

using namespace gfak;
using namespace std;
/**
 * Merges a set of GFA files by:
 * 1. Coordinate their ID spaces (assume each one is an independent subgraph)
 * 2. Add both to a single GFAKluge instance. This will provide a sort.
 * 3. Push out the new GFA files to stdout.
 */

int main(int argc, char** argv){
    vector<string> g_files;
    bool block_order = false;

    if (argc == 1){
        cerr << "gfa_merge <gfa1.gfa> <gfa2.gfa> [<gfa3> ... <gfaN>]" << endl;
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
        c = getopt_long(argc, argv, "h", long_options, &option_index);
        if (c == -1){
            break;
        }

        switch (c){

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

    while(optind < argc){
        g_files.push_back(argv[optind]);
        optind++;
    }

    cerr << "Merging " << g_files.size() << " graphs..." << endl;

    // This does the same thing as IDs,
    // Just uses more memory...
    GFAKluge base;
    base.gfa_2_ize();
    for (auto gfi : g_files){
        gfak::GFAKluge gg;
        gg.parse_gfa_file(gfi);
        base.merge(gg);
    }
    base.gfa_2_ize();
    base.set_version(2.0);
    cout << base.to_string();
    cerr << "Done." << endl;


}
