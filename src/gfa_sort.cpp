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
        cerr << "gfa_sort [-b ]  <GFA_File> >> my_sorted_gfa.gfa" << endl;
        exit(0);
    }

    int c;
    optind = 1;
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
                gfa_file = optarg;
                break;

            case '?':
            case 'h':
                cerr << "gfa_sort [-b --block-order ] -i <GFA_File> >> my_sorted_gfa.gfa" << endl;
                exit(0);

            case 'b':
                block_order = true;
                break;

            default:
                abort();
        }
    }
    gfa_file = argv[optind];

    GFAKluge gg;
    gg.parse_gfa_file(gfa_file);

    if (block_order){
        cout << gg.block_order_string();
    }
    else{
        cout << gg;
    }
    
    return 0;
}
