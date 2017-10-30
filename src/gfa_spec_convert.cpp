#include <getopt.h>
#include <string>
#include <iostream>
#include "gfakluge.hpp"


using namespace std;
using namespace gfak;

void spec_help(char** argv){
    cout << "gfa_spec_convert: convert GFA 0.1 -> 1.0 -> 2.0, walks <-> paths, etc." << endl
        << "Usage: " << argv[0] << " [OPTIONS] <GFA_FILE> > gfa_out.gfa" << endl
        << "Options: " << endl
        << "  -w / --walks   Output paths as walks, but maintain version (NOT SPEC COMPLIANT)." << endl
        << "  -p / --paths   Output walks as paths, but maintain version." << endl
        << "  -s / --spec [0.1, 1.0, 2.0]   Convert the input GFA file to specification [0.1, 1.0, or 2.0]." << endl
        << "                                NB: not all GFA specs are backward/forward compatible, so a subset of the GFA may be used." << endl
        << "  -b / --block-order   Output GFA in block order [HLSP / HLSW | HSEFGUO]."
        << endl;
}

int main(int argc, char** argv){
    string gfa_file = "";
    bool block_order = false;
    double spec_version = 0.1;
    bool use_paths = true;

    if (argc < 2){
        spec_help(argv);
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
            {"paths", no_argument, 0, 'p'},
            {"walks", no_argument, 0, 'w'},
            {"spec", required_argument, 0, 's'},
            {0,0,0,0}
        };
    
        int option_index = 0;
        c = getopt_long(argc, argv, "hbpws:", long_options, &option_index);
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
            case 's':
                spec_version = stod(optarg);
                break;
            case 'w':
                use_paths = false;
                break;
            case 'p':
                use_paths = true;
                break;

            default:
                abort();
        }
    }
    gfa_file = argv[optind];

    GFAKluge gg;
    gg.parse_gfa_file(gfa_file);

    if (use_paths){
        gg.set_version(1.0);
    }
    else{
        gg.set_version(0.1);
        gg.set_walks(true);
    }

    if (spec_version == 0.1){
        gg.set_version(0.1);
    }
    else if (spec_version == 1.0){
        gg.set_version(1.0);
    }
    else if (spec_version == 2.0){
        gg.set_version(2.0);
    }
    else if (spec_version != 0.0){
        cerr << "Invalid specification number: " << spec_version << endl
        << "Please provide one of [0.1, 1.0, 2.0]." << endl;
        exit(22);
    }

    

    if (block_order){
        cout << gg.block_order_string();
    }
    else{
        cout << gg.to_string();
    }

	return 0;
}
