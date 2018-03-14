#include <getopt.h>
#include <cstring>
#include "gfakluge.hpp"

using namespace std;
using namespace gfak;


int gfakluge_package_help(char** argv){
    cerr << "GFAK: manipulate GFA files from the command line." << endl <<
        "Usage: " << argv[0] << " {command} [OPTIONS] <gfa_file> [...]" << endl <<
        "commands: " << endl <<
        "   convert: Convert between GFA 0.1 <-> 1.0 <-> 2.0" << endl <<
        "   diff:    Determine whether two GFA files have identical graphs" << endl <<
        "   extract: Convert the S lines of a GFA file to FASTA format." << endl <<
        "   ids:     Coordinate the ID spaces of multiple GFA graphs."  << endl <<
        "   concat:  Merge GFA graphs (without ID collisions)." << endl <<
        "   sort:    Print a GFA file in HSLP / HSEFGUO order." << endl <<
        "   stats:   Get assembly statistics (e.g. N50) for a GFA file." << endl <<
        "   subset:  Extract the subgraph between two IDs in a graph." << endl <<
        endl;
    return 1;

}

void extract_help(char** argv){
    cerr << argv[0] <<  " extract: extract a FASTA file from GFA" << endl
        << "Usage: " << argv[0] << "extract [ -p ] <GFA_FILE> > file.fa " << endl
            << "-p / --include-paths  include paths in output (might be contigs?? We nevr know)" << endl
        << endl;
}

void diff_help(char** argv){
    cerr << argv[0] << " diff: determine whether two GFA files are the same." << endl
      << "Usage: " << argv[0] << " diff [options] <GFA_File_1> <GFA_File_2>" << endl;
}

void convert_help(char** argv){
    cerr << argv[0] << " convert: convert a file between the various GFA formats." << endl
        << "Usage: " << argv[0] << " convert [options] <GFA_File>" << endl
        << "Options: " << endl
        << "  -s / --spec [0.1, 1.0, 2.0]   Convert the input GFA file to specification [0.1, 1.0, or 2.0]." << endl
        << "                                NB: not all GFA specs are backward/forward compatible, so a subset of the GFA may be used." << endl
        << "  -w / --walks   Output paths as walks, but maintain version (NOT SPEC COMPLIANT)." << endl
        << "  -p / --paths   Output walks as paths, but maintain version." << endl
        << "  -b / --block-order   Output GFA in block order [HSLP / HSLW | HSEFGUO]."
        << endl; 
}

void stats_help(char** argv){
    cerr << argv[0] << " stats: print assembly / graph stats for a GFA file." << endl
        << "Usage: " << argv[0] " stats [options] <GFA_File>" << endl
        << "   -a / --assembly  print assembly statistics (N50, N90, L50, L90)." << endl
        << "   -s / --all       print all graph statistics." << endl
        << "   -l / --length    print the total sequence length (in S lines)." << endl
        << "   -n / --num-nodes print the number of nodes." << endl
        << "   -e / --num-edges print the number of edges." << endl
        << "   -p / --paths     print some path statistics." << endl
        << endl;
}

/**
 * Output a fasta file from input GFA
 */
int extract_main(int argc, char** argv){
    string gfa_file = "";
    bool include_paths = false;

    if (argc < 3){
        cerr << "No GFA file given as input." << endl << endl; 
        extract_help(argv);
        exit(0);
    }

    int c;
    optind = 2;
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
                extract_help(argv);
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

/**
 * Return whether two GFA files are identical
 * in "structure" (number of links, segments, etC)
 */
int diff_main(int argc, char** argv){
    if (argc < 4){
        cerr << "diff requires two GFA files as input." << endl << endl;
        diff_help(argv);
        exit(0);
    }

    int c;
    optind = 2;
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
                diff_help(argv);
                exit(0);
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
        return -1;
    }
    else{
        // Check if all node IDs are the same.

        // Check if the edges of every node ID are the same in both graphs.

        // Check that the graphs both have identical paths.

        // Check fragments and other GFA2 specific items.

    }
    // Everything matches up, return 0
    cerr << "The graphs appear to be identical." << endl;
    
    return 0;
}

/**
 *  Convert a GFA file to another GFA version,
 *  FASTA, or Cytoscape
 */
int convert_main(int argc, char** argv){
    string gfa_file = "";
    bool block_order = false;
    double spec_version = 0.1;
    bool use_paths = true;

    if (argc < 3){
        cerr << "No GFA file provided. Please provide a GFA file to convert" << endl;
        convert_help(argv);
        exit(0);
    }

    int c;
    optind = 2;
    while (true){
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"block-order", no_argument, 0, 'b'},
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
            case '?':
            case 'h':
                convert_help(argv);
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


    gg.set_walks(!use_paths);

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

int ids_main(int argc, char** argv){

}

int merge_main(int argc, char** argv){

}

int sort_main(int argc, char** argv){

}

int stats_main(int argc, char** argv){
    
}


int subset_main(int argc, char** argv){

}


int main(int argc, char** argv){
    
    if (argc < 2){
       cerr << "No command provided. Please provide a second argument to GFA kluge to tell it what to do." << endl;
       return gfakluge_package_help(argv);
    }

    if (strcmp(argv[1], "convert") == 0){
        return convert_main(argc, argv);
    }
    else if (strcmp(argv[1], "diff") == 0){
        return diff_main(argc, argv);
    }
    else if (strcmp(argv[1], "extract") == 0){
        return extract_main(argc, argv);
    }
    else if (strcmp(argv[1], "ids") == 0){
        return ids_main(argv, argv); 
    }
    else if (strcmp(argv[1], "merge") == 0){
        return merge_main(argc, argv);
    }
    else if (strcmp(argv[1], "sort") == 0){
        return sort_main(argc, argv);
    }
    else if (strcmp(argv[1], "stats") == 0){
        return stats_main(argc, argv);
    }
    else if (strcmp(argv[1], "subset") == 0){
        return subset_main(argc, argv);
    }
    else {
        cerr << "No command " << '"' << argv[1] << '"' << endl;
        cerr << "Please use one of the valid gfak commands:" << endl;
        gfakluge_package_help(argv);
    }
    return 0;
}
