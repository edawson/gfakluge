#include <getopt.h>
#include <cstring>
#include "gfakluge.hpp"
#include "gfa_builder.hpp"

using namespace std;
using namespace gfak;

void print_version_help(){
    cerr << "gfak (GFAKluge library) 1.1.1" << endl
        << "Copyright (C) 2015 Eric T. Dawson" << endl
        << "Licensed under the MIT License (https://opensource.org/licenses/MIT)" << endl
        << "This is free, open-source software: you are free to modify and redistribute it." << endl
        << "All code is provided as-is without warranty, as described in the MIT license." << endl
        << endl
        << "Written by Eric T. Dawson" << endl;
}

int gfakluge_package_help(char** argv){
    cerr << "GFAK: manipulate GFA files from the command line." << endl <<
        "Usage: " << argv[0] << " {command} [OPTIONS] <gfa_file> [...]" << endl <<
        "commands: " << endl <<
        "   convert: Convert between GFA 0.1 <-> 1.0 <-> 2.0" << endl <<
        "   diff:    Determine whether two GFA files have identical graphs" << endl <<
        "   extract: Convert the S lines of a GFA file to FASTA format." << endl <<
        "   fillseq: Add sequences from a FASTA file to S lines." << endl <<
        "   ids:     Coordinate the ID spaces of multiple GFA graphs."  << endl <<
        "   concat:  Merge GFA graphs (without ID collisions)." << endl <<
        "   sort:    Print a GFA file in HSLP / HSEFGUO order." << endl <<
        "   stats:   Get assembly statistics (e.g. N50) for a GFA file." << endl <<
        "   subset:  Extract the subgraph between two IDs in a graph." << endl <<
        "   trim:    Remove elements from a GFA graph." << endl <<
        endl;
    return 1;

}

void trim_help(char** argv){
    cerr << argv[0] <<  " trim: remove elements from a GFA graph." << endl
        << "Usage: " << argv[0] << " trim [OPTIONS] <GFA_FILE> " << endl
            << " -l / --length  <INT>  Remove segments (and their edges) if their sequence length is less than <INT>." << endl
            << " -n / --no-ambiguous   Remove segments which have ambiguous bases (i.e. non-ATGC) in their sequence." << endl
            << " -v / --version        print GFAK version and exit." << endl
        << endl;
}


void extract_help(char** argv){
    cerr << argv[0] <<  " extract: extract a FASTA file from GFA" << endl
        << "Usage: " << argv[0] << " extract [ -p ] <GFA_FILE> > file.fa " << endl
            << " -p / --include-paths  include paths in output" << endl
            << " -v / --version        print GFAK version and exit." << endl
        << endl;
}

void fillseq_help(char** argv){
    cerr << argv[0] << " fillseq: fill in S(equence) line sequences from a FASTA file." << endl
        << "Usage: " << argv[0] << " fillseq [options] -f <fasta.fa> <GFA_FILE> >> gfa_filled.gfa." << endl
        << "Options:" << endl
        << " -f / --fasta <f.fa>   {REQUIRED} a FASTA file containing sequences, with the GFA IDs as FASTA IDs." << endl
        << "                      Multiple FASTA files may be passed." << endl
        << " -S / --spec <SPEC>   Output in GFA version <SPEC>" << endl
        << endl;
}

void diff_help(char** argv){
    cerr << argv[0] << " diff: determine whether two GFA files are the same." << endl
      << "Usage: " << argv[0] << " diff [options] <GFA_File_1> <GFA_File_2>" << endl
            << " -v / --version        print GFAK version and exit." << endl
            << endl;
}

void ids_help(char** argv){
    cerr << argv[0] << " ids: coordinate the ID spaces of multiple GFA files." << endl
    << "Usage: " << argv[0] << " ids <GFA_FILE_1> .... <GFA_FILE_N>" << endl
    << "options: " << endl
    << "   -s / --start-ids   Start the relabeling process from <n_id:e_id:p_id>" << endl
    << "   -S / --spec <SPEC>    Output GFA specification version <X>." << endl
    << "   -b / --block-order Output block-order (HSLCP) GFA." << endl
    << "   -v / --version        print GFAK version and exit." << endl
    << endl;
}

void convert_help(char** argv){
    cerr << argv[0] << " convert: convert a file between the various GFA formats." << endl
        << "Usage: " << argv[0] << " convert [options] <GFA_File>" << endl
        << "Options: " << endl
        << "  -S / --spec <SPEC> [one of 0.1, 1.0, 2.0]   Convert the input GFA file to specification [0.1, 1.0, or 2.0]." << endl
        << "                                NB: not all GFA specs are backward/forward compatible, so a subset of the GFA may be used." << endl
        << "  -w / --walks   Output paths as walks, but maintain version (NOT SPEC COMPLIANT)." << endl
        << "  -p / --paths   Output walks as paths, but maintain version." << endl
        << "  -b / --block-order   Output GFA in block order [HSLP / HSLW | HSEFGUO]." << endl
        << "  -v / --version       print GFAK version and exit." << endl
        << "  -f / --fasta         print the S (sequence) elements in FASTA format." << endl
        << endl; 
}

void merge_help(char** argv){
    cerr << argv[0] << " merge: merge multiple GFA files into one structure." << endl
        << "Usage: " << argv[0] << " merge [options] <GFA_FILE_1> ... <GFA_FILE_N>" << endl
        << "Options: " << endl
        << "  -S / --spec <SPEC> [one of 0.1, 1.0, 2.0]   Convert the input GFA file to specification [0.1, 1.0, or 2.0]." << endl
        << "                                NB: not all GFA specs are backward/forward compatible, so a subset of the GFA may be used." << endl
        << "  -b / --block-order   Output GFA in block order [HSLP / HSLW | HSEFGUO]." << endl
        << "  -v / --version        print GFAK version and exit." << endl
        << endl; 
}

void build_help(char** argv){
    cerr << argv[0] << " build: generate a GFA variation graph from a FASTA reference and a VCF." << endl
    << "Usage: " << argv[0] << " build [options] -r <FASTA> -v <VCF>" << endl
    << "Options: " << endl
    << "  -m / --max-node-length <Int>  Maximum size for a node in basepairs (default: 1000)." << endl
    << "  -f / --fasta  <FASTA>     A fasta file to use for constructing the graph backbone." << endl
    << "  -v / --vcf <VCF>              A VCF file containing variants to put in the graph." << endl 
    << endl;
}

void sort_help(char** argv){
    cerr << argv[0] << " sort: sort a GFA file." << endl
        << "Usage: " << argv[0] << " sort [options] <GFA_File>" << endl
        << "Options:" << endl
        << "  -S / --spec <SPEC> [one of 0.1, 1.0, 2.0]   Convert the input GFA file to specification [0.1, 1.0, or 2.0]." << endl
        << "  -v / --version        print GFAK version and exit." << endl
        << endl;
}

void stats_help(char** argv){
    cerr << argv[0] << " stats: print assembly / graph stats for a GFA file." << endl
        << "Usage: " << argv[0] << " stats [options] <GFA_File>" << endl
        << "Options:" << endl
        << "   -a / --assembly  print assembly statistics (N50, N90, L50, L90)." << endl
        << "   -A / --all       print all graph statistics." << endl
        << "   -l / --length    print the total sequence length (in S lines)." << endl
        << "   -n / --num-nodes print the number of nodes." << endl
        << "   -e / --num-edges print the number of edges." << endl
        << "   -p / --paths     print some path statistics." << endl
        << "   -v / --version        print GFAK version and exit." << endl
        << endl;
}

void subset_help(char** argv){
    cerr << argv[0] << " subset: extract a subset of a GFA graph between two node ids." << endl
    << "Usage: " << argv[0] << " subset [options] <gfa_file>" << endl
    << "Options:" << endl
    << "  -s / --start-id  <n_id> Start ID of subgraph." << endl
    << "  -e / --end-id    <n_id> End ID of subgraph." << endl
    << "  -b / --block-order Output GFA in block order." << endl
    << "  -S / --spec <SPEC>   GFA specification version for output." << endl
    << "  -v / --version        print GFAK version and exit." << endl
    << endl;
}

/**
 * Trim segments (and their edges) from a graph
 */
int trim_main(int argc, char** argv){
    string gfa_file;

    int minlen = 0;
    bool no_amb = false;
    bool trim_paths = false;

    if (argc < 3){
        cerr << "No GFA file given as input." << endl << endl;
        trim_help(argv);
        exit(1);
    }
    int c;
    optind = 2;
    while (true){
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"trim-paths", no_argument, 0, 'p'},
            {"length", required_argument, 0, 'l'},
            {"no-ambiguous", no_argument, 0, 'N'},
            {"version", no_argument, 0, 'v'},
            {0,0,0,0}
        };
    
        int option_index = 0;
        c = getopt_long(argc, argv, "hvnpl:", long_options, &option_index);
        if (c == -1){
            break;
        }

        switch (c){
            case '?':
            case 'h':
                trim_help(argv);
                exit(0);
            case 'v':
                print_version_help();
                exit(0);
            case 'l':
                minlen = atoi(optarg);
                break;
            case 'n':
                no_amb = true;
                break;
            case 'p':
                trim_paths = true;
                cerr << "Path trimming is not yet implemented" << endl;
                exit(1);
            default:
                abort();
        }
    }

    if (optind > argc){
        cerr << "Error: no GFA file provided." << endl;
        exit(1);
    }
    gfa_file = argv[optind];
    
    GFAKluge gg;
    gg.parse_gfa_file(gfa_file);
    gg.trim_seqs(minlen, no_amb);
    cout << gg;

    return 0;
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
        exit(1);
    }

    int c;
    optind = 2;
    while (true){
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"include-paths", no_argument, 0, 'p'},
            {"version", no_argument, 0, 'v'},
            {0,0,0,0}
        };
    
        int option_index = 0;
        c = getopt_long(argc, argv, "hpv", long_options, &option_index);
        if (c == -1){
            break;
        }

        switch (c){
            case '?':
            case 'h':
                extract_help(argv);
                exit(0);
            case 'v':
                print_version_help();
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
            << name_seq.second.sequence << endl;
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

int build_main(int argc, char** argv){

    if (argc < 3){
        build_help(argv);
        exit(1);
    }
    
    std::string fasta_file;
    std::string vcf_file;
    char* insertion_fasta = NULL;
    int max_node_size = 128;

    double spec_version = 2.0;
    int c;
    optind = 2;
    while (true){
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"fasta", required_argument, 0, 'f'},
            {"spec", required_argument, 0, 'S'},
            {"vcf", required_argument, 0, 'v'},
            {"max-node-size", required_argument, 0, 'm'},
            {0,0,0,0}
        };
    
        int option_index = 0;
        c = getopt_long(argc, argv, "hm:i:f:S:v:", long_options, &option_index);
        if (c == -1){
            break;
        }

        switch (c){
            case '?':
            case 'h':
                build_help(argv);
                exit(0);
            case 'v':
                vcf_file = string(optarg);
                break;
            case 'i':
                insertion_fasta = optarg;
                break;
            case 'S':
                spec_version = stod(optarg);
                break;
            case 'f':
                fasta_file = string(optarg);
                break;
            case 'm':
                max_node_size = atoi(optarg);
                break;
            default:
                abort();
        }
    }
    gfak::GFAKluge gg;
    gg.set_version(spec_version);

    construct_gfa( (char*) fasta_file.c_str(), (char*) vcf_file.c_str(), (char*) insertion_fasta, gg, max_node_size);

    return 0;
}

int fillseq_main(int argc, char** argv){
    
    if (argc < 3){
        cerr << "fillseq requires a GFA file." << endl << endl;
        fillseq_help(argv);
        exit(1);
    }

    string fasta_file;
    string gfa_file;

    double spec_version = 0.0;
    bool block_order = false;

    int c;
    optind = 2;
    while (true){
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"fasta", required_argument, 0, 'f'},
            {"spec", required_argument, 0, 'S'},
            {"block-order", no_argument, 0, 'b'},
            {"version", no_argument, 0, 'v'},
            {0,0,0,0}
        };
    
        int option_index = 0;
        c = getopt_long(argc, argv, "hbf:S:v", long_options, &option_index);
        if (c == -1){
            break;
        }

        switch (c){
            case '?':
            case 'h':
                extract_help(argv);
                exit(0);
            case 'b':
                block_order = true;
                break;
            case 'v':
                print_version_help();
                exit(0);
            case 'S':
                spec_version = stod(optarg);
                break;
            case 'f':
                fasta_file = string(optarg);
                break;
            default:
                abort();
        }
    }

    if (optind >= argc){
        cerr << "Error: no GFA file provided." << endl;
        exit(1);
    }
    
    gfa_file = argv[optind];
    if (fasta_file.empty()){
        cerr << "Error: no FASTA file provided." << endl;
        exit(1);
    }
    GFAKluge gg;
    gg.parse_gfa_file(gfa_file);
    gg.fill_sequences(fasta_file.c_str());


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
        exit(1);
    }
    

    if (block_order){
        cout << gg.block_order_string();
    }
    else{
        cout << gg.to_string();
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
        exit(1);
    }

    int c;
    optind = 2;
    while (true){
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"version", no_argument, 0, 'v'},
            {0,0,0,0}
        };
    
        int option_index = 0;
        c = getopt_long(argc, argv, "hv", long_options, &option_index);
        if (c == -1){
            break;
        }

        switch (c){
            case '?':
            case 'h':
                diff_help(argv);
                exit(0);
            case 'v':
                print_version_help();
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
    else if (e_1.size() != e_2.size()){
        cerr << "Graphs have different numbers of edges" << endl;
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
    double spec_version = 2.0;
    bool use_paths = true;
    bool make_fasta = false;

    if (argc < 3){
        cerr << "No GFA file provided. Please provide a GFA file to convert" << endl;
        convert_help(argv);
        exit(1);
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
            {"spec", required_argument, 0, 'S'},
            {"version", no_argument, 0, 'v'},
            {"fasta", no_argument, 0, 'f'},
            {0,0,0,0}
        };
    
        int option_index = 0;
        c = getopt_long(argc, argv, "hvfbpwS:", long_options, &option_index);
        if (c == -1){
            break;
        }

        switch (c){
            case '?':
            case 'v':
                print_version_help();
                exit(0);
            case 'h':
                convert_help(argv);
                exit(0);

            case 'b':
                block_order = true;
                break;
            case 'S':
                spec_version = stod(optarg);
                break;
            case 'w':
                use_paths = false;
                break;
            case 'p':
                use_paths = true;
                break;
            case 'f':
                make_fasta = true;
                break;

            default:
                abort();
        }
    }
    gfa_file = argv[optind];

    GFAKluge gg;
    gg.parse_gfa_file(gfa_file);

    gg.set_walks(!use_paths);

    if (make_fasta){
        for (auto s : gg.get_name_to_seq()){
            cout << s.second.as_fasta_record() << endl; 
        }
        exit(0);
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
        exit(1);
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
    vector<string> g_files;
    bool block_order = false;
    string start_string;
    double spec = 0.0;

    if (argc == 1){
        cerr << "gfa_ids -s S_ID:E_ID:F_ID:GA_ID:GR_ID [options] <gfa1> <gfa2> ... [gfaN]" << endl;
        exit(1);
    }

    int c;
    optind = 2;
    while (true){
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"start-ids", required_argument, 0, 's'},
            {"spec", required_argument, 0, 'S'},
            {"blocker-order", no_argument, 0, 'b'},
            {"version", no_argument, 0, 'v'},
            {0,0,0,0}
        };

        int option_index = 0;
        c = getopt_long(argc, argv, "hvbS:s:", long_options, &option_index);
        if (c == -1){
            break;
        }

        switch (c){
            case 'v':
                print_version_help();
                exit(0);

            case '?':
            case 'h':
                ids_help(argv);
                exit(0);

            case 's':
                start_string = optarg;
                break;

            case 'b':
                block_order = true;
                break;
            case 'S':
                spec = stod(optarg);
                break;

            default:
                abort();
        }
    }

    while (optind < argc){
        g_files.push_back(argv[optind]);
        optind++;
    }

    int processed = 0;
    
    for (auto gfi : g_files){
        // get previous ID
        // if it is greater than the minimum ID in this gg,
        // increment all IDs in gg by the prev_id.
        // output this updated gg
        GFAKluge gg;
        gg.parse_gfa_file(gfi);
        gg.re_id(start_string);
        tuple<uint64_t, uint64_t, uint64_t, uint64_t, uint64_t> max_ids = gg.max_ids();
        stringstream xl;
        xl << std::get<0>(max_ids) << ":" << std::get<1>(max_ids) << ":" <<
            std::get<2>(max_ids) << ":" << std::get<3>(max_ids) << ":" <<
            std::get<4>(max_ids);
        start_string = xl.str();
        if (spec != 0){
            gg.set_version(spec);
        }
        cout << (block_order ? gg.block_order_string() : gg.to_string());
        ++processed;
        cerr << "Processed " << processed << " graphs..." << endl;
    }
    cerr << "Done." << endl;

    return 0;

}

int merge_main(int argc, char** argv){
    bool block_order = false;
    double spec = 0.0;
    vector<string> g_files;

    if (argc == 1){
        merge_help(argv);
        exit(1);
    }

    int c;
    optind = 2;
    while (true){
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"spec", required_argument, 0, 'S'},
            {"version", no_argument, 0, 'v'},
            {0,0,0,0}
        };
    
        int option_index = 0;
        c = getopt_long(argc, argv, "hvS:", long_options, &option_index);
        if (c == -1){
            break;
        }

        switch (c){

            case '?':
            case 'h':
                merge_help(argv);
                exit(0);
                break;
            case 'v':
                print_version_help();
                exit(0);
            case 'S':
                spec = stod(optarg);
                break;

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
    if (spec != 0){
        base.set_version(spec);
    }
    if (block_order){
        cout << base.block_order_string();
    }
    else{
        cout << base.to_string();

    }
    return 0;
}

int sort_main(int argc, char** argv){
    string gfa_file = "";
    bool block_order = true;
    double spec_version = 0.0;

    if (argc <= 2){
        sort_help(argv);
        exit(1);
    }

    int c;
    optind = 2;
    while (true){
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"spec", required_argument, 0, 'S'},
            {"version", no_argument, 0, 'v'},
            {0,0,0,0}
        };
    
        int option_index = 0;
        c = getopt_long(argc, argv, "hvS:", long_options, &option_index);
        if (c == -1){
            break;
        }

        switch (c){

            case '?':
            case 'h':
                sort_help(argv);
                exit(0);
                break;
            case 'v':
                print_version_help();
                exit(0);
                break;

            case 'S':
                spec_version = stod(optarg);
                break;

            default:
                abort();
        }
    }
    gfa_file = argv[optind];

    GFAKluge gg;
    gg.parse_gfa_file(gfa_file);

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
        exit(1);
    }

    //cout << gg.block_order_string();
    gg.output_to_stream(cout, block_order);
    
    return 0;
}

int stats_main(int argc, char** argv){
    string gfa_file = "";
    bool show_nodes = false;
    bool show_edges = false;
    bool show_length = false;
    bool assembly_stats = false;
    bool show_paths = false;
    bool all = true;

    if (argc <= 2){
        stats_help(argv);
        exit(1);
    }

    optind = 2;
    int c;
    while (true){
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"num-nodes", no_argument, 0, 'n'},
            {"num-edges", no_argument, 0, 'e'},
            {"length", no_argument, 0, 'l'},
            {"all", no_argument, 0, 'A'},
            {"paths", no_argument, 0, 'p'},
            {"assembly", no_argument, 0, 'a'},
            {"version", no_argument, 0, 'v'},

            {0,0,0,0}
        };
    
        int option_index = 0;
        c = getopt_long(argc, argv, "hpaAnel", long_options, &option_index);
        if (c == -1){
            break;
        }

        switch (c){

            case '?':
            case 'h':
                // nodes, edges, all stats, edges, paths
                stats_help(argv);
                exit(0);
            case 'v':
                print_version_help();
                exit(0);
            case 'a':
                assembly_stats = true;
                all = false;
                break;
            case 'e':
                show_edges = true;
                all = false;
                break;
            case 'n':
                show_nodes = true;
                all = false;
                break;
            case 'l':
                show_length = true;
                all = false;
                break;
            case 'A':
                all = true;
                break;
            case 'p':
                show_paths = true;
                all = false;
                break;

            default:
                abort();
        }
    }


    if (all){
        show_nodes = true;
        show_length = true;
        show_edges = true;
        assembly_stats = true;
    }
    gfa_file = argv[optind];
    GFAKluge gg;
    gg.parse_gfa_file(gfa_file);
    if (show_nodes){
        int num_nodes = gg.get_name_to_seq().size();
        cout << "Number of nodes: " << num_nodes << endl;
    }
    if (show_edges){
        
        uint64_t num_edges = 0;
        uint64_t num_links = 0;
        uint64_t num_contains = 0;
        // map<string, vector<link_elem> >::iterator it;
        // map<string, vector<link_elem> > s_to_l = gg.get_seq_to_link();
        // for (it = s_to_l.begin(); it != s_to_l.end(); it++){
        //     num_edges += (it->second).size();
        // }
        map<string, vector<edge_elem>> c_edge_map = gg.get_seq_to_edges();
        for (auto s : gg.get_name_to_seq()){
            vector<edge_elem> cvec = c_edge_map[s.first];
            for (auto e = cvec.begin(); e != cvec.end(); e++){
                num_edges++;
                int t = e->determine_type();
                if (t == 1){
                    num_links++;
                }
                else if (t == 2){
                    num_contains++;
                }
            }
        }
        
        cout << "Number of edges: " << num_edges << endl;
        cout << "Number of links: " << num_links << endl;
        cout << "Number of containments: " << num_contains << endl;
        
    }
    if (show_length){
        //This one's exciting. Let's iterate over the sequence elements and sum
        //the length of their sequence.
        map<string, sequence_elem, custom_key> my_seqs = gg.get_name_to_seq();
        map<string, sequence_elem, custom_key>::iterator it;
        int64_t total_len = 0;

        for (it = my_seqs.begin(); it != my_seqs.end(); it++){
            total_len += (it->second).length;
        }
        cout << "Total graph length in basepairs: " << total_len << endl;
    }
    
    if (show_paths){
        
    }

    if (assembly_stats){
        cout << "N50: " <<  (uint64_t) gg.get_N50() << endl;
        cout << "N90: " << (uint64_t) gg.get_N90() << endl;
        cout << "L50: " << (uint64_t) gg.get_L50() << endl;
        cout << "L90: " << (uint64_t) gg.get_L90() << endl;
    }
   

    
    return 0;
}


int subset_main(int argc, char** argv){
    vector<string> g_files;
    bool block_order = false;
    double spec = 0;
    uint64_t start_id = 0;
    uint64_t end_id = UINT64_MAX;

    if (argc == 1){
        subset_help(argv);
        exit(1);
    }

    int c;
    optind = 2;
    while (true){
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"block-order", no_argument, 0, 'b'},
            {"spec", required_argument, 0, 'S'},
            {"end-id", required_argument, 0, 'e'},
            {"start-id", required_argument, 0, 's'},
            {"version", no_argument, 0, 'v'},
            {0,0,0,0}
        };
    
        int option_index = 0;
        c = getopt_long(argc, argv, "vs:e:hS:b", long_options, &option_index);
        if (c == -1){
            break;
        }

        switch (c){

            case '?':
            case 'h':
                subset_help(argv);
                exit(0);
            case 'S':
                spec = stod(optarg);
                break;
            case 's':
                start_id = stoul(optarg);
                break;
            case 'e':
                end_id = stoul(optarg);
                break;
            case 'b':
                block_order = true;
                break;
            case 'v':
                print_version_help();
                exit(0);

            default:
                abort();
        }
    }

    if (argc < 3){
        cerr << "No GFA file provided" << endl;
        subset_help(argv);
        exit(1);
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

        if (spec != 0.0){
            outg.set_version(spec);
        }
        cout << outg.to_string() << endl;

    }
    return 0;
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
        return ids_main(argc, argv); 
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
    else if (strcmp(argv[1], "--version") == 0){
        print_version_help();
    }
    else if (strcmp(argv[1], "--help") == 0){
        gfakluge_package_help(argv);
    }
    else if (strcmp(argv[1], "fillseq") == 0){
        return fillseq_main(argc, argv);
    }
    else if (strcmp(argv[1], "build") == 0){
        return build_main(argc, argv);
    }
    else if (strcmp(argv[1], "trim") == 0){
        return trim_main(argc, argv);
    }
    else {
        cerr << "No command " << '"' << argv[1] << '"' << endl;
        cerr << "Please use one of the valid gfak commands:" << endl;
        gfakluge_package_help(argv);
    }
    return 0;
}
