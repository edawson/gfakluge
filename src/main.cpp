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

void ids_help(char** argv){
    cerr << argv[0] << " ids: coordinate the ID spaces of multiple GFA files." << endl
    << "Usage: " << argv[0] << " ids <GFA_FILE_1> .... <GFA_FILE_N>" << endl
    << "options: " << endl
    << "   -s / --start-ids   Start the relabeling process from <n_id:e_id:p_id>" << endl
    << "   -S / --spec <X>    Output GFA specification version <X>." << endl
    << "   -b / --block-order Output block-order (HSLCP) GFA." << endl;
}

void convert_help(char** argv){
    cerr << argv[0] << " convert: convert a file between the various GFA formats." << endl
        << "Usage: " << argv[0] << " convert [options] <GFA_File>" << endl
        << "Options: " << endl
        << "  -S / --spec [0.1, 1.0, 2.0]   Convert the input GFA file to specification [0.1, 1.0, or 2.0]." << endl
        << "                                NB: not all GFA specs are backward/forward compatible, so a subset of the GFA may be used." << endl
        << "  -w / --walks   Output paths as walks, but maintain version (NOT SPEC COMPLIANT)." << endl
        << "  -p / --paths   Output walks as paths, but maintain version." << endl
        << "  -b / --block-order   Output GFA in block order [HSLP / HSLW | HSEFGUO]."
        << endl; 
}

void sort_help(char** argv){
    cerr << argv[0] << " sort: sort a GFA file." << endl
        << "Usage: " << argv[0] << " stats [options] <GFA_File>" << endl
        << "Options:" << endl
        << "  -S / --spec [0.1, 1.0, 2.0]   Convert the input GFA file to specification [0.1, 1.0, or 2.0]." << endl
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
        << endl;
}

void subset_help(char** argv){
    cerr << argv[0] << " subset: extract a subset of a GFA graph between two node ids." << endl
    << "Usage: " << argv[0] << " subset [options] <gfa_file>" << endl
    << "Options:" << endl
    << "  -S / --spec <X>   GFA specification version for output." << endl
    << "  -s / --start-id  <n_id> Start ID of subgraph." << endl
    << "  -e / --end-id    <n_id> End ID of subgraph." << endl
    << "  -b / --block-order Output GFA in block order." << endl
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
            {"spec", required_argument, 0, 'S'},
            {0,0,0,0}
        };
    
        int option_index = 0;
        c = getopt_long(argc, argv, "hbpwS:", long_options, &option_index);
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
            case 'S':
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
    vector<string> g_files;
    bool block_order = false;
    string start_string;
    double spec = 0.0;

    if (argc == 1){
        cerr << "gfa_ids [-b -s S_ID:E_ID:F_ID:GA_ID:GR_ID ] <gfa1> <gfa2> ... [gfaN]" << endl;
        exit(0);
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
            {0,0,0,0}
        };

        int option_index = 0;
        c = getopt_long(argc, argv, "hbS:s:", long_options, &option_index);
        if (c == -1){
            break;
        }

        switch (c){
            case 'i':
                g_files.push_back( optarg );
                break;

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

}

int merge_main(int argc, char** argv){

}

int sort_main(int argc, char** argv){
    string gfa_file = "";
    bool block_order = true;
    double spec = 0.0;

    if (argc == 1){
        sort_help(argv);
        exit(0);
    }

    int c;
    optind = 2;
    while (true){
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"spec", no_argument, 0, 'S'},
            {0,0,0,0}
        };
    
        int option_index = 0;
        c = getopt_long(argc, argv, "hS:", long_options, &option_index);
        if (c == -1){
            break;
        }

        switch (c){

            case '?':
            case 'h':
                sort_help(argv);
                exit(1);
                break;

            case 'S':
                spec = stod(optarg);
                break;

            default:
                abort();
        }
    }
    gfa_file = argv[optind];

    GFAKluge gg;
    gg.parse_gfa_file(gfa_file);

    if (spec != 0){
        gg.set_version(spec);
    }
   
    cout << gg.block_order_string();
    
    
    return 0;
}

int stats_main(int argc, char** argv){
    string gfa_file = "";
    bool show_nodes = false;
    bool show_edges = false;
    bool show_containments = false;
    bool show_alignments = false;
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
                exit(3);

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
        exit(0);
    }

    int c;
    optind = 2;
    while (true){
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"block-order", no_argument, 0, 'b'},
            {"spec-version", required_argument, 0, 'S'},
            {"end-id", required_argument, 0, 'e'},
            {"start-id", required_argument, 0, 's'},
            {0,0,0,0}
        };
    
        int option_index = 0;
        c = getopt_long(argc, argv, "s:e:hS:b", long_options, &option_index);
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

        if (spec != 0.0){
            outg.set_version(spec);
        }
        cout << outg.to_string() << endl;
    }
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
    else {
        cerr << "No command " << '"' << argv[1] << '"' << endl;
        cerr << "Please use one of the valid gfak commands:" << endl;
        gfakluge_package_help(argv);
    }
    return 0;
}
