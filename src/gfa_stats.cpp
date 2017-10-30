#include "gfakluge.hpp"
#include <getopt.h>
#include <string>
#include <iostream>


using namespace std;
using namespace gfak;

void stats_help(){
    cerr << "gfa_stats [ -a -n -e -s -l ] <GFA_File> " << endl;
    cerr << "  -a / --assembly    print assembly statistics (N50, L50, N90, L90)," << endl
         << "  -s / --all         print all graph statistics." << endl
         << "  -l / --length      print total sequence length." << endl
         << "  -n / --num-nodes   print number of nodes." << endl
         << "  -e / --num-edges   print number of edges." << endl
         << "  -p / --paths       print some path statistics." << endl
         << endl;
    exit(0);
}

int main(int argc, char** argv){
    string gfa_file = "";
    bool show_nodes = false;
    bool show_edges = false;
    bool show_containments = false;
    bool show_alignments = false;
    bool show_length = false;
    bool assembly_stats = false;
    bool show_paths = false;

    if (argc == 1){
        stats_help();
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
            {"paths", no_argument, 0, 'p'},
            {"assembly", no_argument, 0, 'a'},
            {0,0,0,0}
        };
    
        int option_index = 0;
        c = getopt_long(argc, argv, "hpanesl", long_options, &option_index);
        if (c == -1){
            break;
        }

        switch (c){

            case '?':
            case 'h':
                // nodes, edges, all stats, edges, paths
                stats_help();
                exit(3);

            case 'a':
                assembly_stats = true;
                break;

            case 'e':
                show_edges = true;
                break;
            case 'n':
                show_nodes = true;
                break;
            case 'l':
                show_length = true;
                break;
            case 's':
                show_nodes = true;
                show_length = true;
                show_edges = true;
                break;

            
            default:
                abort();
        }
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
        cout << "N50: " << gg.get_N50() << endl;
        cout << "N90: " << gg.get_N90() << endl;
        cout << "L50: " << gg.get_L50() << endl;
        cout << "L90: " << gg.get_L90() << endl;
    }
   

    
    return 0;
}
