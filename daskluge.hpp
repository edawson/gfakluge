
#include <string>
#include <sstream>
#include <fstream>
#include <istream>
#include <map>
#include <vector>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <iostream>

using namespace std;
namespace gfak{

    struct custom_key{
        bool isdigit(const string s) const{
            //const char* s_str = s.c_str();
            string::const_iterator it;
            for (it = s.begin(); it != s.end(); it++){
                if (!std::isdigit(*it)){
                    return false;
                }
            }
            return true;
        }
        bool operator() (const string lhs, const string rhs) const{
            if (isdigit(lhs) && isdigit(rhs)){
                return atol(lhs.c_str()) < atol(rhs.c_str());
            }
            else{
                return lhs < rhs;
            }
        }
    };

    struct header_elem{
        std::string key;
        std::string type;
        std::string val;
    };

    struct opt_elem{
        std::string key;
        std::string type;
        std::string val;
    };

    struct annotation_elem{
        std::string key;
        std::string info;
    };

    struct segment_elem{
        std::string id;
        int len;
        std::string sequence;
    };

    struct fragment_elem{
        std::string segment_id;
        std::string id;
        int begin;
        int end;
    };

    struct edge_elem{
        std::string id;
        std::string segment_source_id;
        std::string segment_sink_id;
        int source_begin;
        int source_end;
        int sink_begin;
        int sink_end;
        string CIGAR;
    };

    struct gap_elem{
        std::string edge_id;
        std::string segment_first_id;
        std::string segment_second_id;
        int distance;
        int var;
    };

    struct group_elem{
        std::string id;
    };

    struct path_elem{
        std::string id;
        bool is_collapsed;
        vector<string> seg_ids;
        vector<string> edge_ids;
    };

    class GFAKluge{
        friend std::ostream& operator<<(std::ostream& os, GFAKluge& g);

        public:
        daskluge();
        ~daskluge();
        bool parse_das_file(std::string filename);
        bool parse_das_file(std::istream& gfa_stream);

        vector<header_elem> header_lines;

        //TODO: we should enforce graph structure,
        //i.e.:
        //1. Throw errors on dangling edges
        //2. Guarantee contained elements fall in valid sequences
        //3. Perhaps links and containeds should be added using methods like
        //  add_contained(contained_elem c)
        // to guarantee that they are actually added by their source.
        void add_link(string seq_name, link_elem);
        void add_link(sequence_elem s, link_elem l);
        void add_gap();
        void add_gap();
        void add_edge();
        void add_edge();
        void add_fragment();
        void add_fragment();
        void add_segment();
        void add_segment();
        void add_path(std::string pathname, path_elem path);
        void add_path();
        void collapse_paths();

        void set_version();

        vector<link_elem> get_links();
        vector<link_elem> get_links(string key);

        vector<edge_elem> get_edges();
        vector<edge_elem>(string key);

        vector<segment_elem> get_segments();
        vector<segment_elem> get_segments(string key);
        vector<segment_elem> get_neighboring_segments(string key, bool isFront);

        vector<path_elem> get_paths();
        vector< <vector<path_elem> > > get_path_vectors();
        vector<path_elem> get_paths_of_segment( string key);
        vector<path_elem> get_collapsed_paths();
        vector<path_elem> geT_collapsed_paths(vector<path_elem> paths);

        vector<header_elem> get_header();
        map<tag_type_t, string> get_header();

        map<string, sequence_elem, custom_key> get_name_to_seq();
        map<std::string, vector<link_elem> > get_seq_to_links();
        map<string, vector<path_elem> > get_seq_to_paths();
        map<string, vector<edge_elem> get_segment_to_edges();

        // TODO check whether writing to file is functional
        // Perhaps a write_gfa_file(string filename) method too?
        std::string to_string();
        std::string block_order_string();

        private:
        vector<header_elem> header;
        map<string, *header_elem> header_map;

        //Since we can't compare sequence elements for hashing,
        // we cheat and use their names (which are only sort of guaranteed to be
        // unique.
        std::vector<std::string> split(string s, char delim);
        std::string join(std::vector<std::string> splits, std::string glue);
        string header_string(map<string, header_elem>& opts);
        string opt_string(vector<opt_elem> opts);

    };

};
