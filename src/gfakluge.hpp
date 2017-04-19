#ifndef GFAK_HPP
#define GFAK_HPP
#include <string>
#include <sstream>
#include <fstream>
#include <istream>
#include <algorithm>
#include <map>
#include <vector>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <iomanip>

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

    struct comment_elem{
        std::vector<string> comments;
    };

    struct header_elem{
        std::string key;
        std::string type;
        std::string val;
        // TODO make the val a vector<string>; this should handle every case
        // essentially, a single header_elem will encode all lines beginning with that key.
        // e.g. multiple program lines.
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

    struct path_elem{
        std::string name;
        vector<string> segment_names;
        vector<bool> orientations;
        vector<string> overlaps;
    };

    struct walk_elem{
	std::string source_name;
	std::string name;
	long rank;
	bool is_reverse;
	std::string cigar = "*";
    };

    struct das_path_elem{
        std::string name;
        std::vector<string> segments;
        std::vector<string> overlaps;
    };

	struct alignment_elem{
				std::string source_name;
				int position;
				std::string ref;
				bool source_orientation_forward;
				int length;
                map<string, string> opt_fields;
	};

    struct sequence_elem{
        std::string sequence;
        std::string name;
        vector<opt_elem> opt_fields;
        long id;

    };

    struct link_elem{
        std::string source_name;
        std::string sink_name;
        bool source_orientation_forward;
        bool sink_orientation_forward;
        std::string cigar;
        map<string, string> opt_fields;
    };

    struct contained_elem{
        std::string source_name;
        std::string sink_name;
        bool source_orientation_forward;
        bool sink_orientation_forward;
        int pos;
        std::string cigar;
        map<string, string> opt_fields;
    };


    class GFAKluge{
			friend std::ostream& operator<<(std::ostream& os, GFAKluge& g);

        public:
            GFAKluge();
            ~GFAKluge();
            bool parse_gfa_file(std::string filename);
            bool parse_gfa_file(std::istream& gfa_stream);

            //TODO: we should enforce graph structure,
            //i.e.:
            //1. Throw errors on dangling edges
            //2. Guarantee contained elements fall in valid sequences
            //3. Perhaps links and containeds should be added using methods like
            //  add_contained(contained_elem c)
            // to guarantee that they are actually added by their source.
            void add_link(string seq_name, link_elem);
            void add_link(sequence_elem s, link_elem l);
            void add_contained(string seq_name, contained_elem c);
            void add_contained(sequence_elem s, contained_elem c);
            void add_alignment(string s, alignment_elem a);
            void add_alignment(sequence_elem s, alignment_elem a);
            void add_sequence(sequence_elem s);
            void add_path(std::string pathname, path_elem path);
            void add_walk(std::string walkname, walk_elem w);
            double get_version();
            void set_version(double version);
            void set_version();
            void set_walks(bool ws);

            vector<link_elem> get_links(sequence_elem seq);
            vector<link_elem> get_links(string seq_name);

            vector<contained_elem> get_contained(string seq_name);
            vector<contained_elem> get_contained(sequence_elem seq);

            vector<alignment_elem> get_alignments(string seq_name);
            vector<alignment_elem> get_alignments(sequence_elem seq);

            map<string, header_elem> get_header();
            map<string, sequence_elem, custom_key> get_name_to_seq();
            map<std::string, vector<link_elem> > get_seq_to_link();
            map<std::string, vector<contained_elem> > get_seq_to_contained();
            map<std::string, vector<alignment_elem> > get_seq_to_alignment();
            map<string, vector<walk_elem> > get_seq_to_walks();
            map<string, path_elem> get_name_to_path();

            void paths_as_walks();
            void walks_as_paths();

            // TODO check whether writing to file is functional
            // Perhaps a write_gfa_file(string filename) method too?
            std::string to_string();
            std::string block_order_string();


        private:
            bool use_walks = false;
            // Store whether we've already gone walks->paths and paths->walks
            bool normalized_paths = false;
            bool normalized_walks = false;
            double version = 0.0;
            map<std::string, header_elem> header;
            map<std::string, vector<contained_elem> > seq_to_contained;
            map<std::string, vector<link_elem> > seq_to_link;
            //Since we can't compare sequence elements for hashing,
            // we cheat and use their names (which are only sort of guaranteed to be
            // unique.
            map<string, sequence_elem, custom_key> name_to_seq;
            std::vector<std::string> split(string s, char delim);
            std::string join(std::vector<std::string> splits, std::string glue);
            map<string, vector<alignment_elem> > seq_to_alignment;
            string header_string(map<string, header_elem>& opts);
            string opt_string(vector<opt_elem> opts);
            
            map<string, vector<walk_elem> > seq_to_walks;
	        map<string, path_elem> name_to_path;
    };

};
#endif
