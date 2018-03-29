#ifndef GFAK_HPP
#define GFAK_HPP

#include <string>
#include <sstream>
#include <istream>
#include <algorithm>
#include <map>
#include <vector>
#include <sstream>
#include <cstdlib>
#include <iostream>
#include <bitset>

namespace gfak{
    
    // ( <header> | <segment> | <fragment> | <edge> | <gap> | <group> )

    enum gfa_line_types {HEADER_LINE,SEGMENT_LINE,FRAGMENT_LINE,EDGE_LINE,GAP_LINE,GROUP_LINE};
    
    struct custom_key{
        bool isdigit(const std::string &s) const {
            //const char* s_str = s.c_str();
            std::string::const_iterator it;
            for (it = s.begin(); it != s.end(); it++){
                if (!std::isdigit(*it)){
                    return false;
                }
            }
            return true;
        }
        bool operator() (const std::string &lhs, const std::string &rhs) const{
            if (isdigit(lhs) && isdigit(rhs)) {
                return atol(lhs.c_str()) < atol(rhs.c_str());
            }
            else{
                return lhs < rhs;
            }
        }
    };

    struct comment_elem {
        std::vector<std::string> comments;
    };

    struct header_elem {
        std::string key;
        std::string type;
        std::string val;
        // TODO make the val a vector<std::string>; this should handle every case
        // essentially, a single header_elem will encode all lines beginning with that key.
        // e.g. multiple program lines.
    };

		struct opt_elem {
        std::string key;
        std::string type;
        std::string val;
        std::string to_string() {
            std::stringstream st;
            st << key << ":" << type << ":" << val;
            return st.str();
        }
    };

    struct annotation_elem{
        std::string key;
        std::string info;
    };

    struct path_elem{
        std::string name;
        std::vector<std::string> segment_names;
        std::vector<bool> orientations;
        std::vector<std::string> overlaps;
        std::map<std::string, opt_elem> opt_fields;
        std::string to_string() {
            return "";
        }
    };

    struct walk_elem{
	std::string source_name;
	std::string name;
	long rank;
	bool is_reverse;
	std::string cigar = "*";
    };


	struct alignment_elem{
				std::string source_name;
				int position;
				std::string ref;
				bool source_orientation_forward;
				int length;
        std::map<std::string, std::string> opt_fields;
	};

    struct sequence_elem{
        std::string sequence;
        std::string name;
        uint64_t length = UINT64_MAX;
        std::vector<opt_elem> opt_fields;
        long id;
        std::string to_string_2(){
            std::stringstream st;
            st << "S" << "\t" << name << "\t" << length << "\t" << sequence;
            if (opt_fields.size() > 0){
                for (auto i : opt_fields){
                    st << "\t" << i.to_string();
                }
            }
            return st.str();
        }
        std::string to_string_1(){
            std::stringstream st;
            st << "S" << "\t" << name << "\t" << sequence;
            if (opt_fields.size() > 0){
                for (auto i : opt_fields){
                    st << "\t" << i.to_string();
                }
            }
            return st.str();
        }
    };

    struct link_elem{
        std::string source_name;
        std::string sink_name;
        bool source_orientation_forward;
        bool sink_orientation_forward;
        std::string cigar;
        std::map<std::string, opt_elem> opt_fields;
        std::string to_string(){
            std::stringstream st;
            st << "L" 
               << "\t" << source_name << "\t" 
               << (source_orientation_forward ? "+" : "-") << "\t" << sink_name << "\t" 
               << (sink_orientation_forward ? "+" : "-") << "\t" << cigar; 

             if (opt_fields.size() > 0){
                for (auto i : opt_fields){
                    st << "\t" << i.second.to_string();
                }
            }
            return st.str();
        }
    };

    struct contained_elem {
        std::string source_name;
        std::string sink_name;
        bool source_orientation_forward;
        bool sink_orientation_forward;
        uint64_t pos;
        std::string cigar;
        std::map<std::string, opt_elem> opt_fields;
    };


    // These elements, along with the <length> field in sequence_elem,
    // are all that's needed to parse to GFA2
    struct edge_elem{
        std::string id = "*";
        // 0: unset, 1: link, 2: containment, 3: generic edge or both C/L set (impossible)
        int type;
        std::string source_name;
        std::string sink_name;
        bool source_orientation_forward;
        bool sink_orientation_forward;
        std::bitset<4> ends;
        uint64_t source_begin = 0;
        uint64_t source_end = 0;
        uint64_t sink_begin = 0;
        uint64_t sink_end = 0;
        std::string alignment;
        std::map<std::string, opt_elem> tags;
        int determine_type(){
                    if (type == 1 || type == 2){
                        return type;
                    }
                    // Determine if an edge is a link or a containment and fill in the right fields.
                    //else if (ends.test(0) && ends.test(1) && !ends.test(2) && !ends.test(3) ){
                    else if (ends.test(0) && ends.test(1) && sink_begin == 0 && sink_end == 0 ){
                        type = 1;
                        return 1;
                    }
                    else if (!ends.test(0) && !ends.test(2)){
                        type = 2;
                        return 2;
                    }
                    else{
                        type = 3;
                        return 3;
                    }
                
        }
        std::string to_string_2(){
            std::stringstream st;
            st << "E" << "\t" << id << "\t" <<
                source_name << (source_orientation_forward ? "+" : "-") <<
                 "\t" << sink_name << (sink_orientation_forward ? "+" : "-") << "\t" <<
                std::to_string(source_begin);
                if (ends.test(0)){
                     st << "$";
                }
                st << "\t" << std::to_string(source_end);
                if (ends.test(1)){
                     st << "$";
                }
                st << "\t" << std::to_string(sink_begin);
                if (ends.test(2)){
                    st << "$";
                }
                st << "\t" << std::to_string(sink_end);
                if (ends.test(3)){
                    st << "$";
                }
                st << "\t" << alignment;
                for (auto i : tags){
                    st << "\t" << i.second.to_string();
                }

                return st.str();
        }
        std::string to_string_1() {
            int t = determine_type();
            if (t > 2 || t == 0){
                std::cerr << "warning: unexpressable edge \"" << to_string_2()  << "\"" << std::endl
                     << "will not appear in stdout." << std::endl;
            }
            std::stringstream st;
            st << (t == 1 ? "L" : "C") << "\t" << source_name << "\t" << 
                (source_orientation_forward ? "+" : "-") << "\t" <<
                sink_name << "\t" << (sink_orientation_forward ? "+" : "-") <<
                    (t == 2 ? ("\t" + std::to_string(source_begin)) : "") << "\t" << 
                alignment;
            for (auto i : tags){
                st << "\t" << i.second.to_string();
            }
            return st.str();
        }
    };

    struct gap_elem{
        std::string id;
        std::string source_name;
        std::string sink_name;
        int32_t distance;
        std::map<std::string, opt_elem> tags;
        std::string to_string_2(){
            std::stringstream st;
            st << "G" << "\t" << id << "\t" <<
                source_name << "\t" << sink_name << "\t" <<
                distance;
            for (auto t : tags){
                st << "\t" << t.second.to_string();
            }
            return st.str();
        }

    };

    struct fragment_elem{
        std::string id;
        std::string ref;
        bool ref_orientation = true;
        uint32_t seg_begin;
        uint32_t seg_end;
        uint32_t frag_begin;
        uint32_t frag_end;
        std::bitset<4> ends;
        std::string alignment;
        std::map<std::string, opt_elem> tags;
        std::string to_string_2(){
            std::stringstream st;
            st << "F" << "\t" << id << "\t" <<
                ref << (ref_orientation ? "+" : "-") << "\t" <<
                seg_begin << (ends[0] ? "$" : "") << "\t" <<
                seg_end << (ends[1] ? "$" : "") << "\t" <<
                frag_begin << (ends[2] ? "$" : "") << "\t" <<
                frag_end << (ends[3] ? "$" : "") << "\t" <<
                alignment;
                if (tags.size() > 0){
                    for (auto i : tags){
                        st << "\t" << i.second.to_string();
                    }
                }
                return st.str();
        }
    };

    struct group_elem{
        std::string id;
        bool ordered = false;
        std::vector<std::string> items;
        std::vector<bool> orientations;
        std::map<std::string, opt_elem> tags;
        std::string to_string_2(){
            std::stringstream st;
            st << (ordered ? "O" : "U") << "\t" << id << "\t";
            st << items[0] << (ordered ? (orientations[0] ? "+" : "-") : "" );
            for (int i = 1; i < items.size(); ++i){
                st << " " << items[i] << (ordered ? (orientations[i] ? "+" : "-") : "");
            }
            for (auto i : tags){
                st << "\t" << i.second.to_string();
            }
            return st.str();
        }
    };

    int determine_line_type(char* line);


    class GFAKluge{
			friend std::ostream& operator<<(std::ostream& os, GFAKluge& g);

        public:
            GFAKluge();
            ~GFAKluge();
            bool parse_gfa_file(const std::string &filename);
            bool parse_gfa_file(std::istream& gfa_stream);

            //TODO: we should enforce graph structure,
            //i.e.:
            //1. Throw errors on dangling edges
            //2. Guarantee contained elements fall in valid sequences
            //3. Perhaps links and containeds should be added using methods like
            //  add_contained(contained_elem c)
            // to guarantee that they are actually added by their source.
 			void add_link(std::string seq_name, link_elem l);
            void add_link(sequence_elem s, link_elem l);

            /**
             * GFA2.0 handlers
             */
	        void add_edge(std::string seqname, edge_elem e);
            void add_edge(sequence_elem s, edge_elem e);

	        void add_fragment(std::string seqname, fragment_elem f);
            void add_fragment(sequence_elem s, fragment_elem f);

            void add_gap(gap_elem g);

            void add_group(group_elem g);

            void groups_as_paths();

            // Minimize memory consumption by converting everything
            // to GFA2 compatible containers. We can then convert
            // this if we want to go back to GFA1.
            //void compact();




            /** End GFA2.0. Begin 1.0 / 0.1 **/

	        void add_contained(std::string seq_name, contained_elem c);
            void add_contained(sequence_elem s, contained_elem c);

        void add_alignment(std::string s, alignment_elem a);
            void add_alignment(sequence_elem s, alignment_elem a);

            void add_sequence(sequence_elem s);
            void add_path(std::string pathname, path_elem path);
            void add_walk(std::string walkname, walk_elem w);

            double get_version();
            void set_version(double version);
            void set_version();
            void set_walks(bool ws);

	        std::vector<link_elem> get_links(sequence_elem seq);
        std::vector<link_elem> get_links(std::string seq_name);

        std::vector<contained_elem> get_contained(std::string seq_name);
        std::vector<contained_elem> get_contained(sequence_elem seq);

        std::vector<alignment_elem> get_alignments(std::string seq_name);
        std::vector<alignment_elem> get_alignments(sequence_elem seq);

        std::map<std::string, header_elem> get_header();
        std::map<std::string, sequence_elem, custom_key> get_name_to_seq();
        std::map<std::string, std::vector<link_elem> > get_seq_to_link();
        std::map<std::string, std::vector<contained_elem> > get_seq_to_contained();
        std::map<std::string, std::vector<alignment_elem> > get_seq_to_alignment();
        std::map<std::string, std::vector<walk_elem> > get_seq_to_walks();
        std::map<std::string, path_elem> get_name_to_path();

            // GFA2 getters
        std::map<std::string, std::vector<edge_elem>> get_seq_to_edges();
        std::map<std::string, std::vector<fragment_elem>> get_seq_to_fragments();
        std::map<std::string, std::vector<gap_elem>> get_seq_to_gaps();
        std::map<std::string, group_elem> get_groups();

            void paths_as_walks();
            void walks_as_paths();

            void gfa_2_ize();
            void gfa_1_ize();

            // TODO check whether writing to file is functional
            // Perhaps a write_gfa_file(string filename) method too?
            std::string to_string();
            std::string to_string_2();
            std::string block_order_string();
            std::string block_order_string_2();

            // ID manipulators
        std::tuple<uint64_t, uint64_t, uint64_t, uint64_t, uint64_t> max_ids();
            std::string max_ids_string();
            void re_id(std::tuple<uint64_t, uint64_t, uint64_t, uint64_t, uint64_t>& new_mx);
            void re_id(std::string new_mx_str);

            void merge(GFAKluge& gg);

            /** Assembly stats **/
            double get_N50();
            double get_N90();

            int get_L50();
            int get_L90();
            // uint64_t num_contigs();
            // double simple_connectivity() // reports avg edges / sequence
            // double weighted_connectivity() // weight areas of high connectivity higher
            //      behaves kind of like [(N_edges)^2 / (N_seqs)^2 if N_edges > 2]
            
            
        private:
            bool use_walks = false;
            // Store whether we've already gone walks->paths and paths->walks
            bool normalized_paths = false;
            bool normalized_walks = false;
            
            bool one_compat = false;
            bool two_compat = false;

            uint64_t next_set_or_path_id = 0;
            uint64_t base_seq_id = 0;
            uint64_t base_edge_id = 0;
            uint64_t base_gap_id = 0;
            uint64_t base_frag_id = 0;
            uint64_t base_group_id = 0;
            

            double version = 0.0;
        std::map<std::string, header_elem> header;
        std::map<std::string, std::vector<contained_elem> > seq_to_contained;
        std::map<std::string, std::vector<link_elem> > seq_to_link;
            //Since we can't compare sequence elements for hashing,
            // we cheat and use their names (which are only sort of guaranteed to be
            // unique.
        std::map<std::string, sequence_elem, custom_key> name_to_seq;
        std::vector<std::string> split(std::string s, char delim);
        bool string_is_number(std::string s);
        std::string join(std::vector<std::string> splits, std::string glue);
        std::map<std::string, std::vector<alignment_elem> > seq_to_alignment;
        std::string header_string(std::map<std::string, header_elem>& opts);
        std::string opt_string(std::vector<opt_elem> opts);
            
        std::map<std::string, std::vector<walk_elem> > seq_to_walks;
        std::map<std::string, path_elem> name_to_path;

            /** GFA 2.0 containers **/
        std::map<std::string, std::vector<fragment_elem> > seq_to_fragments;
        std::map<std::string, std::vector<edge_elem> > seq_to_edges;
        std::map<std::string, std::vector<gap_elem> > seq_to_gaps;
        std::map<std::string, group_elem> groups;

            double n50 = -1;
            double l50 = -1;
            double n90 = -1;
            double l90 = -1;
            double simple_connectivity = -1;
            double weighted_connectivity = -1;



            

    };

};
#endif
