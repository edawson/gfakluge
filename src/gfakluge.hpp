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
#include <sys/stat.h>

#include "tinyfa.hpp"
#include "pliib.hpp"

namespace gfak{

    // ( <header> | <segment> | <fragment> | <edge> | <gap> | <group> )

    enum gfa_line_types {HEADER_LINE,SEGMENT_LINE,FRAGMENT_LINE,EDGE_LINE,GAP_LINE,GROUP_LINE,PATH_LINE};

    inline int detect_line_type(const string& line){
        char c = line[0];
        if (c == 'S'){
            return SEGMENT_LINE;
        }
        else if (c == 'E' || c == 'L' || c == 'C'){
            return EDGE_LINE;
        }
        else if (c == 'P'){
            return PATH_LINE;
        }
        else if (c == 'F'){
            return FRAGMENT_LINE;
        }
        else if (c == 'G'){
            return GAP_LINE;
        }
        else if (c == 'O' || c == 'U'){
            return GROUP_LINE;
        }
        else{
            cerr << "Unknown GFA line type: " << line << endl;
            return -1;
        }
        
    };
    

    // Provides the proper sorting behavior,
    // where number-based keys get sorted in numerical order
    // and strings are in lexicographic order.
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

    // Contains basic stats about a GFA graph.
    struct gfak_stats_t{
        uint64_t num_nodes = 0;
        uint64_t num_edges = 0;
        uint64_t num_fragments = 0;
        uint64_t num_gaps = 0;
        uint64_t num_paths = 0;
        
        double N50 = 0.0;
        double N90 = 0.0;
    };

    // Contains a GFA 1/2 header line,
    // Which begin with an 'H' and hold an array of <key:type:value> triples.
    struct header_elem {
        std::string key;
        std::string type;
        std::string val;

        std::string to_string() const{
            std::ostringstream st;
            st << 'H' << '\t' << key << ':' << type << ':' << val << '\n';
            return st.str();
        }
    };

    // Encodes a <key:type:value> triple, which are used as
    // annotations for other GFA line types.
    struct opt_elem {
        std::string key;
        std::string type;
        std::string val;
        std::string to_string() const {
            std::stringstream st;
            st << key << ":" << type << ":" << val;
            return st.str();
        }
    };

    /** 
     * Encodes an annotation line.
    *  These begin with an 'x', and are not part of the official
    *  community GFA spec, but have appeared in early GFA files.
    */
    struct annotation_elem{
        std::string key;
        std::string info;
    };

    /**
     *  Encodes a 'P' (path) line from GFA 0.1-1.0.
     *  These describe a trace through sequence_elem elements (i.e. nodes) in a graph
     *  These are analogous and convertible to ordered groups in GFA2.
     */
    struct path_elem{
        std::string name;
        std::vector<std::string> segment_names;
        // Orientation of each segment occurrence in the path. True for forward (+), false for reverse (-)
        std::vector<bool> orientations;
        std::vector<std::string> overlaps;
        std::map<std::string, opt_elem> opt_fields;
        path_elem(){
            std::string name = "";
            std::vector<std::string> segment_names;
            std::vector<bool> orientations;
            std::vector<std::string> overlaps;
            std::map<std::string, opt_elem> opt_fields;
        };

        /**
         *  Adds a GFA 0.1-style path_element (a "Walk") to a
         *  GFA1-style path container (which is a collection of these elements).
         */
        void add_ranked_segment( const int& rank, const string& seg_name, const bool& ori, const string& overlap, vector<opt_elem> opts){
            int corrected_rank = rank;
            if (rank == 0){
                corrected_rank = this->segment_names.size() + 1;
            }
            segment_names.insert( segment_names.begin() + corrected_rank - 1, seg_name);
            orientations.insert( orientations.begin() + corrected_rank - 1, ori);
            overlaps.insert( overlaps.begin() + corrected_rank - 1, overlap);
        };
        /**
         *  Writes a path to a string in GFA1 format.
         */
        std::string to_string_1() const {
            std::ostringstream st;
            std::vector<std::string> p_segs;
            for (int i = 0; i < segment_names.size(); ++i){
                p_segs.push_back(segment_names[i] + (orientations[i] ? "+" : "-") );
            }
            st << 'P' << '\t' << name << '\t' << pliib::join(p_segs, ",") << '\t' << pliib::join(overlaps, ",");
            std::vector<std::string> opt_strs (opt_fields.size());
            int count = 0;
            for (auto op : opt_fields){
                opt_strs[count++] = op.second.to_string(); 
            }
            st << '\t' << pliib::join(opt_strs, "\t");
            return st.str();
        }
        /**
         *  Writes a path to a string in GFA2 format,
         *  which is identical to the output format for an ordered group
         */
        std::string to_string_2() const{
            stringstream st;
            std::vector<std::string> p_segs;
            for (int i = 0; i < segment_names.size(); ++i){
                p_segs.push_back(segment_names[i] + (orientations[i] ? "+" : "-") );
            }
            st << "O" << '\t' << name << '\t' << pliib::join(p_segs, ",");
            return st.str();
        }

        /**
         *  Writes a path as GFA0.1-style walks to an outstream
         */
        void write_as_walks(std::ostream& os){
            std::stringstream st;
            int32_t rank = 0;
            for (int i = 0; i < this->segment_names.size(); ++i){
                ++rank;
                st << 'W' << '\t' << this->segment_names[i] << '\t' << 
                    this->name << '\t' << rank << '\t' <<
                     (this->orientations[i] ? "+" : "-");
                if (this->overlaps.size() == this->segment_names.size()){
                   st << '\t' << overlaps[i];
                }
                st << endl;
                os << st.str();
                st.str("");
            }

        }

    };


    /**
     *  Represents the non-spec alignment line element.
     * These occasionally appear in early GFA 0.1.
     */
    struct alignment_elem{
        std::string source_name;
        int position;
        std::string ref;
        bool source_orientation_forward;
        int length;
        std::map<std::string, std::string> opt_fields;
    };

    /**
     *  Represents a portion of sequence (i.e. a node)
     *  in a graph across all GFA versions.
     */
    struct sequence_elem{
        std::string sequence = "*";
        std::string name = "*";
        uint64_t length = UINT64_MAX;
        std::vector<opt_elem> opt_fields;
        uint32_t id = 0;
        std::string to_string_2() const{
            std::ostringstream st;

            st << "S" << "\t" << name << "\t" << length << "\t" << sequence;
            if (opt_fields.size() > 0){
                for (auto i : opt_fields){
                    st << "\t" << i.to_string();
                }
            }
            return st.str();
        }
        std::string to_string_1() const{
            std::ostringstream st;
            st << "S" << "\t" << name << "\t" << sequence;
            if (opt_fields.size() > 0){
                for (auto i : opt_fields){
                    st << "\t" << i.to_string();
                }
            }
            return st.str();
        }
        /**
         *   Write the sequence_elems name and sequence
         *   in FASTA format.
         */
        std::string as_fasta_record() const{
            std::ostringstream st;
            st << '>' << ' ' << name << endl
                << sequence << endl;
            return st.str();
        };
    };

    struct link_elem{
        std::string source_name;
        std::string sink_name;
        bool source_orientation_forward;
        bool sink_orientation_forward;
        std::string cigar;
        std::map<std::string, opt_elem> opt_fields;
        std::string to_string() const{
            std::ostringstream st;
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


    /**
     *  Represents an edge line in GFA.
     *  Can represent "links", which are end-to-end, or "containments", which
     *  represent overlaps. All GFA versions use this representation - the link_elem
     *  and containment_elem structs get converted to this internally and exist
     *  for convenience.
     */
    struct edge_elem{

        edge_elem(){
            
        }
        edge_elem(const link_elem& l){
           source_name = l.source_name;
           sink_name = l.sink_name;
           type = 1;
           source_orientation_forward = l.source_orientation_forward;
           sink_orientation_forward = l.source_orientation_forward;
           alignment = l.cigar;
           tags = l.opt_fields;
           ends.set(0,1);
           ends.set(1,1);
           ends.set(2,0);
           ends.set(3,0);
           id = "*";
        }
        edge_elem(const contained_elem& c){
           source_name = c.source_name;
           sink_name = c.sink_name;
           type = 2;
           source_orientation_forward = c.source_orientation_forward;
           sink_orientation_forward = c.source_orientation_forward;
           alignment = c.cigar;
           tags = c.opt_fields;
           ends.set(0,0);
           ends.set(1,0);
           ends.set(2,0);
           ends.set(3,0);
           id = "*";
        }
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
        std::string to_string_2() const{
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
            std::string to_string_2() const {
                std::stringstream st;
                st << "G" << "\t" << id << "\t" <<
                    source_name << "\t" << sink_name << "\t" <<
                    distance;
                for (auto t : tags){
                    st << "\t" << t.second.to_string();
                }
                return st.str();
            }
            std::string to_string() const{
                return to_string_2();
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
            std::string to_string_2() const{
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
            std::string to_string() const{
                return to_string_2();
            }
        };

        struct group_elem{
            std::string id;
            bool ordered = false;
            std::vector<std::string> items;
            std::vector<bool> orientations;
            std::map<std::string, opt_elem> tags;
            std::string to_string_1() const{
                std::stringstream st;
                if (!ordered){
                    return "";
                }

                st << "P" << "\t" << id << "\t";
                st << items[0] << (ordered ? (orientations[0] ? "+" : "-") : "" );
                for (size_t i = 1; i < items.size(); ++i){
                    st << " " << items[i] << (ordered ? (orientations[i] ? "+" : "-") : "");
                }
                for (auto i : tags){
                    st << "\t" << i.second.to_string();
                }

                return st.str();
            }
            std::string to_string_2() const{
                std::stringstream st;
                st << (ordered ? "O" : "U") << "\t" << id << "\t";
                st << items[0] << (ordered ? (orientations[0] ? "+" : "-") : "" );
                for (size_t i = 1; i < items.size(); ++i){
                    st << " " << items[i] << (ordered ? (orientations[i] ? "+" : "-") : "");
                }
                for (auto i : tags){
                    st << "\t" << i.second.to_string();
                }
                return st.str();
            }
        std::string to_walk_string(){
            if (ordered){
                return "";
            }
            int rank = 0;
            stringstream st;
            for (size_t i = 0; i < items.size(); ++i){
                st << "W" << items[i] << '\t' << ++rank << '\t' << orientations[i] << "*" << endl;
            }
            return st.str();
        }
        };

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
            
            /**
             * GFA2.0 handlers
             */
            void add_sequence(sequence_elem s);
            void add_edge(const std::string& seqname, const edge_elem& e);
            void add_edge(const sequence_elem& s, const edge_elem& e);

            void add_fragment(const std::string& seqname, const fragment_elem& f);
            void add_fragment(const sequence_elem& s, const fragment_elem& f);

            void add_gap(const gap_elem& g);

            void add_group(const group_elem& g);

            void groups_as_paths();

            /** End GFA2.0. Begin 1.0 / 0.1 **/

            void add_contained(std::string seq_name, contained_elem c);
            void add_contained(sequence_elem s, contained_elem c);

            void add_alignment(std::string s, alignment_elem a);
            void add_alignment(sequence_elem s, alignment_elem a);


            void add_path(std::string pathname, path_elem path);
            void add_walk(std::string pathname, const int& rank, const string& segname, const bool& ori, const string& overlap, vector<opt_elem> opts);
            void add_link(const std::string& seq_name, const link_elem& l);
            void add_link(const sequence_elem& s, const link_elem& l);


            /** Versioning functions **/
            double get_version();
            void set_version(double version);
            void set_version();
            void set_walks(bool ws);

            /** Getter methods for elements, to keep users out of our data structures **/

            std::vector<contained_elem> get_contained(std::string seq_name);
            std::vector<contained_elem> get_contained(sequence_elem seq);

            std::vector<alignment_elem> get_alignments(std::string seq_name);
            std::vector<alignment_elem> get_alignments(sequence_elem seq);

            std::map<std::string, header_elem> get_header();
            std::map<std::string, sequence_elem, custom_key> get_name_to_seq();
            std::map<std::string, std::vector<link_elem> > get_seq_to_link();
            std::map<std::string, std::vector<contained_elem> > get_seq_to_contained();
            std::map<std::string, std::vector<alignment_elem> > get_seq_to_alignment();
            std::map<std::string, path_elem> get_name_to_path();

            // GFA2 getters
            std::map<std::string, std::vector<edge_elem>> get_seq_to_edges();
            std::map<std::string, std::vector<fragment_elem>> get_seq_to_fragments();
            std::map<std::string, std::vector<gap_elem>> get_seq_to_gaps();
            std::map<std::string, group_elem> get_groups();

            /** Convert paths to walks and walks to paths **/
            void gfa_2_ize();
            void gfa_1_ize();
            void compatibilize();

            // TODO check whether writing to file is functional
            // Perhaps a write_gfa_file(string filename) method too?
            std::string to_string();
            std::string to_string_2();
            std::string block_order_string();
            std::string block_order_string_2();
            void output_to_stream(std::ostream& os, bool output_block_order = false);

            /** Methods for folks that want streaming output. **/
            void write_element(std::ostream& os, const sequence_elem& s);
            void write_element(std::ostream& os, edge_elem e);
            void write_element(std::ostream& os, const fragment_elem& f);
            void write_element(std::ostream& os, const group_elem& g);
            void write_element(std::ostream& os, const gap_elem& g);
            void write_element(std::ostream& os, const header_elem& h);

            std::string header_string();
            
            // ID manipulators
            std::tuple<uint64_t, uint64_t, uint64_t, uint64_t, uint64_t> max_ids();
            std::string max_ids_string();
            void re_id(std::tuple<uint64_t, uint64_t, uint64_t, uint64_t, uint64_t>& new_mx);
            void re_id(std::string new_mx_str);

            /** Merge two GFA graphs **/
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
            
            void fill_sequences(const char* fasta_file);

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
            bool string_is_number(std::string s);
            std::string join(const std::vector<std::string>& splits, const std::string& glue);
            std::map<std::string, std::vector<alignment_elem> > seq_to_alignment;
            std::string header_string(std::map<std::string, header_elem>& opts);
            std::string opt_string(std::vector<opt_elem> opts);

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
