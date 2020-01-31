#ifndef GFAK_HPP
#define GFAK_HPP

#include <string>
#include <sstream>
#include <istream>
#include <algorithm>
#include <functional>
#include <map>
#include <vector>
#include <sstream>
#include <cstdlib>
#include <iostream>
#include <bitset>
#include <unordered_set>
#include <sys/stat.h>
#include <cstdio>
#include <sys/types.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
#include <cassert>


#include "tinyFA.hpp"
#include "pliib.hpp"

namespace gfak{

    // ( <header> | <segment> | <fragment> | <edge> | <gap> | <group> )

    enum gfa_line_types {HEADER_LINE,SEGMENT_LINE,FRAGMENT_LINE,EDGE_LINE,GAP_LINE,GROUP_LINE,PATH_LINE,LINK_LINE,CONTAINED_LINE,WALK_LINE};

    static inline int determine_line_type(const char* line){
        if (line[0] == 'H'){
            return HEADER_LINE;
        }
        else if (line[0] == 'S'){
            return SEGMENT_LINE;
        }
        else if (line[0] == 'E'){
            return EDGE_LINE;
        }
        else if (line[0] == 'L'){
            return LINK_LINE;
        }
        else if (line[0] == 'C'){
            return CONTAINED_LINE;
        }
        else if (line[0] == 'U' || line[0] == 'O'){
            return GROUP_LINE;
        }
        else if (line[0] == 'P'){
            return PATH_LINE;
        }
        else if (line[0] == 'W'){
            return WALK_LINE;
        }
        else if (line[0] == 'G'){
            return GAP_LINE;
        }
        else if (line[0] == 'F'){
            return FRAGMENT_LINE;
        }
        else{
            return -1;
        }
    }


    // Provides the proper sorting behavior,
    // where number-based keys get sorted in numerical order
    // and std::strings are in lexicographic order.
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
        std::uint64_t num_nodes = 0;
        std::uint64_t num_edges = 0;
        std::uint64_t num_fragments = 0;
        std::uint64_t num_gaps = 0;
        std::uint64_t num_paths = 0;

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
        }

        /**
         *  Adds a GFA 0.1-style path_element (a "Walk") to a
         *  GFA1-style path container (which is a collection of these elements).
         */
        void add_ranked_segment( const int& rank, const std::string& seg_name, const bool& ori, const std::string& overlap, std::vector<opt_elem> opts){
            size_t corrected_rank = rank;
            if (rank == 0){
                corrected_rank = this->segment_names.size() + 1;
            }
            segment_names.insert( segment_names.begin() + corrected_rank - 1, seg_name);
            orientations.insert( orientations.begin() + corrected_rank - 1, ori);
            overlaps.insert( overlaps.begin() + corrected_rank - 1, overlap);
        }
        /**
         *  Writes a path to a std::string in GFA1 format.
         */
        std::string to_string_1() const {
            std::ostringstream st;
            std::vector<std::string> p_segs;
            for (size_t i = 0; i < segment_names.size(); ++i){
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
         *  Writes a path to a std::string in GFA2 format,
         *  which is identical to the output format for an ordered group
         */
        std::string to_string_2() const{
            std::ostringstream st;
            std::vector<std::string> p_segs;
            for (size_t i = 0; i < segment_names.size(); ++i){
                p_segs.push_back(segment_names[i] + (orientations[i] ? "+" : "-") );
            }
            st << "O" << '\t' << name << '\t' << pliib::join(p_segs, ",");
            return st.str();
        }

        /**
         *  Writes a path as GFA0.1-style walks to an outstream
         */
        void write_as_walks(std::ostream& os){
            std::ostringstream st;
            int32_t rank = 0;
            for (size_t i = 0; i < this->segment_names.size(); ++i){
                ++rank;
                st << 'W' << '\t' << this->segment_names[i] << '\t' << 
                    this->name << '\t' << rank << '\t' <<
                    (this->orientations[i] ? "+" : "-");
                if (this->overlaps.size() == this->segment_names.size()){
                    st << '\t' << overlaps[i];
                }
                st << std::endl;
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
        uint64_t id = 0;
        std::string sequence = "*";
        std::string name = "*";
        uint64_t length = UINT64_MAX;
        std::vector<opt_elem> opt_fields;
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
            st << '>' << ' ' << name << std::endl
                << sequence;
            return st.str();
        }
    };

    /** Represents a link (a non-contained edge) between two sequence_elems
     *  Note: these are no longer stored internally; instead, they get converted to
     *  edge_elems.
     */
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

    /**
     * Represents a containment edge, which defines how one sequence element is
     * contained within another. Like link_elems, these are no longer stored internally,
     * but are converted to edge_elems.
     */
    struct contained_elem {
        std::string source_name;
        std::string sink_name;
        bool source_orientation_forward;
        bool sink_orientation_forward;
        std::uint64_t pos;
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
        // 0: unset, 1: link, 2: containment, 3: generic edge or both C/L set (which should be impossible)
        int type;
        std::string source_name;
        std::string sink_name;
        bool source_orientation_forward;
        bool sink_orientation_forward;
        std::bitset<4> ends;
        std::uint64_t source_begin = 0;
        std::uint64_t source_end = 0;
        std::uint64_t sink_begin = 0;
        std::uint64_t sink_end = 0;
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
            std::ostringstream st;
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
            std::ostringstream st;
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

        /** Represents a GAP type line ('G'), which is specific to GFA2 */
        struct gap_elem{
            std::string id;
            std::string source_name;
            std::string sink_name;
            std::int32_t distance;
            std::map<std::string, opt_elem> tags;
            std::string to_string_2() const {
                std::ostringstream st;
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

        /** Represents a fragment, which is a portion of sequence contained within
         * a sequence. These are GFA2-specific and are stored relative to their sequence_elem
         */
        struct fragment_elem{
            std::string id;
            std::string ref;
            bool ref_orientation = true;
            std::uint32_t seg_begin;
            std::uint32_t seg_end;
            std::uint32_t frag_begin;
            std::uint32_t frag_end;
            std::bitset<4> ends;
            std::string alignment;
            std::map<std::string, opt_elem> tags;
            std::string to_string_2() const{
                std::ostringstream st;
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

        /** Represents a group of sequence_elems, which could be an ordered path
         *  or an unordered collection.
         */
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
                std::ostringstream st;
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
            /** Convert a group (which must be ordered) to a GFA0.1 style walk std::string. */
            std::string to_walk_string(){
                if (!ordered){
                    return "";
                }
                int rank = 0;
                std::ostringstream st;
                for (size_t i = 0; i < items.size(); ++i){
                    st << "W" << items[i] << '\t' << ++rank << '\t' << orientations[i] << "*" << std::endl;
                }
                return st.str();
            }
        };

        inline size_t mmap_open(const std::string& filename, char*& buf, int& fd) {
            fd = -1;
            assert(!filename.empty());
            // open in binary mode as we are reading from this interface
            fd = open(filename.c_str(), O_RDWR);
            if (fd == -1) {
                assert(false);
            }
            struct stat stats;
            if (-1 == fstat(fd, &stats)) {
                assert(false);
            }
            size_t fsize = stats.st_size;
            if (!(buf =
                        (char*) mmap(NULL,
                            fsize,
                            PROT_READ | PROT_WRITE,
                            MAP_SHARED,
                            fd,
                            0))) {
                assert(false);
            }
            madvise((void*)buf, fsize, POSIX_MADV_WILLNEED | POSIX_MADV_SEQUENTIAL);
            return fsize;
        }

        inline void mmap_close(char*& buf, int& fd, size_t fsize) {
            if (buf) {
                munmap(buf, fsize);
                buf = 0;
            }
            if (fd) {
                close(fd);
                fd = 0;
            }
        }

        class GFAKluge{
            inline friend std::ostream& operator<<(std::ostream& os, GFAKluge& g){
                g.gfa_1_ize();
                g.gfa_2_ize();

                //os << g.to_string();
                g.output_to_stream(os);
                return os;
            }

            private:
            bool use_walks = false;
            // Store whether we've already gone walks->paths and paths->walks
            bool normalized_paths = false;
            bool normalized_walks = false;

            bool one_compat = false;
            bool two_compat = false;

            std::uint64_t next_set_or_path_id = 0;
            std::uint64_t base_seq_id = 0;
            std::uint64_t base_edge_id = 0;
            std::uint64_t base_gap_id = 0;
            std::uint64_t base_frag_id = 0;
            std::uint64_t base_group_id = 0;


            double version = 0.0;
            std::map<std::string, header_elem> header;
            std::map<std::string, std::vector<contained_elem> > seq_to_contained;
            std::map<std::string, std::vector<link_elem> > seq_to_link;
            std::map<std::string, std::vector<alignment_elem> > seq_to_alignment;
            std::map<std::string, path_elem> name_to_path;
            //Since we can't compare sequence elements for hashing,
            // we cheat and use their names (which are only sort of guaranteed to be
            // unique.
            std::map<std::string, sequence_elem, custom_key> name_to_seq;

            inline bool string_is_number(std::string s){
                bool ret = true;
                std::string::iterator it;
                for (it = s.begin(); it != s.end(); it++){
                    if (!isdigit(*it)){
                        ret = false;
                        return ret;
                    }
                }
                return !(s.empty());
            }

            inline std::string header_string(std::map<std::string, header_elem>& headers){
                std::string ret = "H";
                std::map<std::string, header_elem>::iterator it;
                for (it = headers.begin(); it != headers.end(); it++){
                    ret += "\t";
                    header_elem h = it->second;
                    std::string t[] = {h.key, h.type, h.val};
                    std::vector<std::string> temp = std::vector<std::string> (t, t + sizeof(t) / sizeof(std::string));
                    ret += pliib::join(temp, ":");
                }
                return ret;

            }
            inline std::string opt_string(std::vector<opt_elem> opts){
                std::string ret = "";
                for (size_t i = 0; i < opts.size(); i++){
                    opt_elem o = opts[i];
                    if (i != 0){
                        ret += "\t";
                    }
                    std::string t [] = {o.key, o.type, o.val};
                    std::vector<std::string> temp = std::vector<std::string> (t, t + sizeof(t) / sizeof(std::string));
                    ret += pliib::join(temp, ":");
                }
                return ret;
            }



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

            public:
            GFAKluge(){
                std::map<std::string, header_elem> header;
                std::map<std::string, std::vector<contained_elem>, custom_key > seq_to_contained;
                std::map<std::string, std::vector<link_elem>, custom_key > seq_to_link;
                std::map<std::string, std::vector<alignment_elem> , custom_key> seq_to_alignment;
                //Since we can't compare sequence elements for hashing,
                // we cheat and use their names (which are sort of guaranteed to be
                // unique.
                std::map<std::string, sequence_elem, custom_key> name_to_seq;
                std::map<std::string, path_elem> name_to_path;
                //cout << fixed << setprecision(2);


                /** GFA 2.0 containers **/
                std::map<std::string, std::vector<fragment_elem> > seq_to_fragments;
                std::map<std::string, std::vector<edge_elem> > seq_to_edges;
                std::map<std::string, std::vector<gap_elem> > seq_to_gaps;
                std::map<std::string, group_elem> groups;
            }
            ~GFAKluge(){

            }

            inline double detect_version_from_file(const char* filename){
                std::ifstream gfi;
                gfi.open(filename, std::ifstream::in);
                if (!gfi.good()){
                    std::cerr << "Couldn't open GFA file " << filename << "." << std::endl;
                    exit(1);
                }
                std::string line;
                while (getline(gfi, line)){
                    if (determine_line_type(line.c_str()) == HEADER_LINE){
                        std::vector<std::string> tokens = pliib::split(line, '\t');
                        header_elem h;
                        std::vector<std::string> line_tokens = pliib::split(tokens[1], ':');
                        //TODO this is not well implemented
                        // GFA places no guarantees on header format
                        h.key = line_tokens[0];
                        h.type = line_tokens[1];
                        h.val = line_tokens[2];
                        if (h.key.compare("VN") == 0){
                            set_version(stod(h.val));
                            break;
                        }
                        header[h.key] = h;
                    } 
                }
            }

            inline void for_each_sequence_line_in_file(const char* filename, std::function<void(gfak::sequence_elem)> func){
	        std::ifstream gfi;
                gfi.open(filename, std::ifstream::in);
                if (!gfi.good()){
                    std::cerr << "Couldn't open GFA file " << filename << "." << std::endl;
                    exit(1);
                }
		for_each_sequence_line_in_file(gfi, func);
 	
	    };

            inline void for_each_sequence_line_in_file(std::istream& gfi, std::function<void(gfak::sequence_elem)> func){
                std::string line;
                while (getline(gfi, line)){
                    if (determine_line_type(line.c_str()) == SEGMENT_LINE){
                        std::vector<std::string> tokens = pliib::split(line, '\t');
                        sequence_elem s;
                        int tag_index = 3;
                        s.name = tokens[1];

                        if (this->version >= 2.0 || string_is_number(tokens[2])){
                            s.length = stoi(tokens[2]); 
                            s.sequence = tokens[3];
                            tag_index = 4;
                        }
                        else{

                            s.sequence = tokens[2];
                            if (tokens[2] == "*"){
                                s.length = UINT64_MAX;
                            }
                            else{
                                s.length = s.sequence.length();
                            }
                            tag_index = 3;
                        }

                        size_t i;
                        if (tokens.size() > 3){
                            for (i = tag_index; i < tokens.size(); i++){
                                //opt fields are in key:type:val format
                                std::vector<std::string> opt_field = pliib::split(tokens[i], ':');
                                opt_elem o;
                                o.key = opt_field[0];
                                o.type = opt_field[1];
                                o.val = pliib::join(std::vector<std::string> (opt_field.begin() + 2, opt_field.end()), ":");
                                s.opt_fields.push_back(o);
                                if (o.key == "LN" && s.length == UINT64_MAX){
                                    s.length = stoul(o.val);
                                }
                            }
                        }
                        func(s);
                    }
                    else if (determine_line_type(line.c_str()) == HEADER_LINE){
                        std::vector<std::string> tokens = pliib::split(line, '\t');
                        header_elem h;
                        std::vector<std::string> line_tokens = pliib::split(tokens[1], ':');
                        //TODO this is not well implemented
                        // GFA places no guarantees on header format
                        h.key = line_tokens[0];
                        h.type = line_tokens[1];
                        h.val = line_tokens[2];
                        if (h.key.compare("VN") == 0){
                            set_version(stod(h.val));
                        }
                        header[h.key] = h;
                    } 
                }
            }


            inline void for_each_edge_line_in_file(char* filename, std::function<void(gfak::edge_elem)> func){
                std::ifstream gfi;
                gfi.open(filename, std::ifstream::in);
                if (!gfi.good()){
                    std::cerr << "Couldn't open GFA file " << filename << "." << std::endl;
                    exit(1);
                }

		for_each_edge_line_in_file(gfi, func);
	    };

            inline void for_each_edge_line_in_file(std::istream& gfi, std::function<void(gfak::edge_elem)> func){

                std::string line;
                while (getline(gfi, line)){
                    if (determine_line_type(line.c_str()) == EDGE_LINE){
                        std::vector<std::string> tokens = pliib::split(line, '\t');
                        edge_elem e;
                        e.id = tokens[1];

                        std::string x = tokens[2];
                        e.source_name = x.substr(0, x.length() - 1);
                        e.source_orientation_forward = (x.back() == '+');

                        x = tokens[3];
                        e.sink_name = x.substr(0, x.length() - 1);
                        e.sink_orientation_forward = (x.back() == '+');

                        x = tokens[4];
                        e.ends.set(0, (x.back() == '$' ? 1 : 0));
                        e.source_begin = (e.ends.test(0) ? stoul(x.substr(0, x.length() - 1)) : stoul(x));

                        x = tokens[5];
                        e.ends.set(1, (x.back() == '$' ? 1 : 0));
                        e.source_end = (e.ends.test(1) ? stoul(x.substr(0, x.length() - 1)) : stoul(x));

                        x = tokens[6];
                        e.ends.set(2, (x.back() == '$' ? 1 : 0));
                        e.sink_begin = (e.ends.test(2) ? stoul(x.substr(0, x.length() - 1)) : stoul(x));


                        x = tokens[7];
                        e.ends.set(3, (x.back() == '$' ? 1 : 0));
                        e.sink_end = (e.ends.test(3) ? stoul(x.substr(0, x.length() - 1)) : stoul(x));

                        e.alignment = tokens[8];

                        if (tokens.size() > 9){
                            for (size_t i = 9; i < tokens.size(); i++){
                                //opt fields are in key:type:val format
                                std::vector<std::string> opt_field = pliib::split(tokens[i], ':');
                                opt_elem o;
                                o.key = opt_field[0];
                                o.type = opt_field[1];
                                o.val = pliib::join(std::vector<std::string> (opt_field.begin() + 2, opt_field.end()), ":");
                                e.tags[o.key] = o;

                            }
                        }

                        func(e);
                    }
                    else if (determine_line_type(line.c_str()) == LINK_LINE){
                        std::vector<std::string> tokens = pliib::split(line, '\t');
                        edge_elem e;
                        e.type = 1;
                        e.source_name = tokens[1];
                        e.sink_name = tokens[3];
                        //TODO: search the input std::strings for "-" and "+" and set using ternary operator
                        e.source_orientation_forward = tokens[2] == "+" ? true : false;
                        e.ends.set(0, 1);
                        e.ends.set(1,1);
                        e.ends.set(2,0);
                        e.ends.set(3, 0);
                        e.sink_orientation_forward = tokens[4] == "+" ? true : false;
                        if (tokens.size() >= 6){
                            e.alignment = tokens[5];
                        }
                        else{
                            e.alignment = "*";
                        }

                        if (tokens.size() >= 7){
                            for (size_t i = 6; i < tokens.size(); i++){
                                //opt fields are in key:type:val format
                                std::vector<std::string> opt_field = pliib::split(tokens[i], ':');
                                opt_elem o;
                                o.key = opt_field[0];
                                o.type = opt_field[1];
                                o.val = pliib::join(std::vector<std::string> (opt_field.begin() + 2, opt_field.end()), ":");
                                e.tags[o.key] = o;

                            }
                        }

                        func(e);
                    }
                    else if (determine_line_type(line.c_str()) == CONTAINED_LINE){
                        std::vector<std::string> tokens = pliib::split(line, '\t');
                        contained_elem c;

                        edge_elem e(c);
                        func(e);
                    }
                    else if (determine_line_type(line.c_str()) == HEADER_LINE){
                        std::vector<std::string> tokens = pliib::split(line, '\t');
                        header_elem h;
                        std::vector<std::string> line_tokens = pliib::split(tokens[1], ':');
                        //TODO this is not well implemented
                        // GFA places no guarantees on header format
                        h.key = line_tokens[0];
                        h.type = line_tokens[1];
                        h.val = line_tokens[2];
                        if (h.key.compare("VN") == 0){
                            set_version(stod(h.val));
                        }
                        header[h.key] = h;
                    } 
                }
            }



            // Per-element parsing of paths, only supports GFA 1.0
            //inline void for_each_path_element_in_file(const char* filename, std::function<void(const std::string&, const std::string&, bool, const std::string&)> func){
	inline void for_each_path_element_in_file(const char* filename, std::function<void(const std::string&, const std::string&, bool, const std::string&, bool, bool)> func){
                int gfa_fd = -1;
                char* gfa_buf = nullptr;
                size_t gfa_filesize = mmap_open(filename, gfa_buf, gfa_fd);
                if (gfa_fd == -1) {
                    std::cerr << "Couldn't open GFA file " << filename << "." << std::endl;
                    exit(1);
                }
                std::string line;
                size_t i = 0;
                //bool seen_newline = true;
                while (i < gfa_filesize) {
                    std::string path_name;
                    // scan forward
                    if (i > 0 && gfa_buf[i-1] == '\n' && gfa_buf[i] == 'P') {
                        // path line
                        // scan forward to find name
                        i += 2;
                        while (gfa_buf[i] != '\t') {
                            path_name.push_back(gfa_buf[i++]);
                        }
                        ++i; // get to path id/orientation description
                        size_t j = i;
                        while (gfa_buf[++j] != '\t');
                        ++j; // skip over delimiter
                        // now j points to the overlaps
                        while (gfa_buf[i] != '\t' && gfa_buf[j] != '\n' && j+1 != gfa_filesize) {
                            std::string id;
                            char c = gfa_buf[i];
                            while (c != ',' && c != '\t' && c != '+' && c != '-') {
                                id.push_back(c);
                                c = gfa_buf[++i];
                            }
                            bool is_rev = gfa_buf[i++]=='-';
                            ++i; // skip over delimiter
                            c = gfa_buf[j];
                            std::string overlap;
                            while (c != ',' && c != '\t' && c != '\n') {
                                overlap.push_back(c);
                                if (j+1 == gfa_filesize) break;
                                c = gfa_buf[++j];
                            }
                            ++j; // skip over delimiter
                            func(path_name, id, is_rev, overlap, false, false);
                        }
                    }
                    ++i;
                }
                mmap_close(gfa_buf, gfa_fd, gfa_filesize);
            }

	   /** 
            inline void for_each_path_line_in_file(const char* filename, std::function<void(gfak::path_elem)> func){
                int gfa_fd = -1;
                char* gfa_buf = nullptr;
                size_t gfa_filesize = mmap_open(filename, gfa_buf, gfa_fd);
                if (gfa_fd == -1) {
                    std::cerr << "Couldn't open GFA file " << filename << "." << std::endl;
                    exit(1);
                }
                std::string line;
                size_t i = 0;
                bool seen_newline = true;
                while (i < gfa_filesize) {
                    std::string path_name;
                    // scan forward
                    if (i > 0 && gfa_buf[i-1] == '\n' && gfa_buf[i] == 'P') {
                        // path line
                        // scan forward to find name
                        i += 2;
                        while (gfa_buf[i] != '\t') {
                            path_name.push_back(gfa_buf[i++]);
                        }
                        ++i; // get to path id/orientation description
                        size_t j = i;
                        while (gfa_buf[j++] != '\t');
                        // now j points to the overlaps
 			char b = gfa_buf[i], c = gfa_buf[j];
                        bool is_empty = true;
                        while (b != '\t' && j != gfa_filesize) {
                            string id;
                            if (b == ',') b = gfa_buf[++i];
                            while (b != ',' && b != '\t' && b != '+' && b != '-') {
                                id.push_back(b);
                                b = gfa_buf[++i];
                            }
                            bool is_rev = b=='-';
                            b = gfa_buf[++i];
                            string overlap;
                            if (c == ',') c = gfa_buf[++j];
                            while (c != ',' && c != '\t' && c != '\n' && c != ' ' && j != gfa_filesize) {
				        overlap.push_back(c);
                                c = gfa_buf[++j];
                            }
                            func(path_name, id, is_rev, overlap, false, false);
                            is_empty = false;
                        }
                        if (is_empty) {
                            func(path_name, 0, false, "", true, false);
                        }
                    }
                    ++i;
                }
                mmap_close(gfa_buf, gfa_fd, gfa_filesize);
            };
**/
	    
            inline void for_each_path_line_in_file(const char* filename, std::function<void(gfak::path_elem)> func){
	        std::ifstream gfi;
                gfi.open(filename, std::ifstream::in);
                if (!gfi.good()){
                    std::cerr << "Couldn't open GFA file " << filename << "." << std::endl;
                    exit(1);
                }
		for_each_path_line_in_file(gfi, func);	
	    };

            // Only supports GFA 1.0 style paths
            inline void for_each_path_line_in_file(std::istream& gfi, std::function<void(gfak::path_elem)> func){

                std::string line;
                while (getline(gfi, line)){
                    if (determine_line_type(line.c_str()) == PATH_LINE){
                        std::vector<std::string> tokens = pliib::split(line, '\t');
                        path_elem p;
                        // Parse a GFA 1.0 path element
                        p.name = tokens[1];
                        std::vector<std::string> ids_and_orientations;
                        pliib::split(tokens[2], ',', ids_and_orientations);
                        p.segment_names.resize(ids_and_orientations.size());
                        p.orientations.resize(ids_and_orientations.size());
                        for (size_t t = 0; t < ids_and_orientations.size(); ++t){
                            std::string x = ids_and_orientations[t];
                            bool orientation = ((x.back()) == '+' || x.front() == '+');
                            std::string id = x.substr(0, x.length() - 1);
                            p.segment_names[t] = id;
                            p.orientations[t] = orientation;
                        }
                        if (tokens.size() > 3){
                            std::vector<std::string> spltz = pliib::split(tokens[3], ',');
                            for (auto z : spltz){
                                p.overlaps.push_back(z);
                            }
                            //p.overlaps.assign( spltz.begin(), spltz.end());
                        }
                        else{
                            //p.overlaps.assign(p.segment_names.size(), "*");
                            for (auto z : p.segment_names){
                                p.overlaps.push_back("*");
                            }
                        }
                        func(p);
                    }
                    else if (determine_line_type(line.c_str()) == HEADER_LINE){
                        std::vector<std::string> tokens = pliib::split(line, '\t');
                        header_elem h;
                        std::vector<std::string> line_tokens = pliib::split(tokens[1], ':');
                        //TODO this is not well implemented
                        // GFA places no guarantees on header format
                        h.key = line_tokens[0];
                        h.type = line_tokens[1];
                        h.val = line_tokens[2];
                        if (h.key.compare("VN") == 0){
                            set_version(stod(h.val));
                        }
                        header[h.key] = h;
                    } 
                }
            }

            // Only supports GFA 2.0 style paths (i.e. groups, both ordered and unordered)
            inline void for_each_ordered_group_line_in_file(const char* filename, std::function<void(gfak::group_elem)> func){
                std::ifstream gfi;
                gfi.open(filename, std::ifstream::in);
                if (!gfi.good()){
                    std::cerr << "Couldn't open GFA file " << filename << "." << std::endl;
                    exit(1);
                }
                std::string line;
                while (getline(gfi, line)){
                    if (determine_line_type(line.c_str()) == PATH_LINE){
                        std::vector<std::string> tokens = pliib::split(line, '\t');
                        group_elem g;
                        g.ordered = false;
                        g.id = tokens[0];
                        if (g.id == "*"){
                            g.id = std::to_string(++base_group_id);
                        }
                        g.items = pliib::split(tokens[2], ' ');
                        if (tokens.size() > 8){
                            for (size_t i = 9; i < tokens.size(); i++){
                                //opt fields are in key:type:val format
                                std::vector<std::string> opt_field = pliib::split(tokens[i], ':');
                                opt_elem o;
                                o.key = opt_field[0];
                                o.type = opt_field[1];
                                o.val = pliib::join(std::vector<std::string> (opt_field.begin() + 2, opt_field.end()), ":");
                                g.tags[o.key] = o;

                            }
                        }
                        func(g);
                    }
                    else if (determine_line_type(line.c_str()) == HEADER_LINE){
                        std::vector<std::string> tokens = pliib::split(line, '\t');
                        header_elem h;
                        std::vector<std::string> line_tokens = pliib::split(tokens[1], ':');
                        //TODO this is not well implemented
                        // GFA places no guarantees on header format
                        h.key = line_tokens[0];
                        h.type = line_tokens[1];
                        h.val = line_tokens[2];
                        if (h.key.compare("VN") == 0){
                            set_version(stod(h.val));
                        }
                        header[h.key] = h;
                    } 
                }
            }


            // GFAKluge does not, but maybe should, enforce graph structure,
            // i.e.:
            // 1. Throw errors on dangling edges
            // 2. Guarantee contained elements fall in valid sequences


            /**
             * GFA2.0 handlers
             */
            inline void add_sequence(sequence_elem s){
                name_to_seq[s.name] = s;
            }
            inline void add_edge(const std::string& seqname, const edge_elem& e){
                seq_to_edges[seqname].push_back(e);
            }
            inline void add_edge(const sequence_elem& s, const edge_elem& e){
                add_edge(s.name, e);
            }

            inline void add_fragment(const std::string& seqname, const fragment_elem& f){
                seq_to_fragments[seqname].push_back(f);
            }

            inline void add_fragment(const sequence_elem& s, const fragment_elem& f){
                add_fragment(s.name, f);
            }

            inline void add_gap(const gap_elem& g){
                seq_to_gaps[g.source_name].push_back(g);
            }

            inline void add_group(const group_elem& g){
                this->groups[g.id] = g;
            }

            // Convert group_elems (i.e. U or O lines) to path_elems (P lines)
            inline void groups_as_paths(){
                for (auto g : groups){
                    if (g.second.ordered){
                        path_elem p;
                        p.name = g.first;
                        p.segment_names = g.second.items;
                        if (this->name_to_seq.size() > 0){
                            p.overlaps.resize(p.segment_names.size());
                            for (size_t i = 0; i < p.segment_names.size(); ++i){
                                int len = 0;
                                if (name_to_seq.find(p.segment_names[i]) != name_to_seq.end()){
                                    std::string s(name_to_seq.at(p.segment_names[i]).sequence);
                                    len = s.length();
                                }
                                p.overlaps[i].assign(std::to_string(len) + "M");
                            }
                        }
                        p.orientations = g.second.orientations;
                        p.opt_fields = g.second.tags;
                        add_path(p.name, p);
                    }
                    else{
                        std::cerr << "Group " << g.first << " is unordered; skipping adding it to the paths." << std::endl;
                    }
                }
            }

            /** End GFA2.0. Begin 1.0 / 0.1 **/

            /** Add a containment line to the GFAKluge object
             *  Behind the scenes, an edge_elem is created and
             *  used to represent the contained_elem, which
             *  are provided as syntactic sugar.
             *  N.B.: these are stored relative to the sequence_elem
             *  in which they are contained.
             */
            inline void add_contained(std::string seq_name, contained_elem c){
                seq_to_contained[seq_name].push_back(c);
            }
            inline void add_contained(sequence_elem s, contained_elem c){
                edge_elem e(c);
                seq_to_edges[s.name].push_back(e);
            }

            /**
             * Add an alignment_elem to the GFAKluge object.
             * These are stored by the sequence_elem to which they align.
             */
            inline void add_alignment(std::string s, alignment_elem a){
                seq_to_alignment[s].push_back(a);
            }
            inline void add_alignment(sequence_elem s, alignment_elem a){
                seq_to_alignment[s.name].push_back(a);
            }


            /**
             * Functions for adding paths or walks (which are single elements in an ordered path)
             */
            inline void add_path(std::string pathname, path_elem path){
                name_to_path[pathname] = path;
            }
            inline void add_walk(std::string pathname, const int& rank, const std::string& segname, const bool& ori, const std::string& overlap, std::vector<opt_elem> opts){
                if (name_to_path.find(pathname) == name_to_path.end()){
                    path_elem p;
                    p.name = pathname;
                    add_path(p.name, p);
                }
                name_to_path.at(pathname).add_ranked_segment( rank, segname, ori, overlap, opts);

            }

            /**
             *  Add a link_elem to represent a link between two sequence_elems.
             *  These are stored relative to their source sequence_elem.
             *  N.B.: An edge_elem is created internally to represent the link
             *  and no link_elem is stored. They're simply provided as syntactic sugar.
             */
            inline void add_link(const std::string& seq_name, const link_elem& link){
                edge_elem e(link);
                seq_to_edges[seq_name].push_back(e); 
            }
            inline void add_link(const sequence_elem& s, const link_elem& link){
                edge_elem e(link);
                seq_to_edges[s.name].push_back(e);
            }


            /** Versioning functions **/
            inline double get_version(){
                return this->version;
            }
            inline void set_version(double version){
                header_elem verz;
                verz.key = "VN";
                verz.type="Z";
                this->version = version;
                verz.val = std::to_string((double) this->version).substr(0,3);
                this->header[verz.key] = verz;
                //gfa_1_ize();
                //gfa_2_ize();
            }
            inline void set_version(){
                header_elem verz;
                verz.key = "VN";
                verz.type="Z";
                verz.val = std::to_string((double) this->version).substr(0,3);

                this->header[verz.key] = verz;
            }
            // Use walks, rather than paths, for outputting GFA below v2.0.
            inline void set_walks(bool ws){
                this->use_walks = ws;
            }

            /** Methods for folks that want streaming output.
             *  Writes a single element to an ostream, using either
             *  that element's to_string_1() or to_string_2() function
             *  depending on the version set in the GFAKluge object.
             */
            inline void write_element(std::ostream& os, const sequence_elem& s) {
                if (this->version >= 2.0){
                    os << s.to_string_2();
                }
                else{
                    os << s.to_string_1();
                }
            }
            inline void write_element(std::ostream& os, edge_elem e){
                if (this->version >= 2.0){
                    os << e.to_string_2();
                }
                else{
                    os << e.to_string_1();
                }
            }
            inline void write_element(std::ostream& os, const fragment_elem& f){
                os << f.to_string();
            }
            inline void write_element(std::ostream& os, const group_elem& g){
                if (this->version >= 2.0){
                    os << g.to_string_2();
                }
                else{
                    os << g.to_string_1();
                }
            }
            inline void write_element(std::ostream& os, const gap_elem& g){
                os << g.to_string();
            }
            inline void write_element(std::ostream& os, const header_elem& h){
                os << h.to_string();
            }

            /** Writes the header to a std::string. */
            inline std::string header_string(){
                std::stringstream st;
                for (auto h : get_header()){
                    write_element(st, h.second);
                }
                return st.str();
            }


            /** Getter methods for elements, to keep users out of our data structures
             *  All of these return a copy of the backing structure in the GFAKluge object. 
             */
            inline std::vector<contained_elem> get_contained(std::string seq_name){
                return seq_to_contained[seq_name];
            }
            inline std::vector<contained_elem> get_contained(sequence_elem seq){
                std::string seq_name = seq.name;
                return seq_to_contained[seq_name];
            }

            inline std::vector<alignment_elem> get_alignments(std::string seq_name){
                return seq_to_alignment[seq_name];
            }
            inline std::vector<alignment_elem> get_alignments(sequence_elem seq){
                return get_alignments(seq.name);
            }

            inline std::map<std::string, header_elem> get_header(){
                return header;
            }
            inline std::map<std::string, sequence_elem, custom_key> get_name_to_seq(){
                return name_to_seq;
            }
            inline std::map<std::string, std::vector<link_elem> > get_seq_to_link(){
                for (auto s : name_to_seq){
                    for (auto e = seq_to_edges[s.first].begin(); e != seq_to_edges[s.first].end(); e++){
                        if (e->determine_type() == 1){
                            link_elem l;
                            l.source_name = e->source_name;
                            l.sink_name = e->sink_name;
                            l.cigar = e->alignment;
                            l.source_orientation_forward = e->source_orientation_forward;
                            l.sink_orientation_forward = e->sink_orientation_forward;
                            l.opt_fields = e->tags;
                            add_link(s.second, l);
                        }

                    }
                }

                return seq_to_link;
            }
            inline std::map<std::string, std::vector<contained_elem> > get_seq_to_contained(){
                for (auto s : name_to_seq){
                    for (auto e = seq_to_edges[s.first].begin(); e != seq_to_edges[s.first].end(); e++){
                        if (e->determine_type() == 2){
                            contained_elem c;
                            c.source_name = e->source_name;
                            c.sink_name = e->sink_name;
                            c.cigar = e->alignment;
                            c.pos = e->source_begin;
                            c.source_orientation_forward = e->source_orientation_forward;
                            c.sink_orientation_forward = e->sink_orientation_forward;
                            c.opt_fields = e->tags;
                            add_contained(s.second, c);
                        }

                    }
                }
                return seq_to_contained;
            }
            inline std::map<std::string, std::vector<alignment_elem> > get_seq_to_alignment(){
                return seq_to_alignment;
            }
            inline std::map<std::string, path_elem> get_name_to_path(){
                return name_to_path;
            }

            // GFA2 getters
            inline std::map<std::string, std::vector<edge_elem>> get_seq_to_edges(){
                return seq_to_edges;
            }
            inline std::map<std::string, std::vector<fragment_elem>> get_seq_to_fragments(){
                return seq_to_fragments;
            }
            inline std::map<std::string, std::vector<gap_elem>> get_seq_to_gaps(){
                return seq_to_gaps;
            }
            inline std::map<std::string, group_elem> get_groups(){
                return groups;
            }

            /** 
             * Convert between GFA1 and GFA2 representations internally.
             * paths/groups are conveted to walks and walks to paths,
             * and edge_elems get their type field set.
             *   **/
            inline void gfa_2_ize(){
                if (!two_compat){
                    // Fix S line length field if needed.
                    for (auto s : name_to_seq){
                        s.second.length = (s.second.sequence != "*" ? (uint64_t) s.second.sequence.length() : s.second.length);
                        // Make an edge for each link
                        for (auto l : seq_to_link[s.first]){
                            edge_elem e;
                            e.type = 1;
                            e.source_name = l.source_name;
                            e.sink_name = l.sink_name;
                            e.source_orientation_forward = l.source_orientation_forward;
                            e.sink_orientation_forward = l.sink_orientation_forward;
                            e.ends.set(0,1);
                            e.ends.set(1,1);
                            e.ends.set(2,0);
                            e.ends.set(3,0);
                            e.alignment = l.cigar;
                            e.tags = l.opt_fields;
                            seq_to_edges[e.source_name].push_back(e);
                        }
                        // Make an edge for each containment
                        for (auto c : seq_to_contained[s.first]){
                            edge_elem e;
                            e.type = 2;
                            e.source_name = c.source_name;
                            e.sink_name = c.sink_name;
                            e.source_orientation_forward = c.source_orientation_forward;
                            e.sink_orientation_forward = c.sink_orientation_forward;
                            e.alignment = c.cigar;

                            std::string overlap = "";
                            int i = 0;
                            while (std::isdigit(e.alignment[i])){
                                overlap += e.alignment[0];
                                ++i;
                            }
                            e.source_end = stoi(overlap) + e.source_begin;
                            e.sink_end = stoi(overlap);

                            if (e.source_end == s.second.length){
                                e.ends.set(0,1);
                                e.ends.set(1,1);
                            }
                            if(e.sink_end == name_to_seq[e.sink_name].length){
                                e.ends.set(2,1);
                                e.ends.set(3,1);
                            }

                            e.tags = c.opt_fields;
                            seq_to_edges[e.source_name].push_back(e);
                        }
                    }
                    // Paths -> ordered groups
                    //walks_as_paths();
                    for (auto p : name_to_path){
                        group_elem g;
                        g.id = p.first;
                        g.ordered = true;
                        g.items = p.second.segment_names;
                        g.orientations = p.second.orientations;

                        g.tags = p.second.opt_fields;
                        add_group(g);
                    }
                    this->two_compat = true;
                }

            }
            inline void gfa_1_ize(){
                if (!one_compat){


                    /**
                     * Set base_edge_id to max_seq_id + 1
                     * Swap edges to links
                     * Ordered groups -> paths,
                     * warn about missing sets
                     * warn about missing fragments (we could output them as seqs)
                     * warn about missing gaps (we could remove them)
                     * 
                     */
                    std::string k = name_to_seq.rbegin()->first;
                    if (string_is_number(k)){
                        base_edge_id = std::stoul(k) + 1;
                    }
                    else{
                        base_edge_id = name_to_seq.size();
                    }
                    groups_as_paths();
                    for (auto s = name_to_seq.begin(); s != name_to_seq.end(); s++){
                        if (s->second.sequence != "*"){
                            s->second.length = s->second.sequence.length();
                        }
                        for (auto e = seq_to_edges[s->first].begin(); e != seq_to_edges[s->first].end(); e++){
                            int t = e->determine_type();
                            if (e->id == "*"){
                                e->id = std::to_string(++base_edge_id);
                            }
                            if (t == 1){
                                e->ends.set(0,1);
                                e->ends.set(1,1);
                                e->ends.set(2,0);
                                e->ends.set(3,0);            
                                e->source_begin = s->second.length;
                                e->source_end = s->second.length;
                                e->sink_begin = 0;
                                e->sink_end = 0;
                                if (e->id == "*"){
                                    e->id = std::to_string(++base_edge_id);
                                }
                            }
                            else if (t == 2){
                                sequence_elem sink = name_to_seq[e->sink_name];
                                e->ends.set(0,(e->source_begin == s->second.length));
                                e->ends.set(1,(e->source_end == s->second.length));
                                e->ends.set(2,(e->sink_begin == sink.length));
                                e->ends.set(3, (e->sink_end == sink.length));
                            }
                            else{
                                std::cerr << "Skipping edge not expressable in GFA2: \"" << e->to_string_2() << "\"" << std::endl;
                            }
                        }
                    }
                    this->one_compat = true;
                }
            }
            // Wraps both gfa_1_ize() and gfa_2_ize()
            inline void compatibilize(){
                gfa_1_ize();
                gfa_2_ize();
            }

            /**
             * Output the entire GFAKluge object as a std::string,
             * outputting the GFA version set on read or by set_version()
             */
            inline std::string to_string(){

                gfa_1_ize();
                gfa_2_ize();
                if (this->version >= 2.0){
                    return to_string_2();
                }

                std::stringstream ret;
                //First print header lines.
                if (header.size() > 0){
                    ret << header_string(header) + "\n";
                }

                if (name_to_path.size() > 0 && this->version >= 1.0){
                    std::map<std::string, path_elem>::iterator pt;
                    for (pt = name_to_path.begin(); pt != name_to_path.end(); ++pt){
                        std::stringstream pat;
                        pat << "P" << "\t";
                        pat << pt->second.name << "\t";
                        std::vector<std::string> ovec;
                        for (size_t seg_ind = 0; seg_ind < pt->second.segment_names.size(); seg_ind++){
                            std::stringstream o_str;
                            o_str << pt->second.segment_names[seg_ind] << (pt->second.orientations[seg_ind] ? "+" : "-");
                            ovec.push_back(o_str.str());
                        }
                        pat << pliib::join(ovec, ",");
                        if (pt->second.overlaps.size() > 0){
                            pat << "\t" << pliib::join(pt->second.overlaps, ",");
                        }
                        pat << "\n";
                        ret << pat.str();
                    }
                }
                for (auto s : name_to_seq){
                    ret << s.second.to_string_1() << std::endl;
                    for (auto e : seq_to_edges[s.first]){
                        ret << e.to_string_1() << std::endl;;
                    }
                    /**
                     *  NB: There are no Fragments in GFA1, so we don't output them.
                     *  We also don't output annotation lines as they're out of spec.
                     *  Also, we check at the function start if we're outputting GFA2,
                     *  so we shouldn't have to do any checks past that point.
                     */

                }


                return ret.str();
            }
            // Force GFA2 std::string output
            inline std::string to_string_2(){
                this->gfa_2_ize();

                std::stringstream ret;
                // Header
                if (header.size() > 0){
                    ret << header_string(header) << std::endl;
                }
                for (auto p : groups){
                    ret << p.second.to_string_2() << std::endl;
                }
                for (auto s : name_to_seq){
                    ret << s.second.to_string_2() << std::endl;
                    for (auto f : seq_to_fragments[s.first]){
                        ret << f.to_string_2() << std::endl;
                    }
                    for (auto e : seq_to_edges[s.first]){
                        ret << e.to_string_2() << std::endl;
                    }
                    for (auto g : seq_to_gaps[s.first]){
                        ret << g.to_string_2() << std::endl;
                    }
                }
                return ret.str();
            }

            /**
             * Output a block_ordered GFA std::string representing the entire
             * GFAKluge object.
             */
            inline std::string block_order_string(){

                this->gfa_1_ize();
                this->gfa_2_ize();

                if (version >= 2.0){
                    return block_order_string_2();
                }
                std::stringstream ret;

                //First print header lines.
                if (header.size() > 0){
                    ret << header_string(header) + "\n";
                }

                std::map<std::string, sequence_elem>::iterator st;

                for (st = name_to_seq.begin(); st != name_to_seq.end(); st++){
                    ret << st->second.to_string_1() << std::endl;
                }



                for (st = name_to_seq.begin(); st != name_to_seq.end(); st++){
                    for (auto e : seq_to_edges[st->first]){
                        if (e.type == 1){
                            ret << e.to_string_1() << std::endl;
                        }

                    }
                }
                for (st = name_to_seq.begin(); st != name_to_seq.end(); st++){
                    for (auto e : seq_to_edges[st->first]){
                        if (e.type == 2){
                            ret << e.to_string_1() << std::endl;
                        }

                    }
                }

                if (name_to_path.size() > 0 && this->version == 1.0){
                    std::map<std::string, path_elem>::iterator pt;
                    for (pt = name_to_path.begin(); pt != name_to_path.end(); ++pt){
                        std::stringstream pat;
                        pat << "P\t" << pt->second.name << "\t";
                        std::vector<std::string> ovec;
                        for (size_t oi = 0; oi < pt->second.segment_names.size(); oi++){
                            std::stringstream o_str;
                            o_str << pt->second.segment_names[oi] << (pt->second.orientations[oi] ? "+" : "-");
                            ovec.push_back(o_str.str());
                        }
                        pat << pliib::join(ovec, ",");
                        if (pt->second.overlaps.size() > 0){
                            pat << "\t" << pliib::join(pt->second.overlaps, ",");
                        }
                        pat << "\n";
                        ret << pat.str();
                    }
                }
                else if (this->version < 1.0){
                    for (auto p : name_to_path){
                        std::stringstream st;
                        p.second.write_as_walks(st);
                    }
                } 
                return ret.str();

            }
            // Force GFA2 std::string output in block order.
            inline std::string block_order_string_2(){
                this->gfa_2_ize();

                std::stringstream ret;
                // Header
                if (header.size() > 0){
                    ret << header_string(header) + "\n";
                }
                // Sequences
                for (auto s : name_to_seq){
                    ret << s.second.to_string_2() << "\n";
                }
                // Fragments
                for (auto s : seq_to_fragments){
                    for (auto f : seq_to_fragments[s.first]){
                        ret << f.to_string_2() << "\n";
                    }

                }
                // Gaps
                for (auto s : name_to_seq){
                    for (auto g : seq_to_gaps[s.first]){
                        ret << g.to_string_2() << "\n";
                    }
                }
                // Edges
                for (auto s : name_to_seq){
                    for (auto e : seq_to_edges[s.first]){
                        ret << e.to_string_2() << "\n";
                    }
                }

                // Paths
                for (auto g : groups){
                    ret << g.second.to_string_2() << "\n";
                }
                return ret.str();
            }

            /**
             * Write the GFAKluge object as GFA0.1/1.0/2.0 to an ostream (e.g. stdout)
             */
            inline void output_to_stream(std::ostream& os, bool output_block_order = false){
                this->gfa_1_ize();
                this->gfa_2_ize();
                if (this->version == 2.0 && output_block_order){
                    // Header
                    if (header.size() > 0){
                        os << header_string(header) + "\n";
                    }
                    // Sequences
                    for (auto s : name_to_seq){
                        os << s.second.to_string_2() << "\n";
                    }
                    // Fragments
                    for (auto s : seq_to_fragments){
                        for (auto f : seq_to_fragments[s.first]){
                            os << f.to_string_2() << "\n";
                        }

                    }
                    // Gaps
                    for (auto s : name_to_seq){
                        for (auto g : seq_to_gaps[s.first]){
                            os << g.to_string_2() << "\n";
                        }
                    }
                    // Edges
                    for (auto s : name_to_seq){
                        for (auto e : seq_to_edges[s.first]){
                            os << e.to_string_2() << "\n";
                        }
                    }

                    // Paths
                    for (auto g : groups){
                        os << g.second.to_string_2() << "\n";
                    }

                }
                else if (this->version == 2.0){
                    if (header.size() > 0){
                        os << header_string(header) << std::endl;
                    }
                    for (auto p : groups){
                        os << p.second.to_string_2() << std::endl;
                    }
                    for (auto s : name_to_seq){
                        os << s.second.to_string_2() << std::endl;
                        for (auto f : seq_to_fragments[s.first]){
                            os << f.to_string_2() << std::endl;
                        }
                        for (auto e : seq_to_edges[s.first]){
                            os << e.to_string_2() << std::endl;
                        }
                        for (auto g : seq_to_gaps[s.first]){
                            os << g.to_string_2() << std::endl;
                        }
                    }
                }
                else if (this->version < 2.0 && output_block_order){
                    //First print header lines.
                    if (header.size() > 0){
                        os << header_string(header) + "\n";
                    }
                    std::map<std::string, sequence_elem>::iterator st;
                    for (st = name_to_seq.begin(); st != name_to_seq.end(); st++){
                        os << st->second.to_string_1() << std::endl;
                    }

                    for (st = name_to_seq.begin(); st != name_to_seq.end(); st++){
                        for (auto e : seq_to_edges[st->first]){
                            if (e.type == 1){
                                os << e.to_string_1() << std::endl;
                            }

                        }
                    }
                    for (st = name_to_seq.begin(); st != name_to_seq.end(); st++){
                        for (auto e : seq_to_edges[st->first]){
                            if (e.type == 2){
                                os << e.to_string_1() << std::endl;
                            }

                        }
                    }

                    if (name_to_path.size() > 0 && this->version >= 1.0){
                        std::map<std::string, path_elem>::iterator pt;
                        for (pt = name_to_path.begin(); pt != name_to_path.end(); ++pt){
                            std::stringstream pat;
                            pat << "P\t" << pt->second.name << "\t";
                            std::vector<std::string> ovec;
                            for (size_t oi = 0; oi < pt->second.segment_names.size(); oi++){
                                std::stringstream o_str;
                                o_str << pt->second.segment_names[oi] << (pt->second.orientations[oi] ? "+" : "-");
                                ovec.push_back(o_str.str());
                            }
                            pat << pliib::join(ovec, ",");
                            if (pt->second.overlaps.size() > 0){
                                pat << "\t" << pliib::join(pt->second.overlaps, ",");
                            }
                            pat << "\n";
                            os << pat.str();
                        }
                    }
                    else if (this->version < 1.0 && name_to_path.size() > 0){
                        std::map<std::string, path_elem>::iterator pt;
                        for (pt = name_to_path.begin(); pt != name_to_path.end(); pt++){
                            pt->second.write_as_walks(os);
                        }
                    }


                }
                else if (this->version < 2.0){
                    //First print header lines.
                    if (header.size() > 0){
                        os << header_string(header) + "\n";
                    }

                    if (name_to_path.size() > 0 && this->version >= 1.0){
                        std::map<std::string, path_elem>::iterator pt;
                        for (pt = name_to_path.begin(); pt != name_to_path.end(); ++pt){
                            std::stringstream pat;
                            pat << "P" << "\t";
                            pat << pt->second.name << "\t";
                            std::vector<std::string> ovec;
                            for (size_t seg_ind = 0; seg_ind < pt->second.segment_names.size(); seg_ind++){
                                std::stringstream o_str;
                                o_str << pt->second.segment_names[seg_ind] << (pt->second.orientations[seg_ind] ? "+" : "-");
                                ovec.push_back(o_str.str());
                            }
                            pat << pliib::join(ovec, ",");
                            if (pt->second.overlaps.size() > 0){
                                pat << "\t" << pliib::join(pt->second.overlaps, ",");
                            }
                            pat << "\n";
                            os << pat.str();
                        }
                    }
                    else if (this->version < 1.0 && name_to_path.size() > 0){
                        std::map<std::string, path_elem>::iterator pt;
                        for (pt = name_to_path.begin(); pt != name_to_path.end(); pt++){
                            pt->second.write_as_walks(os);
                        }
                    }

                    for (auto s : name_to_seq){
                        os << s.second.to_string_1() << std::endl;
                        for (auto e : seq_to_edges[s.first]){
                            os << e.to_string_1() << std::endl;;
                        }
                        /**
                         *  NB: There are no Fragments in GFA1, so we don't output them.
                         *  We also don't output annotation lines as they're out of spec.
                         *  Also, we check at the function start if we're outputting GFA2,
                         *  so we shouldn't have to do any checks past that point.
                         */

                    }


                }

                else{
                    std::cerr << "Invalid version " << this->version << std::endl;
                    exit(9);

                }
            }


            // ID manipulators
            /** Return the highest IDs present in this GFAKluge object
             *  for sequence_elems, edge_elems,
             *  fragment_elems, gap_elems, and group_elems.
             */
            inline std::tuple<uint64_t, uint64_t, uint64_t, uint64_t, uint64_t> max_ids(){
                return std::make_tuple(this->base_seq_id, this->base_edge_id,
                        this->base_frag_id, this->base_gap_id, this->base_group_id);
            }
            /** Return the highest IDs (as above) but as a colon-delimited std::string rather
             *  than a tuple.
             */
            inline std::string max_ids_string(){
                std::tuple<std::uint64_t, std::uint64_t, std::uint64_t, std::uint64_t, std::uint64_t> x = max_ids();
                std::vector<std::string> max_str(5);
                max_str[0] = std::to_string(std::get<0>(x));
                max_str[1] = std::to_string(std::get<1>(x));
                max_str[2] = std::to_string(std::get<2>(x));
                max_str[3] = std::to_string(std::get<3>(x));
                max_str[4] = std::to_string(std::get<4>(x));
                return pliib::join(max_str, ":");
            }
            /** Bump the IDs of sequence-, edge-, fragment-, gap-, and group_elems to 
             *  be greater than new_mx. Useful for concatenating graphs.
             */
            inline void re_id(std::tuple<uint64_t, uint64_t, uint64_t, uint64_t, uint64_t>& new_mx){

                base_seq_id = std::get<0>(new_mx);
                base_edge_id = std::get<1>(new_mx);
                base_frag_id = std::get<2>(new_mx);
                base_gap_id = std::get<3>(new_mx);
                base_group_id = std::get<4>(new_mx);
                std::uint64_t seg_diff = base_seq_id;
                // Segments
                //int num_segs = name_to_seq.size();

                // locally cache name_to_seq,
                // seq_to_edges, seq_to_fragments, groups,
                // and seq_to_gaps before clearing them.
                std::map<std::string, sequence_elem, custom_key> n_s;
                std::map<std::string, std::vector<edge_elem>> s_e;
                std::map<std::string, std::vector<fragment_elem>> s_f;
                std::map<std::string, std::vector<gap_elem>> s_g;
                std::map<std::string, group_elem> g_g;

                for (auto ns : name_to_seq){
                    std::string old_name = ns.second.name;
                    ns.second.id = ++base_seq_id;
                    ns.second.name = std::to_string(ns.second.id);    


                    // Paths



                    //uint64_t edge_count = 0;
                    // Edges
                    for (auto e : seq_to_edges[old_name]){
                        e.source_name = ns.second.name;
                        e.sink_name = std::to_string(seg_diff + stoul(e.sink_name));
                        e.id = std::to_string(++base_edge_id);
                        //add_edge(e.source_name, e);
                        s_e[e.source_name].push_back(e);
                    }
                    // Fragments
                    for (auto f : seq_to_fragments[old_name]){
                        f.id = ns.second.name;
                        //add_fragment(f.id, f);
                        s_f[f.id].push_back(f);
                    }
                    // Gaps
                    for (auto g : seq_to_gaps[old_name]){
                        g.id = std::to_string(++base_gap_id);
                        g.source_name = ns.second.name;
                        g.sink_name = std::to_string(seg_diff + stoul(g.sink_name));
                        s_g[g.source_name].push_back(g);
                        //add_gap(g);
                    }

                    // Groups
                    for (auto g : groups){
                        //We may not want to increment group IDs, as they might represent unique paths!
                        //g.second.id = to_string(++base_group_id);
                        for (size_t i = 0; i < g.second.items.size(); ++i){
                            g.second.items[i] = std::to_string(seg_diff + stoul(g.second.items[i]));
                        }
                        g_g[g.first] = g.second;
                        //add_group(g.second);
                    }

                    //add_sequence(ns.second);

                    n_s[ns.second.name] = ns.second;
                }

                // Clear our internal std::maps
                name_to_seq = n_s;
                seq_to_edges = s_e;
                seq_to_fragments = s_f;
                seq_to_gaps = s_g;
                groups = g_g;

            }
            /** Bump the IDs of each type of GFA2 element so that the lowest ID
             *  for each type is defined by new_mx, which is a colon-delimited std::string.
             */
            inline void re_id(std::string new_mx_str){
                std::vector<uint64_t> starts(5);
                std::vector<std::string> starts_strs = pliib::split(new_mx_str, ':');
                for (size_t i = 0; i < starts_strs.size(); ++i){
                    starts[i] = stoul(starts_strs[i]);
                }
                std::tuple<std::uint64_t, std::uint64_t, std::uint64_t, std::uint64_t, std::uint64_t> n_ids = std::make_tuple(starts[0], starts[1],
                        starts[2], starts[3], starts[4]);
                re_id(n_ids);

            }

            /** Merge two GFA graphs **/
            // TODO check incompatible version numbers
            // TODO Check colliding groups, headers
            inline void merge(GFAKluge& gg){
                std::unordered_set<std::string> seg_ids;
                // Merge headers
                for (auto h : gg.get_header()){
                    header[h.first] = h.second;
                }
                for (auto s : this->get_name_to_seq()){
                    seg_ids.insert(s.first);
                }

                std::map<std::string, sequence_elem, custom_key> ss = gg.get_name_to_seq();
                std::map<std::string, std::vector<fragment_elem>> sf = gg.get_seq_to_fragments();
                std::map<std::string, std::vector<gap_elem>> sg = gg.get_seq_to_gaps();
                std::map<std::string, std::vector<edge_elem>> se = gg.get_seq_to_edges();

                if (this->two_compat){
                    for (auto s : ss){
                        if (!seg_ids.count(s.second.name)){
                            this->add_sequence(s.second);
                            for (auto e : se){
                                seq_to_edges.insert(e);
                                //this->add_edge(e.source_name, e);
                            }
                            for (auto g : sg[s.first]){
                                this->add_gap(g);
                            }
                            for (auto f : sf[s.first]){
                                this->add_fragment(s.first, f);
                            }
                        }
                        else{
                            std::cerr << "WARNING: DUPLICATE IDS " << s.second.name << std::endl <<
                                " will be lost." << std::endl;
                        }
                        for (auto g : gg.get_groups()){
                            this->add_group(g.second);
                        }
                    }
                }
                else if (one_compat){

                }

            }

            /** Assembly stats **/
            // we calculate the N50 based on the 'S' lines,
            // Though in theory an O line might also be a contig
            inline double get_N50(){
                std::vector<double> s_lens;
                uint64_t total_len = 0;
                double n = 0.0;
                for (auto s = name_to_seq.begin(); s != name_to_seq.end(); s++){
                    s_lens.push_back(s->second.length);
                    total_len += s->second.length;
                }
                double avg = total_len * 0.50 ;
                std::sort(s_lens.begin(), s_lens.end());
                //int middle = floor((double) s_lens.size() / 2.0);
                double cumul_size = 0.0;
                for (size_t i = 0; i < s_lens.size(); i++){
                    cumul_size += s_lens[i];
                    if (cumul_size >= avg){
                        n = s_lens[i];
                        break;
                    }
                }
                return n;
            }
            inline double get_N90(){
                std::vector<double> s_lens;
                uint64_t total_len = 0;
                double n = 0.0;
                for (auto s = name_to_seq.begin(); s != name_to_seq.end(); s++){
                    s_lens.push_back(s->second.length);
                    total_len += s->second.length;
                }
                double avg = total_len * 0.90;
                std::sort(s_lens.rbegin(), s_lens.rend());
                //int middle = floor((double) s_lens.size() / 2.0);
                double cumul_size = 0.0;
                for (size_t i = 0; i < s_lens.size(); i++){
                    cumul_size += s_lens[i];
                    if (cumul_size >= avg){
                        n = s_lens[i];
                        break;
                    }
                }
                return n;
            }

            inline int get_L50(){
                std::vector<double> s_lens;
                uint64_t total_len = 0;
                for (auto s = name_to_seq.begin(); s != name_to_seq.end(); s++){
                    s_lens.push_back(s->second.length);
                    total_len += s->second.length;
                }
                double avg = total_len * 0.50;
                std::sort(s_lens.rbegin(), s_lens.rend());
                //int middle = floor((double) s_lens.size() / 2.0);
                double cumul_size = 0.0;
                for (size_t i = 0; i < s_lens.size(); i++){
                    cumul_size += s_lens[i];
                    if (cumul_size >= avg){
                        return i + 1;
                    }
                }
                return -1;
            }

            inline int get_L90(){
                std::vector<double> s_lens;
                uint64_t total_len = 0;
                for (auto s = name_to_seq.begin(); s != name_to_seq.end(); s++){
                    s_lens.push_back(s->second.length);
                    total_len += s->second.length;
                }
                double avg = total_len * 0.90;
                std::sort(s_lens.rbegin(), s_lens.rend());
                //int middle = floor((double) s_lens.size() / 2.0);
                double cumul_size = 0.0;
                for (size_t i = 0; i < s_lens.size(); i++){
                    cumul_size += s_lens[i];
                    if (cumul_size >= avg){
                        return i + 1;
                    }
                }
                return -1;
            }
            // uint64_t num_contigs();
            // double simple_connectivity() // reports avg edges / sequence
            // double weighted_connectivity() // weight areas of high connectivity higher
            //      behaves kind of like [(N_edges)^2 / (N_seqs)^2 if N_edges > 2]

            /** Given the name of a FASTA file,
             *  fill in the sequence field of each sequence_elem with an entry from
             *  that file with a correspondinmg name. If no entry is present, maintain
             *  the "*" placeholder that should be present in that element's sequence field.
             */
            inline void fill_sequences(const char* fasta_file){

                TFA::tiny_faidx_t tf;
                if (TFA::checkFAIndexFileExists(fasta_file)){
                    TFA::parseFAIndex(fasta_file, tf);
                }
                else{
                    std::cerr << "Creating index for " << fasta_file << "." << std::endl;
                    TFA::createFAIndex(fasta_file, tf);
                    std::cerr << "Created index." << std::endl;
                    TFA::writeFAIndex(fasta_file, tf);
                    std::cerr << "Wrote index to file." << std::endl;
                }

                for (std::map<std::string, sequence_elem, custom_key>::iterator it = name_to_seq.begin(); it != name_to_seq.end(); it++){
                    if (tf.hasSeqID(it->second.name.c_str())){
                        char* s;
                        TFA::getSequence(tf, it->second.name.c_str(), s);
                        it->second.sequence.assign(s);
                        // The length field should already be filled, but it
                        // might be good to check.
                        delete [] s;
                    }
                }
            }

            /**
             * Remove any sequence_elems (and any edges connected to them) that have
             * a sequence shorter than a certain length.
             * Returns true if the graph is modified.
             */
            inline bool trim_seqs(const int& minlen = 0, const bool& no_ambiguous = false){
                bool graph_modified = false;

                std::unordered_set<std::string> dropped_seqs;
                std::map<std::string, sequence_elem, custom_key>::iterator n_to_s;
                for (n_to_s = name_to_seq.begin(); n_to_s != name_to_seq.end(); n_to_s++){
                    auto& s = n_to_s->second;
                    if (s.length == UINT64_MAX){
                        std::cerr << "Length unset for sequence " << s.name << "; removing from graph." << std::endl;
                        dropped_seqs.insert(s.name);
                        name_to_seq.erase(n_to_s);
                        graph_modified = true;
                    }
                    else if (minlen > 0 && s.length < (size_t)minlen){
                        dropped_seqs.insert(s.name);
                        name_to_seq.erase(n_to_s);
                        graph_modified = true;
                    }

                    if (no_ambiguous && !pliib::canonical(s.sequence)){
                        dropped_seqs.insert(s.name);
                        name_to_seq.erase(n_to_s);
                        graph_modified = true;
                    }
                }

                std::map<std::string, std::vector<edge_elem>>::iterator s_to_e;
                for (s_to_e = seq_to_edges.begin(); s_to_e != seq_to_edges.end(); s_to_e++){
                    // Remove all edges that have dropped sequences as their source.
                    if (dropped_seqs.find(s_to_e->first) != dropped_seqs.end()){
                        seq_to_edges.erase(s_to_e);
                        graph_modified = true;
                        continue;
                    }
                    // Remove any individual edges that have a dropped sequence as their sink.
                    std::vector<edge_elem> pass_elems;
                    for (auto ee : s_to_e->second){
                        if (dropped_seqs.find(ee.source_name) == dropped_seqs.end()){
                            graph_modified = true;
                            pass_elems.push_back(ee);
                        }
                    }
                    seq_to_edges[s_to_e->first] = pass_elems;
                }


                return graph_modified;
            }





            /** Parse a GFA file to a GFAKluge object. */
            inline bool parse_gfa_file(const std::string &filename) {
                std::ifstream gfi;
                gfi.open(filename.c_str(), std::ifstream::in);
                if (!gfi.good()){
                    std::cerr << "Couldn't open GFA file " << filename << "." << std::endl;
                    exit(1);
                }

                bool ret = parse_gfa_file(gfi);
                gfi.close();
                return ret;

            }
            inline bool parse_gfa_file(std::istream& instream){
                std::string line;
                std::vector<std::string> line_tokens;
                while (getline(instream, line)){
                    std::vector<std::string> tokens = pliib::split(line, '\t');
                    if (tokens[0] == "H"){
                        header_elem h;
                        line_tokens = pliib::split(tokens[1], ':');
                        //TODO this is not well implemented
                        // GFA places no guarantees on header format
                        h.key = line_tokens[0];
                        h.type = line_tokens[1];
                        h.val = line_tokens[2];
                        if (h.key.compare("VN") == 0){
                            set_version(stod(h.val));
                        }

                        header[h.key] = h;
                    }
                    else if (tokens[0] ==  "S"){
                        //TODO: we've got some tokens at the end of the line
                        //that have not been handled yet.
                        int tag_index = 3;
                        sequence_elem s;
                        s.name = tokens[1];

                        if (this->version >= 2.0 || string_is_number(tokens[2])){
                            s.length = stoi(tokens[2]); 
                            s.sequence = tokens[3];
                            tag_index = 4;
                        }
                        else{

                            s.sequence = tokens[2];
                            if (tokens[2] == "*"){
                                s.length = UINT64_MAX;
                            }
                            else{
                                s.length = s.sequence.length();
                            }
                            tag_index = 3;
                        }

                        size_t i;
                        if (tokens.size() > 3){
                            for (i = tag_index; i < tokens.size(); i++){
                                //opt fields are in key:type:val format
                                std::vector<std::string> opt_field = pliib::split(tokens[i], ':');
                                opt_elem o;
                                o.key = opt_field[0];
                                o.type = opt_field[1];
                                o.val = pliib::join(std::vector<std::string> (opt_field.begin() + 2, opt_field.end()), ":");
                                s.opt_fields.push_back(o);
                                if (o.key == "LN" && s.length == UINT64_MAX){
                                    s.length = stoul(o.val);
                                }
                            }
                        }
                        name_to_seq[s.name] = s;
                    }
                    else if (tokens[0] == "E"){
                        edge_elem e;

                        e.id = tokens[1];

                        std::string x = tokens[2];
                        e.source_name = x.substr(0, x.length() - 1);
                        e.source_orientation_forward = (x.back() == '+');

                        x = tokens[3];
                        e.sink_name = x.substr(0, x.length() - 1);
                        e.sink_orientation_forward = (x.back() == '+');

                        x = tokens[4];
                        e.ends.set(0, (x.back() == '$' ? 1 : 0));
                        e.source_begin = (e.ends.test(0) ? stoul(x.substr(0, x.length() - 1)) : stoul(x));

                        x = tokens[5];
                        e.ends.set(1, (x.back() == '$' ? 1 : 0));
                        e.source_end = (e.ends.test(1) ? stoul(x.substr(0, x.length() - 1)) : stoul(x));

                        x = tokens[6];
                        e.ends.set(2, (x.back() == '$' ? 1 : 0));
                        e.sink_begin = (e.ends.test(2) ? stoul(x.substr(0, x.length() - 1)) : stoul(x));


                        x = tokens[7];
                        e.ends.set(3, (x.back() == '$' ? 1 : 0));
                        e.sink_end = (e.ends.test(3) ? stoul(x.substr(0, x.length() - 1)) : stoul(x));

                        e.alignment = tokens[8];

                        if (tokens.size() > 9){
                            for (size_t i = 9; i < tokens.size(); i++){
                                //opt fields are in key:type:val format
                                std::vector<std::string> opt_field = pliib::split(tokens[i], ':');
                                opt_elem o;
                                o.key = opt_field[0];
                                o.type = opt_field[1];
                                o.val = pliib::join(std::vector<std::string> (opt_field.begin() + 2, opt_field.end()), ":");
                                e.tags[o.key] = o;

                            }
                        }
                        add_edge(e.source_name, e);

                    }
                    else if (tokens[0] == "G"){
                        // <- G <gid:opt_id> <sid1:ref> <sid2:ref> <dist:int> (* | <var:int>) <tag>*
                        gap_elem g;
                        g.id = tokens[1];
                        g.source_name = tokens[2];
                        g.sink_name = tokens[3];
                        g.distance = stoul(tokens[4]);
                        add_gap(g);

                    }
                    else if (tokens[0] == "F"){
                        fragment_elem f;
                        f.id = tokens[1];
                        f.ref = tokens[2].substr(0, tokens[2].length() - 1);
                        f.ref_orientation = (tokens[2].back() == '+');
                        f.seg_begin = stoul(tokens[3]);
                        f.seg_end = stoul(tokens[4]);
                        f.frag_begin = stoul(tokens[5]);
                        f.frag_end = stoul(tokens[6]);
                        f.ends.set(0, (tokens[3].back() == '$' ? 1 : 0));
                        f.ends.set(1, (tokens[4].back() == '$' ? 1 : 0));
                        f.ends.set(2, (tokens[5].back() == '$' ? 1 : 0));
                        f.ends.set(3, (tokens[6].back() == '$' ? 1 : 0));
                        f.alignment = tokens[7];
                        if (tokens.size() > 8){
                            for (size_t i = 9; i < tokens.size(); i++){
                                //opt fields are in key:type:val format
                                std::vector<std::string> opt_field = pliib::split(tokens[i], ':');
                                opt_elem o;
                                o.key = opt_field[0];
                                o.type = opt_field[1];
                                o.val = pliib::join(std::vector<std::string> (opt_field.begin() + 2, opt_field.end()), ":");
                                f.tags[o.key] = o;

                            }
                        }
                        add_fragment(f.id, f);
                    }
                    else if (tokens[0] == "O"){
                        group_elem g;
                        g.ordered = true;
                        g.id = tokens[1];
                        if (g.id == "*"){
                            g.id = std::to_string(++base_group_id);
                        }
                        std::vector<std::string> g_ids = pliib::split(tokens[2], ' ');
                        for (size_t i = 0 ; i < g_ids.size(); i++){
                            g.items.push_back(g_ids[i].substr(0, g_ids[i].length() - 1));
                            g.orientations.push_back(g_ids[i].back() == '+');
                        }
                        if (tokens.size() > 8){
                            for (std::size_t i = 9; i < tokens.size(); i++){
                                //opt fields are in key:type:val format
                                std::vector<std::string> opt_field = pliib::split(tokens[i], ':');
                                opt_elem o;
                                o.key = opt_field[0];
                                o.type = opt_field[1];
                                o.val = pliib::join(std::vector<std::string> (opt_field.begin() + 2, opt_field.end()), ":");
                                g.tags[o.key] = o;

                            }
                        }
                        add_group(g);
                    }
                    else if (tokens[0] == "U"){
                        group_elem g;
                        g.ordered = false;
                        g.id = tokens[0];
                        if (g.id == "*"){
                            g.id = std::to_string(++base_group_id);
                        }
                        g.items = pliib::split(tokens[2], ' ');
                        if (tokens.size() > 8){
                            for (size_t i = 9; i < tokens.size(); i++){
                                //opt fields are in key:type:val format
                                std::vector<std::string> opt_field = pliib::split(tokens[i], ':');
                                opt_elem o;
                                o.key = opt_field[0];
                                o.type = opt_field[1];
                                o.val = pliib::join(std::vector<std::string> (opt_field.begin() + 2, opt_field.end()), ":");
                                g.tags[o.key] = o;

                            }
                        }
                        add_group(g);
                    }
                    else if (tokens[0] ==  "L"){
                        // TODO: we need to deal with  where the link is given before
                        // its corresponding sequence in the file. TODO this is probably
                        // now fixed by using the std::string: sequence std::map.
                        edge_elem e;
                        e.type = 1;
                        e.source_name = tokens[1];
                        e.sink_name = tokens[3];
                        //TODO: search the input std::strings for "-" and "+" and set using ternary operator
                        e.source_orientation_forward = tokens[2] == "+" ? true : false;
                        e.ends.set(0, 1);
                        e.ends.set(1,1);
                        e.ends.set(2,0);
                        e.ends.set(3, 0);
                        e.sink_orientation_forward = tokens[4] == "+" ? true : false;
                        if (tokens.size() >= 6){
                            e.alignment = tokens[5];
                        }
                        else{
                            e.alignment = "*";
                        }

                        if (tokens.size() >= 7){
                            for (size_t i = 6; i < tokens.size(); i++){
                                //opt fields are in key:type:val format
                                std::vector<std::string> opt_field = pliib::split(tokens[i], ':');
                                opt_elem o;
                                o.key = opt_field[0];
                                o.type = opt_field[1];
                                o.val = pliib::join(std::vector<std::string> (opt_field.begin() + 2, opt_field.end()), ":");
                                e.tags[o.key] = o;

                            }
                        }

                        add_edge(e.source_name, e);
                    }
                    else if (tokens[0] == "C"){
                        //contained_elem c;
                        edge_elem e;
                        e.type = 2;
                        e.source_name = tokens[1];
                        e.sink_name = tokens[3];
                        e.source_orientation_forward = tokens[2] == "+" ? true : false;
                        e.sink_orientation_forward = tokens[4] == "+" ? true : false;
                        e.sink_begin = 0;
                        e.source_begin = stoul(tokens[5]);
                        e.ends.set(3, 1);
                        if (tokens.size() > 6){
                            e.alignment = tokens[6];
                            std::string overlap = "";
                            int i = 0;
                            while (std::isdigit(e.alignment[i])){
                                overlap += e.alignment[0];
                                ++i;
                            }
                            e.source_end = stoi(overlap) + e.source_begin;
                            e.sink_end = stoi(overlap);
                        }
                        else{
                            e.alignment = "*";
                        }

                        if (tokens.size() >= 8){
                            for (size_t i = 8; i < tokens.size(); i++){
                                //opt fields are in key:type:val format
                                std::vector<std::string> opt_field = pliib::split(tokens[i], ':');
                                opt_elem o;
                                o.key = opt_field[0];
                                o.type = opt_field[1];
                                o.val = pliib::join(std::vector<std::string> (opt_field.begin() + 2, opt_field.end()), ":");
                                e.tags[o.key] = o;

                            }
                        }


                        add_edge(e.source_name, e);
                    }
                    else if (tokens[0] == "W"){
                        std::string pname(tokens[2]);
                        std::string wname(tokens[1]);
                        int rank;
                        bool ori;
                        std::string overlap;
                        std::vector<opt_elem> opts;

                        if (tokens[3].compare("+") == 0 || tokens[3].compare("-") == 0){
                            rank = 0;
                            ori = (tokens[3] == "-" ? false : true);
                            overlap = tokens[4];
                        }
                        else{
                            rank = atoi(tokens[3].c_str());
                            ori = (tokens[4] == "-" ? false : true);
                            overlap = tokens[5];
                        }

                        add_walk(pname, rank, wname, ori, overlap, opts);
                    }
                    else if (tokens[0] == "P"){
                        if (this->version >= 1.0 && this->version < 2.0){
                            // Parse a GFA 1.0 path element
                            path_elem p;
                            p.name = tokens[1];
                            std::vector<std::string> ids_and_orientations;
                            pliib::split(tokens[2], ',', ids_and_orientations);
                            p.segment_names.resize(ids_and_orientations.size());
                            p.orientations.resize(ids_and_orientations.size());
                            for (size_t t = 0; t < ids_and_orientations.size(); ++t){
                                std::string x = ids_and_orientations[t];
                                bool orientation = ((x.back()) == '+' || x.front() == '+');
                                std::string id = x.substr(0, x.length() - 1);
                                p.segment_names[t] = id;
                                p.orientations[t] = orientation;
                            }


                            if (tokens.size() > 3){
                                std::vector<std::string> spltz = pliib::split(tokens[3], ',');
                                for (auto z : spltz){
                                    p.overlaps.push_back(z);
                                }
                                //p.overlaps.assign( spltz.begin(), spltz.end());
                            }
                            else{
                                //p.overlaps.assign(p.segment_names.size(), "*");
                                for (auto z : p.segment_names){
                                    p.overlaps.push_back("*");
                                }
                            }

                            name_to_path[p.name] = p;
                        }
                        else if (this->version < 1.0){
                            std::string pname(tokens[2]);
                            std::string wname(tokens[1]);
                            int rank;
                            bool ori;
                            std::string overlap;
                            std::vector<opt_elem> opts;

                            if (tokens[3].compare("+") == 0 || tokens[3].compare("-") == 0){
                                // TODO: this is a hack
                                // if no rank is present, insert at the end
                                rank = 0;
                                ori = (tokens[3].compare("+") == 0);
                                overlap = tokens[4];
                            }
                            else{
                                rank = atoi(tokens[3].c_str());
                                ori = (tokens[4] == "+");
                                overlap = tokens[5];
                            }
                            add_walk(pname, rank, wname, ori, overlap, opts);
                        }
                        else{
                            std::cerr << "Cannot parse; version of GFA is too new. Version: " << this->version << std::endl;
                            exit(-1);
                        }

                    }
                    else if (tokens[0] == "x"){
                        annotation_elem x;
                        x.key = tokens[1];
                        x.info = tokens[2];
                    }
                    else if (tokens[0] == "a"){
                        alignment_elem a;
                        a.source_name = tokens[1];
                        a.position = atoi(tokens[2].c_str());
                        a.ref = tokens[3];
                        a.source_orientation_forward = tokens[4] == "+" ? true : false;
                        a.length = atoi(tokens[5].c_str());
                        add_alignment(a.source_name, a);
                    }
                    else if (tokens[0] == "#"){
                        continue;
                    }
                    else{
                        std::cerr << "Unknown line identifier  encountered: " << tokens[0] <<  " . Exiting." << std::endl;
                        exit(1);
                    }

                }
                gfa_1_ize();
                gfa_2_ize();


                return true;

            }
        };

    }
#endif
