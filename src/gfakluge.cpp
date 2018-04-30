#include "gfakluge.hpp"

#include <unordered_set>
#include <fstream>

using namespace std;
namespace gfak{

    int determine_line_type(char* line){
       if (line[0] == 'H'){
            return HEADER_LINE;
       }
       else if (line[0] == 'S'){
            return SEGMENT_LINE;
       }
       else if (line[0] == 'E' || line[0] == 'L' || line[0] == 'C'){
            return EDGE_LINE;
       }
       else if (line[0] == 'U' || line[0] == 'O' || line[0] == 'P' || line[0] == 'W'){
            return GROUP_LINE;
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

    map<string, vector<edge_elem>> GFAKluge::get_seq_to_edges(){
        return seq_to_edges;
    }
    map<string, vector<fragment_elem>> GFAKluge::get_seq_to_fragments(){
        return seq_to_fragments;
    }
    map<string, vector<gap_elem>> GFAKluge::get_seq_to_gaps(){
        return seq_to_gaps;
    }
    map<string, group_elem> GFAKluge::get_groups(){
        return groups;
    }

    void GFAKluge::add_edge(const string& seqname, const edge_elem& e){
        seq_to_edges[seqname].push_back(e);
    }
    void GFAKluge::add_edge(const sequence_elem& s, const edge_elem& e){
        add_edge(s.name, e);
    }

    void GFAKluge::add_fragment(const string& seqname, const fragment_elem& f){
        seq_to_fragments[seqname].push_back(f);
    }

    void GFAKluge::add_fragment(const sequence_elem& s, const fragment_elem& f){
        add_fragment(s.name, f);
    }

    void GFAKluge::add_gap(const gap_elem& g){
        seq_to_gaps[g.source_name].push_back(g);
    }

    void GFAKluge::add_group(const group_elem& g){
        this->groups[g.id] = g;
    }

    void GFAKluge::groups_as_paths(){
        for (auto g : groups){
            if (g.second.ordered){
                path_elem p;
                p.name = g.first;
                p.segment_names = g.second.items;
                p.orientations = g.second.orientations;
                p.opt_fields = g.second.tags;
                add_path(p.name, p);
            }
            else{
                cerr << "Group " << g.first << " is unordered; skipping adding it to the paths." << endl;
            }
        }
    }

    void GFAKluge::gfa_1_ize(){
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
            string k = name_to_seq.rbegin()->first;
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
                        cerr << "Skipping edge not expressable in GFA2: \"" << e->to_string_2() << "\"" << endl;
                    }
                }
            }
            this->one_compat = true;
        }
    }

    void GFAKluge::gfa_2_ize(){
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

                    string overlap = "";
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
            walks_as_paths();
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



    GFAKluge::GFAKluge(){
        map<std::string, header_elem> header;
        map<std::string, vector<contained_elem>, custom_key > seq_to_contained;
        map<std::string, vector<link_elem>, custom_key > seq_to_link;
        map<std::string, vector<alignment_elem> , custom_key> seq_to_alignment;
        //Since we can't compare sequence elements for hashing,
        // we cheat and use their names (which are sort of guaranteed to be
        // unique.
        map<string, sequence_elem, custom_key> name_to_seq;
        map<string, vector<walk_elem>, custom_key > seq_to_walks;
        map<string, path_elem> name_to_path;
        //cout << fixed << setprecision(2);


        /** GFA 2.0 containers **/
        map<std::string, vector<fragment_elem> > seq_to_fragments;
        map<std::string, vector<edge_elem> > seq_to_edges;
        map<std::string, vector<gap_elem> > seq_to_gaps;
        map<string, group_elem> groups;
    }

    GFAKluge::~GFAKluge(){

    }

    double GFAKluge::get_version(){
        return this->version;
    }

    void GFAKluge::set_version(double v){
        header_elem verz;
        verz.key = "VN";
        verz.type="Z";
        this->version = v;
        verz.val = std::to_string((double) this->version).substr(0,3);
        this->header[verz.key] = verz;
        //gfa_1_ize();
        //gfa_2_ize();
    }
    void GFAKluge::set_version(){
        header_elem verz;
        verz.key = "VN";
        verz.type="Z";
        verz.val = std::to_string((double) this->version).substr(0,3);

        this->header[verz.key] = verz;
    }

    bool GFAKluge::string_is_number(string s){
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

    bool GFAKluge::parse_gfa_file(const std::string &filename) {
        ifstream gfi;
        gfi.open(filename.c_str(), std::ifstream::in);
        if (!gfi.good()){
            cerr << "Invalid input stream. Exiting." << endl;
            return false;
        }

        return parse_gfa_file(gfi);

    }

    bool GFAKluge::parse_gfa_file(istream& instream){
        string line;
        vector<string> line_tokens;
        while (getline(instream, line)){
            vector<string> tokens = pliib::split(line, '\t');
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
                //s.id = atol(s.name.c_str());
                int i;
                if (tokens.size() > 3){
                    for (i = tag_index; i < tokens.size(); i++){
                        //opt fields are in key:type:val format
                        vector<string> opt_field = pliib::split(tokens[i], ':');
                        opt_elem o;
                        o.key = opt_field[0];
                        o.type = opt_field[1];
                        o.val = join(vector<string> (opt_field.begin() + 2, opt_field.end()), ":");
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
                
                string x = tokens[2];
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
                    for (int i = 9; i < tokens.size(); i++){
                         //opt fields are in key:type:val format
                        vector<string> opt_field = pliib::split(tokens[i], ':');
                        opt_elem o;
                        o.key = opt_field[0];
                        o.type = opt_field[1];
                        o.val = join(vector<string> (opt_field.begin() + 2, opt_field.end()), ":");
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
                    for (int i = 9; i < tokens.size(); i++){
                         //opt fields are in key:type:val format
                        vector<string> opt_field = pliib::split(tokens[i], ':');
                        opt_elem o;
                        o.key = opt_field[0];
                        o.type = opt_field[1];
                        o.val = join(vector<string> (opt_field.begin() + 2, opt_field.end()), ":");
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
                vector<string> g_ids = pliib::split(tokens[2], ' ');
                for (int i = 0 ; i < g_ids.size(); i++){
                        g.items.push_back(g_ids[i].substr(0, g_ids[i].length() - 1));
                        g.orientations.push_back(g_ids[i].back() == '+');
                }
                if (tokens.size() > 8){
                    for (int i = 9; i < tokens.size(); i++){
                         //opt fields are in key:type:val format
                        vector<string> opt_field = pliib::split(tokens[i], ':');
                        opt_elem o;
                        o.key = opt_field[0];
                        o.type = opt_field[1];
                        o.val = join(vector<string> (opt_field.begin() + 2, opt_field.end()), ":");
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
                    for (int i = 9; i < tokens.size(); i++){
                         //opt fields are in key:type:val format
                        vector<string> opt_field = pliib::split(tokens[i], ':');
                        opt_elem o;
                        o.key = opt_field[0];
                        o.type = opt_field[1];
                        o.val = join(vector<string> (opt_field.begin() + 2, opt_field.end()), ":");
                        g.tags[o.key] = o;
           
                    }
                }
                add_group(g);
            }
            else if (tokens[0] ==  "L"){
                // TODO: we need to deal with  where the link is given before
                // its corresponding sequence in the file. TODO this is probably
                // now fixed by using the string: sequence map.
                edge_elem e;
                e.type = 1;
                e.source_name = tokens[1];
                e.sink_name = tokens[3];
                //TODO: search the input strings for "-" and "+" and set using ternary operator
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
                    for (int i = 6; i < tokens.size(); i++){
                         //opt fields are in key:type:val format
                        vector<string> opt_field = pliib::split(tokens[i], ':');
                        opt_elem o;
                        o.key = opt_field[0];
                        o.type = opt_field[1];
                        o.val = join(vector<string> (opt_field.begin() + 2, opt_field.end()), ":");
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
                    string overlap = "";
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
                    for (int i = 8; i < tokens.size(); i++){
                         //opt fields are in key:type:val format
                        vector<string> opt_field = pliib::split(tokens[i], ':');
                        opt_elem o;
                        o.key = opt_field[0];
                        o.type = opt_field[1];
                        o.val = join(vector<string> (opt_field.begin() + 2, opt_field.end()), ":");
                        e.tags[o.key] = o;
           
                    }
                }


                add_edge(e.source_name, e);
            }
            else if (tokens[0] == "W"){
                walk_elem w;
                w.source_name = tokens[1];
                w.name = tokens[2];
                //TODO check wheteher the next token is rank or direction
                if (tokens[3].compare("+") == 0 || tokens[3].compare("-") == 0){
                    w.rank = 0;
                    w.is_reverse = tokens[3] == "+" ? false : true;
                    w.cigar = tokens[4];

                }
                else{
                    w.rank = atol(tokens[3].c_str());
                    w.is_reverse = tokens[4] == "+" ? false : true;
                    w.cigar = tokens[5];
                }
                add_walk(w.source_name, w);
            }
            else if (tokens[0] == "P"){
                if (this->version >= 1.0 && this->version < 2.0){
                    // Parse a GFA 1.0 path element
                    path_elem p;
                    p.name = tokens[1];
                    vector<string> ids_and_orientations;
                    pliib::split(tokens[2], ',', ids_and_orientations);
                    for (auto x : ids_and_orientations){
                        bool orientation = ((x.back()) == '+' || x.front() == '+');
                        string id = x.substr(0, x.length() - 1);
                        p.segment_names.push_back(id);
                        p.orientations.push_back(orientation);
                    }

                    if (tokens.size() > 3){
                        vector<string> spltz = pliib::split(tokens[3], ',');
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
                    walk_elem w;
                    w.source_name = tokens[1];
                    w.name = tokens[2];
                    //TODO check wheteher the next token is rank or direction
                    if (tokens[3].compare("+") == 0 || tokens[3].compare("-") == 0){
                        w.rank = 0;
                        w.is_reverse = tokens[3] == "+" ? false : true;
                        w.cigar = tokens[4];

                    }
                    else{
                        w.rank = stol(tokens[3].c_str());
                        w.is_reverse = tokens[4] == "+" ? false : true;
                        w.cigar = tokens[5];
                    }
                    add_walk(w.source_name, w);
                }
                else{
                    cerr << "Cannot parse; version of GFA is too new. Version: " << this->version << endl;
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
                cerr << "Unknown line identifier  encountered: " << tokens[0] <<  " . Exiting." << endl;
                exit(1);
            }

        }

        if (! (this->normalized_paths)){
            walks_as_paths();
        }
        if (! ( this->normalized_walks)){
            paths_as_walks();
        }

        gfa_1_ize();
        gfa_2_ize();
        

        return true;

    }

    void GFAKluge::add_sequence(sequence_elem s){
        name_to_seq[s.name] = s;
    }

    void GFAKluge::add_path(string path_name, path_elem p){
        name_to_path[path_name] = p;
    }

    void GFAKluge::add_walk(string seq_name, walk_elem w){
        seq_to_walks[seq_name].push_back(w);
    }

    void GFAKluge::add_link(const sequence_elem& seq, const link_elem& link){
        edge_elem e(link);
        seq_to_edges[seq.name].push_back(e);
    }

    void GFAKluge::add_contained(sequence_elem seq, contained_elem con){
        seq_to_contained[seq.name].push_back(con);
    }

    void GFAKluge::add_link(const string& seq_name, const link_elem& link){
        edge_elem e(link);
        seq_to_edges[seq_name].push_back(e); 
    }

    void GFAKluge::add_alignment(string seq_name, alignment_elem a){
        seq_to_alignment[seq_name].push_back(a);
    }

    void GFAKluge::add_alignment(sequence_elem seq, alignment_elem a){
        seq_to_alignment[seq.name].push_back(a);
    }

    void GFAKluge::add_contained(string seq_name, contained_elem con){
        seq_to_contained[seq_name].push_back(con);
    }

    vector<link_elem> GFAKluge::get_links(const string& seq_name){
        return seq_to_link[seq_name];
    }

    vector<link_elem> GFAKluge::get_links(const sequence_elem& seq){
        string seq_name = seq.name;
        return seq_to_link[seq_name];
    }

    vector<contained_elem> GFAKluge::get_contained(string seq_name){
        return seq_to_contained[seq_name];
    }

    vector<contained_elem> GFAKluge::get_contained(sequence_elem seq){
        string seq_name = seq.name;
        return seq_to_contained[seq_name];
    }

    map<string, sequence_elem, custom_key> GFAKluge::get_name_to_seq(){
        return name_to_seq;
    }

    map<string, vector<link_elem> > GFAKluge::get_seq_to_link(){
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

    map<string, path_elem> GFAKluge::get_name_to_path(){
        return name_to_path;
    }

    map<string, vector<walk_elem> > GFAKluge::get_seq_to_walks(){
        return seq_to_walks;
    }

    map<string, vector<contained_elem> > GFAKluge::get_seq_to_contained(){
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

    map<string, vector<alignment_elem> > GFAKluge::get_seq_to_alignment(){
        return seq_to_alignment;
    }

    map<string, header_elem> GFAKluge::get_header(){
        return header;
    }

    string GFAKluge::join(const vector<string>& splits, const string& glue){
        string ret = "";
        for (int i = 0; i < splits.size(); i++){
            if (i != 0){
                ret += glue;
            }
            ret += splits[i];
        }

        return ret;
    }


    //TODO we should use a string stream here for efficiency.
    string GFAKluge::header_string(map<string, header_elem>& headers){
        string ret = "H";
        map<string, header_elem>::iterator it;
        for (it = headers.begin(); it != headers.end(); it++){
            ret += "\t";
            header_elem h = it->second;
            string t[] = {h.key, h.type, h.val};
            vector<string> temp = vector<string> (t, t + sizeof(t) / sizeof(string));
            ret += join(temp, ":");
        }
        return ret;

    }

    void GFAKluge::set_walks(bool ws){
        this->use_walks = ws;
    }

    void GFAKluge::paths_as_walks(){
        if (!this->normalized_paths && seq_to_walks.empty() && !name_to_path.empty()){
            if (name_to_path.empty() && !this->normalized_walks){
                walks_as_paths();
            }
        
            map<string, path_elem>::iterator it;
            for (it = name_to_path.begin(); it != name_to_path.end(); it++){
                // Create a Walk for every segment in the path's segment list
                for (int i = 0; i < it->second.segment_names.size(); i++){
                    walk_elem w;
                    w.name = it->second.name;
                    w.rank = i;
                    w.source_name = it->second.segment_names[i];
                    w.is_reverse = (it->second.orientations[i]);
                    for (int z = 0; z < it->second.overlaps.size(); z++){
                        w.cigar = it->second.overlaps[z];
                    }
                    for (int z = it->second.overlaps.size(); z < it->second.segment_names.size(); z++ ){
                        w.cigar = "*";
                    }
                    add_walk(w.source_name, w);
                }
            }
            this->normalized_paths = true;
            }
        
    }

    void GFAKluge::walks_as_paths(){
        if (!this->normalized_walks && name_to_path.empty() && !seq_to_walks.empty()){
            if (seq_to_walks.empty() && !this->normalized_paths){
                paths_as_walks();
            }
            map<string, vector<walk_elem> > pathname_to_walk_elems;
            for(auto n_to_s : name_to_seq){
                for (auto w : seq_to_walks[n_to_s.first]){
                    pathname_to_walk_elems[w.name].push_back(w);
                }   
            }

            struct walk_sort_key{
                inline bool operator() (walk_elem first, walk_elem second){
                    return (first.rank != 0 && second.rank != 0) ? (first.rank < second.rank) : false;
                }
            };
            

            for (auto x : pathname_to_walk_elems){
                path_elem p;
                p.name = x.first;
                //cout << x.first << "\t";
                std::sort(x.second.begin(), x.second.end(), walk_sort_key());
                for (auto y : x.second){
                    
                //  cout << y.source_name << "\t";
                    p.segment_names.push_back(y.source_name);
                    p.overlaps.push_back(y.cigar);
                    p.orientations.push_back(y.is_reverse);
                }
                //cout << endl;
                add_path(p.name, p);
            }
        
            this->normalized_walks = true;
        }
            

    }

    string GFAKluge::opt_string(vector<opt_elem> opts){
        string ret = "";
        for (int i = 0; i < opts.size(); i++){
            opt_elem o = opts[i];
            if (i != 0){
                ret += "\t";
            }
            string t [] = {o.key, o.type, o.val};
            vector<string> temp = vector<string> (t, t + sizeof(t) / sizeof(string));
            ret += join(temp, ":");
        }
        return ret;
    }

    std::string GFAKluge::block_order_string(){
        
        this->gfa_1_ize();
        this->gfa_2_ize();

        if (version >= 2.0){
            return block_order_string_2();
        }
			stringstream ret;
			int i;
			//First print header lines.
			if (header.size() > 0){
					ret << header_string(header) + "\n";
			}
            if (!normalized_paths){
                walks_as_paths();
            }
            if (!normalized_walks){
                paths_as_walks();
            }
			map<std::string, sequence_elem>::iterator st;

			for (st = name_to_seq.begin(); st != name_to_seq.end(); st++){
				ret << st->second.to_string_1() << endl;
			}



			for (st = name_to_seq.begin(); st != name_to_seq.end(); st++){
                for (auto e : seq_to_edges[st->first]){
                    if (e.type == 1){
                        ret << e.to_string_1() << endl;
                    }

                }
            }
            for (st = name_to_seq.begin(); st != name_to_seq.end(); st++){
                for (auto e : seq_to_edges[st->first]){
                    if (e.type == 2){
                        ret << e.to_string_1() << endl;
                    }

                }
            }
            
			for (st = name_to_seq.begin(); st != name_to_seq.end(); st++){
                if (this->version < 1.0 && seq_to_walks[st->first].size() > 0){
                    for (i = 0; i < seq_to_walks[st->first].size(); i++){
								stringstream pat;
								pat << (use_walks ? "W" : "P") << "\t" + seq_to_walks[st->first][i].source_name << "\t";
								pat << seq_to_walks[st->first][i].name << "\t";
								if (!(seq_to_walks[st->first][i].rank ==  0L)){
										pat << seq_to_walks[st->first][i].rank << "\t";
								}
								pat << (seq_to_walks[st->first][i].is_reverse ? "-" : "+");
								pat << "\t";
								pat << seq_to_walks[st->first][i].cigar + "\n";
								ret << pat.str();
					}
                }
			}
            if (name_to_path.size() > 0 && this->version == 1.0){
                map<string, path_elem>::iterator pt;
                for (pt = name_to_path.begin(); pt != name_to_path.end(); ++pt){
                    stringstream pat;
                    pat << "P\t" << pt->second.name << "\t";
                    vector<string> ovec;
                    for (int oi = 0; oi < pt->second.segment_names.size(); oi++){
                        stringstream o_str;
                        o_str << pt->second.segment_names[oi] << (pt->second.orientations[oi] ? "-" : "+");
                        ovec.push_back(o_str.str());
                    }
                    pat << join(ovec, ",");
                    if (pt->second.overlaps.size() > 0){
                        pat << "\t" << join(pt->second.overlaps, ",");
                    }
                    pat << "\n";
                    ret << pat.str();
                }
            }
			return ret.str();

    }

    std::string GFAKluge::to_string_2(){
        this->gfa_2_ize();
        
        stringstream ret;
        // Header
        if (header.size() > 0){
            ret << header_string(header) << endl;
        }
        for (auto p : groups){
            ret << p.second.to_string_2() << endl;
        }
        for (auto s : name_to_seq){
            ret << s.second.to_string_2() << endl;
            for (auto f : seq_to_fragments[s.first]){
                ret << f.to_string_2() << endl;
            }
            for (auto e : seq_to_edges[s.first]){
                ret << e.to_string_2() << endl;
            }
            for (auto g : seq_to_gaps[s.first]){
                ret << g.to_string_2() << endl;
            }
        }
        return ret.str();
    }

    std::string GFAKluge::block_order_string_2(){
        this->gfa_2_ize();

        stringstream ret;
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


    std::string GFAKluge::to_string(){

        gfa_1_ize();
        gfa_2_ize();
        if (this->version >= 2.0){
            return to_string_2();
        }

        stringstream ret;
        int i;
        //First print header lines.
        if (header.size() > 0){
            ret << header_string(header) + "\n";
        }
        
        if (!this->normalized_walks){
            walks_as_paths();
        }
        if (!this->normalized_paths){
            paths_as_walks();
        }
        if (name_to_path.size() > 0 && this->version >= 1.0){
                map<string, path_elem>::iterator pt;
                for (pt = name_to_path.begin(); pt != name_to_path.end(); ++pt){
                    stringstream pat;
                    pat << "P" << "\t";
                    pat << pt->second.name << "\t";
                    vector<string> ovec;
                    for (int seg_ind = 0; seg_ind < pt->second.segment_names.size(); seg_ind++){
                        stringstream o_str;
                        o_str << pt->second.segment_names[seg_ind] << (pt->second.orientations[seg_ind] ? "+" : "-");
                        ovec.push_back(o_str.str());
                    }
                    pat << join(ovec, ",");
                    if (pt->second.overlaps.size() > 0){
                        pat << "\t" << join(pt->second.overlaps, ",");
                    }
                    pat << "\n";
                    ret << pat.str();
                }
        }
        for (auto s : name_to_seq){
            ret << s.second.to_string_1() << endl;
            for (auto e : seq_to_edges[s.first]){
                ret << e.to_string_1() << endl;;
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

    tuple<uint64_t, uint64_t, uint64_t, uint64_t, uint64_t> GFAKluge::max_ids(){
        return std::make_tuple(this->base_seq_id, this->base_edge_id,
             this->base_frag_id, this->base_gap_id, this->base_group_id);
    }

    std::string GFAKluge::max_ids_string(){
        tuple<uint64_t, uint64_t, uint64_t, uint64_t, uint64_t> x = max_ids();
        std::vector<string> max_str(5);
        max_str[0] = std::to_string(std::get<0>(x));
        max_str[1] = std::to_string(std::get<1>(x));
        max_str[2] = std::to_string(std::get<2>(x));
        max_str[3] = std::to_string(std::get<3>(x));
        max_str[4] = std::to_string(std::get<4>(x));
        return join(max_str, ":");
    }

    void GFAKluge::re_id(std::string new_mx_str){
        vector<uint64_t> starts(5);
        vector<string> starts_strs = pliib::split(new_mx_str, ':');
        for (int i = 0; i < starts_strs.size(); ++i){
            starts[i] = stoul(starts_strs[i]);
        }
        tuple<uint64_t, uint64_t, uint64_t, uint64_t, uint64_t> n_ids = std::make_tuple(starts[0], starts[1],
            starts[2], starts[3], starts[4]);
        re_id(n_ids);
        
    }
    void GFAKluge::re_id(tuple<uint64_t, uint64_t,
                               uint64_t, uint64_t,
                               uint64_t>& new_mx){

        base_seq_id = std::get<0>(new_mx);
        base_edge_id = std::get<1>(new_mx);
        base_frag_id = std::get<2>(new_mx);
        base_gap_id = std::get<3>(new_mx);
        base_group_id = std::get<4>(new_mx);
        uint64_t seg_diff = base_seq_id;
        // Segments
        int num_segs = name_to_seq.size();

        // locally cache name_to_seq,
        // seq_to_edges, seq_to_fragments, groups,
        // and seq_to_gaps before clearing them.
        map<string, sequence_elem, custom_key> n_s;
        map<string, vector<edge_elem>> s_e;
        map<string, vector<fragment_elem>> s_f;
        map<string, vector<gap_elem>> s_g;
        map<string, group_elem> g_g;



        for (auto ns : name_to_seq){
            string old_name = ns.second.name;
            ns.second.id = ++base_seq_id;
            ns.second.name = std::to_string(ns.second.id);    

            if (one_compat){
                // Links

                // Contains

                // Paths
            }
        
            if (two_compat && version >= 2.0){
                uint64_t edge_count = 0;
                // Edges
                for (auto e : seq_to_edges[old_name]){
                    e.source_name = ns.second.name;
                    e.sink_name += std::to_string(seg_diff + stoul(e.sink_name));
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
                    for (int i = 0; i < g.second.items.size(); ++i){
                        g.second.items[i] = std::to_string(seg_diff + stoul(g.second.items[i]));
                    }
                    g_g[g.first] = g.second;
                    //add_group(g.second);
                }
                
                //add_sequence(ns.second);
        }
        n_s[ns.second.name] = ns.second;

    }

        // Clear our internal maps
        name_to_seq = n_s;
        seq_to_edges = s_e;
        seq_to_fragments = s_f;
        seq_to_gaps = s_g;
        groups = g_g;
        
    }

    // TODO check incompatible version numbers
    // TODO Check colliding groups, headers
    void GFAKluge::merge(GFAKluge& gg){
        std::unordered_set<string> seg_ids;
        // Merge headers
        for (auto h : gg.get_header()){
            header[h.first] = h.second;
        }
        for (auto s : this->get_name_to_seq()){
            seg_ids.insert(s.first);
        }

        map<string, sequence_elem, custom_key> ss = gg.get_name_to_seq();
        map<string, vector<fragment_elem>> sf = gg.get_seq_to_fragments();
        map<string, vector<gap_elem>> sg = gg.get_seq_to_gaps();
        map<string, vector<edge_elem>> se = gg.get_seq_to_edges();

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
                    cerr << "WARNING: DUPLICATE IDS " << s.second.name << endl <<
                    " will be lost." << endl;
                }
                for (auto g : gg.get_groups()){
                    this->add_group(g.second);
                }
            }
        }
        else if (one_compat){

        }
        
    }

    // we calculate the N50 based on the 'S' lines,
    // Though in theory an O line might also be a contig
    double GFAKluge::get_N50(){
        vector<double> s_lens;
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
        for (int i = 0; i < s_lens.size(); i++){
            cumul_size += s_lens[i];
            if (cumul_size >= avg){
                n = s_lens[i];
                break;
            }
        }
        return n;
    }

    double GFAKluge::get_N90(){
        vector<double> s_lens;
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
        for (int i = 0; i < s_lens.size(); i++){
            cumul_size += s_lens[i];
            if (cumul_size >= avg){
                n = s_lens[i];
                break;
            }
        }
        return n;
    }
    

    int GFAKluge::get_L50(){
        vector<double> s_lens;
        uint64_t total_len = 0;
        double n = 0.0;
        for (auto s = name_to_seq.begin(); s != name_to_seq.end(); s++){
            s_lens.push_back(s->second.length);
            total_len += s->second.length;
        }
        double avg = total_len * 0.50;
        std::sort(s_lens.rbegin(), s_lens.rend());
        //int middle = floor((double) s_lens.size() / 2.0);
        double cumul_size = 0.0;
        for (int i = 0; i < s_lens.size(); i++){
            cumul_size += s_lens[i];
            if (cumul_size >= avg){
                return i + 1;
            }
        }
        return -1;
    }

    int GFAKluge::get_L90(){
        vector<double> s_lens;
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
        for (int i = 0; i < s_lens.size(); i++){
            cumul_size += s_lens[i];
            if (cumul_size >= avg){
                return i + 1;
            }
        }
        return -1;
    }

    void GFAKluge::output_to_stream(std::ostream& os, bool output_block_order){
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
                os << header_string(header) << endl;
            }
            for (auto p : groups){
                os << p.second.to_string_2() << endl;
            }
            for (auto s : name_to_seq){
                os << s.second.to_string_2() << endl;
                for (auto f : seq_to_fragments[s.first]){
                    os << f.to_string_2() << endl;
                }
                for (auto e : seq_to_edges[s.first]){
                    os << e.to_string_2() << endl;
                }
                for (auto g : seq_to_gaps[s.first]){
                    os << g.to_string_2() << endl;
                }
            }
        }
        else if (this->version < 2.0 && output_block_order){
			//First print header lines.
			if (header.size() > 0){
					os << header_string(header) + "\n";
			}
			map<std::string, sequence_elem>::iterator st;

			for (st = name_to_seq.begin(); st != name_to_seq.end(); st++){
				os << st->second.to_string_1() << endl;
			}


            if (!this->normalized_walks){
                walks_as_paths();
            }
            if (!this->normalized_paths){
                paths_as_walks();
            }

			for (st = name_to_seq.begin(); st != name_to_seq.end(); st++){
                for (auto e : seq_to_edges[st->first]){
                    if (e.type == 1){
                        os << e.to_string_1() << endl;
                    }

                }
            }
            for (st = name_to_seq.begin(); st != name_to_seq.end(); st++){
                for (auto e : seq_to_edges[st->first]){
                    if (e.type == 2){
                        os << e.to_string_1() << endl;
                    }

                }
            }
            
			for (st = name_to_seq.begin(); st != name_to_seq.end(); st++){
                if (this->version < 1.0 && seq_to_walks[st->first].size() > 0){
                    for (int i = 0; i < seq_to_walks[st->first].size(); i++){
								stringstream pat;
								pat << (use_walks ? "W" : "P") << "\t" + seq_to_walks[st->first][i].source_name << "\t";
								pat << seq_to_walks[st->first][i].name << "\t";
								if (!(seq_to_walks[st->first][i].rank ==  0L)){
										pat << seq_to_walks[st->first][i].rank << "\t";
								}
								pat << (seq_to_walks[st->first][i].is_reverse ? "-" : "+");
								pat << "\t";
								pat << seq_to_walks[st->first][i].cigar + "\n";
								os << pat.str();
					}
                }
			}
            if (name_to_path.size() > 0 && this->version == 1.0){
                map<string, path_elem>::iterator pt;
                for (pt = name_to_path.begin(); pt != name_to_path.end(); ++pt){
                    stringstream pat;
                    pat << "P\t" << pt->second.name << "\t";
                    vector<string> ovec;
                    for (int oi = 0; oi < pt->second.segment_names.size(); oi++){
                        stringstream o_str;
                        o_str << pt->second.segment_names[oi] << (pt->second.orientations[oi] ? "-" : "+");
                        ovec.push_back(o_str.str());
                    }
                    pat << join(ovec, ",");
                    if (pt->second.overlaps.size() > 0){
                        pat << "\t" << join(pt->second.overlaps, ",");
                    }
                    pat << "\n";
                    os << pat.str();
                }
            }


        }
        else if (this->version < 2.0){
            //First print header lines.
            if (header.size() > 0){
                os << header_string(header) + "\n";
            }
        
            if (!this->normalized_walks){
                walks_as_paths();
            }
            if (!this->normalized_paths){
                paths_as_walks();
            }
            if (name_to_path.size() > 0 && this->version >= 1.0){
                    map<string, path_elem>::iterator pt;
                    for (pt = name_to_path.begin(); pt != name_to_path.end(); ++pt){
                        stringstream pat;
                        pat << "P" << "\t";
                        pat << pt->second.name << "\t";
                        vector<string> ovec;
                        for (int seg_ind = 0; seg_ind < pt->second.segment_names.size(); seg_ind++){
                            stringstream o_str;
                            o_str << pt->second.segment_names[seg_ind] << (pt->second.orientations[seg_ind] ? "+" : "-");
                            ovec.push_back(o_str.str());
                        }
                        pat << join(ovec, ",");
                        if (pt->second.overlaps.size() > 0){
                            pat << "\t" << join(pt->second.overlaps, ",");
                        }
                        pat << "\n";
                        os << pat.str();
                    }
            }

			map<std::string, sequence_elem>::iterator st;
            for (st = name_to_seq.begin(); st != name_to_seq.end(); st++){
                if (this->version < 1.0 && seq_to_walks[st->first].size() > 0){
                    for (int i = 0; i < seq_to_walks[st->first].size(); i++){
								stringstream pat;
								pat << (use_walks ? "W" : "P") << "\t" + seq_to_walks[st->first][i].source_name << "\t";
								pat << seq_to_walks[st->first][i].name << "\t";
								if (!(seq_to_walks[st->first][i].rank ==  0L)){
										pat << seq_to_walks[st->first][i].rank << "\t";
								}
								pat << (seq_to_walks[st->first][i].is_reverse ? "-" : "+");
								pat << "\t";
								pat << seq_to_walks[st->first][i].cigar + "\n";
								os << pat.str();
					}
                }
			}
            for (auto s : name_to_seq){
                os << s.second.to_string_1() << endl;
                for (auto e : seq_to_edges[s.first]){
                    os << e.to_string_1() << endl;;
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
            cerr << "Invalid version " << this->version << endl;
                exit(9);

        }
    }

    void GFAKluge::write_element(std::ostream& os, const sequence_elem& s){
        if (this->version >= 2.0){
            os << s.to_string_2();
        }
        else{
            os << s.to_string_1();
        }
    }
    void GFAKluge::write_element(std::ostream& os, edge_elem e){
        if (this->version >= 2.0){
            os << e.to_string_2();
        }
        else{
            os << e.to_string_1();
        }
    }
    void GFAKluge::write_element(std::ostream& os, const fragment_elem& f){
        os << f.to_string();
    }
    void GFAKluge::write_element(std::ostream& os, const group_elem& g){
        if (this->version >= 2.0){
            os << g.to_string_2();
        }
        else{
            os << g.to_string_1();
        }
    }
    void GFAKluge::write_element(std::ostream& os, const gap_elem& g){
        os << g.to_string();
    }

    // Avoid calling to_string as it murders mem usage
    std::ostream& operator<<(std::ostream& os, GFAKluge& g){
        g.gfa_1_ize();
        g.gfa_2_ize();

        //os << g.to_string();
        g.output_to_stream(os);
        return os;
    }
}
