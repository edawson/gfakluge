#include "gfakluge.hpp"

using namespace std;
namespace gfak{


    void GFAKluge::add_edge(string seqname, edge_elem e){
        seq_to_edges[seqname].push_back(e);
    }
    void GFAKluge::add_edge(sequence_elem s, edge_elem e){
        add_edge(s.name, e);
    }

    void GFAKluge::add_fragment(string seqname, fragment_elem f){
        seq_to_fragments[seqname].push_back(f);
    }

    void GFAKluge::add_fragment(sequence_elem s, fragment_elem f){
        add_fragment(s.name, f);
    }

    void GFAKluge::add_gap(gap_elem g){
        seq_to_gaps[g.source_name].push_back(g);
    }

    void GFAKluge::add_group(group_elem g){
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

    // MISNOMER: also converts edges -> containments
    void GFAKluge::edges_as_links(){
        for (auto s : name_to_seq){
            for (auto e : seq_to_edges[s.first]){
                int t = e.determine_type();
                if (t == 1){
                    link_elem l;
                    l.source_name = e.source_name;
                    l.sink_name = e.sink_name;
                    l.source_orientation_forward = e.source_orientation_forward;
                    l.sink_orientation_forward = e.sink_orientation_forward;
                    l.cigar = e.alignment;
                    l.opt_fields = e.tags;
                    add_link(l.source_name, l);
                }
                else if (t == 2){
                    contained_elem c;
                    c.source_name = e.source_name;
                    c.sink_name = e.sink_name;
                    c.source_orientation_forward = e.source_orientation_forward;
                    c.sink_orientation_forward = e.sink_orientation_forward;
                    c.pos = e.source_begin;
                    c.cigar = e.alignment;
                    c.opt_fields = e.tags;
                    add_contained(c.source_name, c);
                }
                else{
                    cerr << "Skipping edge not expressable in GFA2: \"" << e.to_string_2() << "\"" << endl;
                }
            }
        }
    }

    void GFAKluge::gfa_1_ize(){
        /**
        * Swap edges to links
        * Ordered groups -> paths,
        * warn about missing sets
        * warn about missing fragments (we could output them as seqs)
        * warn about missing gaps (we could remove them)
        * 
        */
        groups_as_paths();
        edges_as_links();
    }

    void GFAKluge::gfa_2_ize(){
        if (!two_compat){
            // Fix S line length field if needed.
            for (auto s : name_to_seq){
                s.second.length = (uint64_t) s.second.sequence.length();
            }
            // L and C lines -> E lines
            for (auto s : seq_to_link){
                sequence_elem s_seq = name_to_seq[s.first];
                for (auto l : s.second){
                    edge_elem e;

                    e.id = std::to_string(++base_edge_id);
                    e.type.set(0, 1);
                    e.source_name = l.source_name;
                    e.source_begin = (l.source_orientation_forward ? s_seq.length : 0);
                    e.source_end = (l.source_orientation_forward ? s_seq.length : 0);
                    e.sink_name = l.sink_name;
                    sequence_elem sink_seq = name_to_seq[l.sink_name];
                    e.sink_begin = (l.sink_orientation_forward ? 0 : sink_seq.length);
                    e.sink_end = (l.sink_orientation_forward ? 0 : sink_seq.length);

                    e.source_orientation_forward = l.source_orientation_forward;
                    e.sink_orientation_forward = l.sink_orientation_forward;
                    e.alignment = l.cigar;
                    e.tags = l.opt_fields;
                    if (l.source_orientation_forward){
                        e.ends.set(0,1);
                        e.ends.set(1,1);
                    }
                    if (!l.sink_orientation_forward){
                        e.ends.set(2,1);
                        e.ends.set(3,1);
                    }
                    
                    
                    
                    add_edge(e.source_name, e);
                }
            }
            for (auto s : seq_to_contained){
                for (auto c : s.second){
                    edge_elem e;
                    e.id = std::to_string(++base_edge_id);
                    e.type.set(1,1);
                    e.source_name = c.source_name;
                    e.sink_name = c.sink_name;
                    e.source_orientation_forward = c.source_orientation_forward;
                    e.sink_orientation_forward = c.sink_orientation_forward;
                    e.source_begin = c.pos;
                    e.source_end = c.pos + name_to_seq[c.sink_name].length;
                    if (c.pos == name_to_seq[c.source_name].length){
                        e.ends.set(0, 1);
                        e.ends.set(1, 1);
                    }
                    e.sink_begin = 0;
                    e.sink_end = name_to_seq[c.sink_name].length;
                    e.ends.set(2, (c.sink_orientation_forward ? 0 : 1));
                    e.ends.set(3, (c.sink_orientation_forward ? 1 : 0));
                
                    e.alignment = c.cigar;
                    e.tags = c.opt_fields;
                    add_edge(e.source_name, e);
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
    }
    void GFAKluge::set_version(){
        header_elem verz;
        verz.key = "VN";
        verz.type="Z";
        verz.val = std::to_string((double) this->version).substr(0,3);

        this->header[verz.key] = verz;
    }

    // Borrow from
    //http://stackoverflow.com/questions/236129/split-a-string-in-c
    // Thanks StackOverflow!
    vector<string> GFAKluge::split(string s, char delim){
        vector<string> ret;
        stringstream sstream(s);
        string temp;
        while(getline(sstream, temp, delim)){
            ret.push_back(temp);
        }
        return ret;

    }

    bool GFAKluge::parse_gfa_file(string filename){
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
            vector<string> tokens = split(line, '\t');
            if (tokens[0] == "H"){
                header_elem h;
                line_tokens = split(tokens[1], ':');
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
                if (this->version >= 2.0){
                   s.length = stoi(tokens[2]); 
                   s.sequence = tokens[3];
                   tag_index = 4;
                }
                else{
                    s.sequence = tokens[2];
                    s.length = s.sequence.length();
                    tag_index = 3;
                }
                //s.id = atol(s.name.c_str());
                int i;
                if (tokens.size() > 3){
                    for (i = tag_index; i < tokens.size(); i++){
                        //opt fields are in key:type:val format
                        vector<string> opt_field = split(tokens[i], ':');
                        opt_elem o;
                        o.key = opt_field[0];
                        o.type = opt_field[1];
                        o.val = join(vector<string> (opt_field.begin() + 2, opt_field.end()), ":");
                        s.opt_fields.push_back(o);
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
                e.ends.set(0, (x.back() == '+' ? 1 : 0));
                e.source_begin = (e.ends.test(0) ? stoul(x.substr(0, x.length() - 1)) : stoul(x));

                x = tokens[5];
                e.ends.set(1, (x.back() == '+' ? 1 : 0));
                e.source_end = (e.ends.test(1) ? stoul(x.substr(0, x.length() - 1)) : stoul(x));

                x = tokens[6];
                e.ends.set(2, (x.back() == '+' ? 1 : 0));
                e.sink_begin = (e.ends.test(2) ? stoul(x.substr(0, x.length() - 1)) : stoul(x));


                x = tokens[7];
                e.ends.set(3, (x.back() == '+' ? 1 : 0));
                e.sink_end = (e.ends.test(3) ? stoul(x.substr(0, x.length() - 1)) : stoul(x));
                
                e.alignment = tokens[8];

                if (tokens.size() > 9){
                    for (int i = 9; i < tokens.size(); i++){
                         //opt fields are in key:type:val format
                        vector<string> opt_field = split(tokens[i], ':');
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
                f.alignment = tokens[7];
                if (tokens.size() > 8){
                    for (int i = 9; i < tokens.size(); i++){
                         //opt fields are in key:type:val format
                        vector<string> opt_field = split(tokens[i], ':');
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
                vector<string> g_ids = split(tokens[2], ' ');
                for (int i = 0 ; i < g_ids.size(); i++){
                        g.items.push_back(g_ids[i].substr(0, g_ids[i].length() - 1));
                        g.orientations.push_back(g_ids[i].back() == '+');
                }
                if (tokens.size() > 8){
                    for (int i = 9; i < tokens.size(); i++){
                         //opt fields are in key:type:val format
                        vector<string> opt_field = split(tokens[i], ':');
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
                g.items = split(tokens[2], ' ');
                if (tokens.size() > 8){
                    for (int i = 9; i < tokens.size(); i++){
                         //opt fields are in key:type:val format
                        vector<string> opt_field = split(tokens[i], ':');
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
                link_elem l;
                l.source_name = tokens[1];
                l.sink_name = tokens[3];
                //TODO: search the input strings for "-" and "+" and set using ternary operator
                l.source_orientation_forward = tokens[2] == "+" ? true : false;
                l.sink_orientation_forward = tokens[4] == "+" ? true : false;
                //l.pos = tokens[0];
                if (tokens.size() >= 6){
                    l.cigar = tokens[5];
                }
                else{
                    l.cigar = "*";
                }
                add_link(l.source_name, l);
            }
            else if (tokens[0] == "C"){
                contained_elem c;
                //TODO fix token indices here
                c.source_name = tokens[1];
                c.sink_name = tokens[3];
                c.source_orientation_forward = tokens[2] == "+" ? true : false;
                c.sink_orientation_forward = tokens[4] == "+" ? true : false;
                c.pos = atoi(tokens[5].c_str());
                if (tokens.size() > 6){
                    c.cigar = tokens[6];
                }
                else{
                    c.cigar = "*";
                }

                add_contained(c.source_name, c);
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
                    vector<string> ids_and_orientations = split(tokens[2], ',');
                    for (auto x : ids_and_orientations){
                        bool orientation = ((x.back()) == '+' || x.front() == '+');
                        string id = x.substr(0, x.length() - 1);
                        p.segment_names.push_back(id);
                        p.orientations.push_back(orientation);
                    }

                    if (tokens.size() > 3){
                        vector<string> spltz = split(tokens[3], ',');
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

    void GFAKluge::add_link(sequence_elem seq, link_elem link){
        seq_to_link[seq.name].push_back(link);
    }

    void GFAKluge::add_contained(sequence_elem seq, contained_elem con){
        seq_to_contained[seq.name].push_back(con);
    }

    void GFAKluge::add_link(string seq_name, link_elem link){
        seq_to_link[seq_name].push_back(link);
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

    vector<link_elem> GFAKluge::get_links(string seq_name){
        return seq_to_link[seq_name];
    }

    vector<link_elem> GFAKluge::get_links(sequence_elem seq){
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
        return seq_to_link;
    }

    map<string, path_elem> GFAKluge::get_name_to_path(){
        return name_to_path;
    }

    map<string, vector<walk_elem> > GFAKluge::get_seq_to_walks(){
        return seq_to_walks;
    }

    map<string, vector<contained_elem> > GFAKluge::get_seq_to_contained(){
        return seq_to_contained;
    }

    map<string, vector<alignment_elem> > GFAKluge::get_seq_to_alignment(){
        return seq_to_alignment;
    }

    map<string, header_elem> GFAKluge::get_header(){
        return header;
    }

    string GFAKluge::join(vector<string> splits, string glue){
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
			if (name_to_seq.size() > 0){
					for (st = name_to_seq.begin(); st != name_to_seq.end(); st++){
							ret << "S\t" + (st->second).name + "\t" + (st->second).sequence;
							if ((st->second).opt_fields.size() > 0){
									ret << "\t" + opt_string((st->second).opt_fields);
							}
							ret << "\n";
					}
				}


						for (st = name_to_seq.begin(); st != name_to_seq.end(); st++){
							if (seq_to_link[st->first].size() > 0){
									for (i = 0; i < seq_to_link[st->first].size(); i++){
											string link = "L\t" + seq_to_link[st->first][i].source_name + "\t";
											link += (seq_to_link[st->first][i].source_orientation_forward ? "+" : "-");
											link += "\t";
											link += seq_to_link[st->first][i].sink_name + "\t";
											link += (seq_to_link[st->first][i].sink_orientation_forward ? "+" : "-");
											link += "\t";
											link += seq_to_link[st->first][i].cigar + "\n";
											ret << link;
									}

							}
						}

							//TODO iterate over contained segments
							for (st = name_to_seq.begin(); st != name_to_seq.end(); st++){
							if (seq_to_contained[st->first].size() > 0){
									for (i = 0; i < seq_to_contained[st->first].size(); i++){
											stringstream cont;
											cont <<  "C" << "\t" << seq_to_contained[st->first][i].source_name << "\t";
											cont << (seq_to_contained[st->first][i].source_orientation_forward ? "+" : "-");
											cont << "\t";
											cont << seq_to_contained[st->first][i].sink_name << "\t";
											cont << (seq_to_contained[st->first][i].sink_orientation_forward ? "+" : "-");
											cont << "\t";
											cont << seq_to_contained[st->first][i].pos << "\t";
											cont << seq_to_contained[st->first][i].cigar << "\n";
											ret << cont.str();
									}
							}
					}

					//TODO iterate over links
					//L    segName1,segOri1,segName2,segOri2,CIGAR      Link
            
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
            if (name_to_path.size() > 0 && this->version >= 1.0){
                map<string, path_elem>::iterator pt;
                for (pt = name_to_path.begin(); pt != name_to_path.end(); ++pt){
                    stringstream pat;
                    pat << "P\t" << pt->second.name << "\t";
                    vector<string> ovec;
                    for (int oi = 0; oi < pt->second.segment_names.size(); oi++){
                        stringstream o_str;
                        o_str << pt->second.segment_names[oi] << (pt->second.orientations[i] ? "-" : "+");
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



			//TODO iterate over annotation lines.


			//Print sequences and links in order, then annotation lines.

			return ret.str();

    }

    std::string GFAKluge::to_string_2(){
        this->gfa_2_ize();
        
        stringstream ret;
        // Header
        if (header.size() > 0){
            ret << header_string(header) + "\n";
        }
        for (auto p : groups){
            ret << p.second.to_string_2() << endl;
        }
        for (auto s : name_to_seq){
            ret << s.second.to_string_2() << "\n";
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


    //TODO this should use stringstream too...
    std::string GFAKluge::to_string(){
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
                        o_str << pt->second.segment_names[seg_ind] << (pt->second.orientations[seg_ind] ? "-" : "+");
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
        if (name_to_seq.size() > 0){
            map<std::string, sequence_elem>::iterator st;
            for (st = name_to_seq.begin(); st != name_to_seq.end(); st++){
                ret << "S\t" + (st->second).name + "\t" + (st->second).sequence;
                if ((st->second).opt_fields.size() > 0){
                    ret << "\t" + opt_string((st->second).opt_fields);
                }
                ret << "\n";

                if (seq_to_walks[st->first].size() > 0 && this->version < 1.0){
                    for (i = 0; i < seq_to_walks[st->first].size(); i++){
                        stringstream pat;
                        pat << (use_walks ? "W" : "P")<< "\t" + seq_to_walks[st->first][i].source_name << "\t";
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
                if (seq_to_link[st->first].size() > 0){
                    for (i = 0; i < seq_to_link[st->first].size(); i++){
                        string link = "L\t" + seq_to_link[st->first][i].source_name + "\t";
                        link += (seq_to_link[st->first][i].source_orientation_forward ? "+" : "-");
                        link += "\t";
                        link += seq_to_link[st->first][i].sink_name + "\t";
                        link += (seq_to_link[st->first][i].sink_orientation_forward ? "+" : "-");
                        link += "\t";
                        link += seq_to_link[st->first][i].cigar + "\n";
                        ret << link;
                    }

                }

                //TODO iterate over contained segments
                if (seq_to_contained[st->first].size() > 0){
                    for (i = 0; i < seq_to_contained[st->first].size(); i++){
                        stringstream cont;
                        cont <<  "C" << "\t" << seq_to_contained[st->first][i].source_name << "\t";
                        cont << (seq_to_contained[st->first][i].source_orientation_forward ? "+" : "-");
                        cont << "\t";
                        cont << seq_to_contained[st->first][i].sink_name << "\t";
                        cont << (seq_to_contained[st->first][i].sink_orientation_forward ? "+" : "-");
                        cont << "\t";
                        cont << seq_to_contained[st->first][i].pos << "\t";
                        cont << seq_to_contained[st->first][i].cigar << "\n";
                        ret << cont.str();
                    }
                }
            }


        }
        //TODO iterate over annotation lines.


        //Print sequences and links in order, then annotation lines.


        return ret.str();
    }

    tuple<uint64_t, uint64_t, uint64_t, uint64_t> GFAKluge::max_ids(){
        return std::make_tuple(this->base_seq_id, this->base_edge_id,
             this->base_frag_id, this->base_gap_id);
    }
    void GFAKluge::re_id(tuple<uint64_t, uint64_t, uint64_t, uint64_t>& new_mx){

        if (one_compat){
            // Segments
            // Links
            // Contains
            // Paths
        }
        
        if (two_compat){
            // Edges
            // Fragments
            // Gaps
            // Groups
        }
        
    }

    std::ostream& operator<<(std::ostream& os, GFAKluge& g){
        os << g.to_string();
        return os;
    }
}
