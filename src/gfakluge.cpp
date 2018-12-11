#include "gfakluge.hpp"

#include <unordered_set>
#include <fstream>

using namespace std;
namespace gfak{




    void GFAKluge::fill_sequences(const char* fasta_file){

        TFA::tiny_faidx_t tf;
        if (TFA::checkFAIndexFileExists(fasta_file)){
            TFA::parseFAIndex(fasta_file, tf);
        }
        else{
            cerr << "Creating index for " << fasta_file << "." << endl;
            TFA::createFAIndex(fasta_file, tf);
            cerr << "Created index." << endl;
            TFA::writeFAIndex(fasta_file, tf);
            cerr << "Wrote index to file." << endl;
        }
        
        for (std::map<string, sequence_elem, custom_key>::iterator it = name_to_seq.begin(); it != name_to_seq.end(); it++){
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

    bool GFAKluge::trim_seqs(const int& min_len, const bool& no_amb){
        bool graph_modified = false;

        unordered_set<string> dropped_seqs;
        map<string, sequence_elem, custom_key>::iterator n_to_s;
        for (n_to_s = name_to_seq.begin(); n_to_s != name_to_seq.end(); n_to_s++){
            auto& s = n_to_s->second;
            if (s.length == UINT64_MAX){
                cerr << "Length unset for sequence " << s.name << "; removing from graph." << endl;
                dropped_seqs.insert(s.name);
                name_to_seq.erase(n_to_s);
                graph_modified = true;
            }
            else if (min_len > 0 && s.length < min_len){
                dropped_seqs.insert(s.name);
                name_to_seq.erase(n_to_s);
                graph_modified = true;
            }

            if (no_amb && !pliib::canonical(s.sequence)){
                dropped_seqs.insert(s.name);
                name_to_seq.erase(n_to_s);
                graph_modified = true;
            }
        }

        map<string, vector<edge_elem>>::iterator s_to_e;
        for (s_to_e = seq_to_edges.begin(); s_to_e != seq_to_edges.end(); s_to_e++){
            // Remove all edges that have dropped sequences as their source.
            if (dropped_seqs.find(s_to_e->first) != dropped_seqs.end()){
                seq_to_edges.erase(s_to_e);
                graph_modified = true;
                continue;
            }
            // Remove any individual edges that have a dropped sequence as their sink.
            vector<edge_elem> pass_elems;
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



    void GFAKluge::compatibilize(){
        gfa_1_ize();
        gfa_2_ize();
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


    string GFAKluge::join(const vector<string>& splits, const string& glue){
        string ret = "";
        for (size_t i = 0; i < splits.size(); i++){
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


    string GFAKluge::opt_string(vector<opt_elem> opts){
        string ret = "";
        for (size_t i = 0; i < opts.size(); i++){
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

			//First print header lines.
			if (header.size() > 0){
					ret << header_string(header) + "\n";
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
            
            if (name_to_path.size() > 0 && this->version == 1.0){
                map<string, path_elem>::iterator pt;
                for (pt = name_to_path.begin(); pt != name_to_path.end(); ++pt){
                    stringstream pat;
                    pat << "P\t" << pt->second.name << "\t";
                    vector<string> ovec;
                    for (size_t oi = 0; oi < pt->second.segment_names.size(); oi++){
                        stringstream o_str;
                        o_str << pt->second.segment_names[oi] << (pt->second.orientations[oi] ? "+" : "-");
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
            else if (this->version < 1.0){
                for (auto p : name_to_path){
                    stringstream st;
                    p.second.write_as_walks(st);
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
        //First print header lines.
        if (header.size() > 0){
            ret << header_string(header) + "\n";
        }
        
        if (name_to_path.size() > 0 && this->version >= 1.0){
                map<string, path_elem>::iterator pt;
                for (pt = name_to_path.begin(); pt != name_to_path.end(); ++pt){
                    stringstream pat;
                    pat << "P" << "\t";
                    pat << pt->second.name << "\t";
                    vector<string> ovec;
                    for (size_t seg_ind = 0; seg_ind < pt->second.segment_names.size(); seg_ind++){
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


            // Paths
            
        
            
            uint64_t edge_count = 0;
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
        for (size_t i = 0; i < s_lens.size(); i++){
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
        for (size_t i = 0; i < s_lens.size(); i++){
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

    int GFAKluge::get_L90(){
        vector<double> s_lens;
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

            if (name_to_path.size() > 0 && this->version >= 1.0){
                map<string, path_elem>::iterator pt;
                for (pt = name_to_path.begin(); pt != name_to_path.end(); ++pt){
                    stringstream pat;
                    pat << "P\t" << pt->second.name << "\t";
                    vector<string> ovec;
                    for (size_t oi = 0; oi < pt->second.segment_names.size(); oi++){
                        stringstream o_str;
                        o_str << pt->second.segment_names[oi] << (pt->second.orientations[oi] ? "+" : "-");
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
            else if (this->version < 1.0 && name_to_path.size() > 0){
                map<string, path_elem>::iterator pt;
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
                    map<string, path_elem>::iterator pt;
                    for (pt = name_to_path.begin(); pt != name_to_path.end(); ++pt){
                        stringstream pat;
                        pat << "P" << "\t";
                        pat << pt->second.name << "\t";
                        vector<string> ovec;
                        for (size_t seg_ind = 0; seg_ind < pt->second.segment_names.size(); seg_ind++){
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
            else if (this->version < 1.0 && name_to_path.size() > 0){
                map<string, path_elem>::iterator pt;
                for (pt = name_to_path.begin(); pt != name_to_path.end(); pt++){
                    pt->second.write_as_walks(os);
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

    std::string GFAKluge::header_string(){
        std::stringstream st;
        for (auto h : get_header()){
            write_element(st, h.second);
        }
        return st.str();
    }

    void GFAKluge::write_element(std::ostream& os, const header_elem& h){
        os << h.to_string();
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
