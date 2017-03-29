#include "gfakluge.hpp"

using namespace std;
namespace gfak{


    GFAKluge::GFAKluge(){
        map<std::string, header_elem> header;
        map<std::string, vector<contained_elem> > seq_to_contained;
        map<std::string, vector<link_elem> > seq_to_link;
        map<std::string, vector<alignment_elem> > seq_to_alignment;
        //Since we can't compare sequence elements for hashing,
        // we cheat and use their names (which are sort of guaranteed to be
        // unique.
        map<string, sequence_elem, custom_key> name_to_seq;
        map<string, vector<path_elem> > seq_to_walks;
        map<string, path_elem> name_to_path;
        //cout << fixed << setprecision(2);

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
        regex forward_re("\\+");
        regex reverse_re("\\-");
        regex id_re("[0-9A-Za-z]*");
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
                sequence_elem s;
                s.name = tokens[1];
                s.sequence = tokens[2];
                //s.id = atol(s.name.c_str());
                int i;
                if (tokens.size() > 3){
                    for (i = 3; i < tokens.size(); i++){
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
                l.cigar = tokens[5];
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
                        bool orientation = regex_search(x, reverse_re);
                        string id = x.substr(0, x.length() - 1);
                        p.segment_names.push_back(id);
                        p.orientations.push_back(orientation);
                    }

                    if (tokens.size() > 3){
                        p.overlaps = split(tokens[3], ',');
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
                        w.rank = atol(tokens[3].c_str());
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

        paths_as_walks();
        walks_as_paths();

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

    map<string, vector<walk_elem> > GFAKluge::paths_as_walks(){
        map<string, path_elem>::iterator it;
        for (it = name_to_path.begin(); it != name_to_path.end(); it++){
            // Create a Walk for every segment in the path's segment list
            for (int i = 0; i < it->second.segment_names.size(); i++){
                walk_elem w;
                w.name = it->second.name;
                w.source_name = it->second.segment_names[i];
                if (!it->second.overlaps.empty()){
                    w.cigar = it->second.overlaps[i];
                    w.is_reverse = (it->second.orientations[i]);
                }
                add_walk(w.source_name, w);
            }
        }
        return seq_to_walks;
    }

    map<string, path_elem> GFAKluge::walks_as_paths(){
        map<string, vector<walk_elem> >::iterator it;
        map<string, vector<string> > pathname_to_seqnames;
        map<string, vector<string> > pathname_to_overlaps;
        map<string, vector<bool> > pathname_to_reverse;
        for (it = seq_to_walks.begin(); it != seq_to_walks.end(); it++){
            for (int i = 0; i < it->second.size(); i++){
                walk_elem w = it->second[i];
                pathname_to_seqnames[w.name].push_back(w.source_name);
                pathname_to_overlaps[w.name].push_back(w.cigar);
                pathname_to_reverse[w.name].push_back(w.is_reverse);
            }
        }
        for (auto iter : pathname_to_seqnames){
            string name = iter.first;
            path_elem p;
            p.name = name;
            p.segment_names = pathname_to_seqnames[name];
            p.orientations = pathname_to_reverse[name];
            for (int i = 0; i < pathname_to_overlaps.size(); i++){
                stringstream ssr;
                string cigar = pathname_to_overlaps[name][i];
                ssr << cigar;
                // if (pathname_to_reverse[name][i]){
                //     ssr << "R";
                // }
                pathname_to_overlaps[name][i] = ssr.str();
            }
            p.overlaps = pathname_to_overlaps[name];
            add_path(p.name, p);

        }
        return name_to_path;

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
			stringstream ret;
			int i;
			//First print header lines.
			if (header.size() > 0){
					ret << header_string(header) + "\n";
			}
            //walks_as_paths();
            //paths_as_walks();
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
                    pat << "P\t" << pt->second.name << "\t" << join(pt->second.segment_names, ",");
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


    //TODO this should use stringstream too...
    std::string GFAKluge::to_string(){
        stringstream ret;
        int i;
        //First print header lines.
        if (header.size() > 0){
            ret << header_string(header) + "\n";
        }

        //paths_as_walks();
        //walks_as_paths();
        if (name_to_path.size() > 0 && this->version >= 1.0){
                map<string, path_elem>::iterator pt;
                for (pt = name_to_path.begin(); pt != name_to_path.end(); ++pt){
                    stringstream pat;
                    pat << "P\t" << pt->second.name << "\t" << join(pt->second.segment_names, ",");
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

    std::ostream& operator<<(std::ostream& os, GFAKluge& g){
        os << g.to_string();
        return os;
    }
}
