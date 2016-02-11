#include "gfakluge.hpp"

using namespace std;
namespace gfak{

    struct map_comp{
        bool operator() (const char& lhs, const char& rhs) const{
            // if ( isdigit(lhs) && isdigit(rhs)){
            //     return (long) lhs > (long) rhs;
            // }
            // else{
                return lhs - rhs;
          //  }
        }
    };


    GFAKluge::GFAKluge(){
        map<std::string, std::string> header;
        map<std::string, vector<contained_elem> > seq_to_contained;
        map<std::string, vector<link_elem> > seq_to_link;
        map<std::string, vector<alignment_elem> > seq_to_alignment;
        //Since we can't compare sequence elements for hashing,
        // we cheat and use their names (which are sort of guaranteed to be
        // unique.
        map<string, sequence_elem, map_comp> name_to_seq;
        map<string, vector<path_elem> > seq_to_paths;
    }

    GFAKluge::~GFAKluge(){

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
							h.val = line_tokens[1];
							header[h.key] = h.val;
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
									for (i = 0; i < tokens.size(); i++){
											//opt fields are in key:val format
											vector<string> key_val = split(tokens[i], ':');
											if (key_val.size() != 2){
													cerr << "WARNING: Unknown pattern in optional field of a sequence entry." << endl;
													cerr << "FIELD WILL BE DISCARDED" << endl;
													continue;
											}
											else{
													s.opt_fields[key_val[0]] = key_val[1];
											}
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
							c.cigar = tokens[6];
							add_contained(c.source_name, c);
					}
          else if (tokens[0] == "P"){
            path_elem p;
            p.name = tokens[2];
            p.source_name = tokens[1];
            //p.rank = ;
            p.is_reverse = tokens[3] == "+" ? false : true;
            p.cigar = tokens[4];
            add_path(p.source_name, p);
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
					else{
							cerr << "Unknown line identifier  encountered: " << tokens[0] <<  " . Exiting." << endl;
					}

			}

			return true;

    }

    void GFAKluge::add_sequence(sequence_elem s){
        name_to_seq[s.name] = s;
    }

    void GFAKluge::add_path(string seq_name, path_elem p){
        seq_to_paths[seq_name].push_back(p);
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

    map<string, sequence_elem> GFAKluge::get_name_to_seq(){
        return name_to_seq;
    }

    map<string, vector<link_elem> > GFAKluge::get_seq_to_link(){
        return seq_to_link;
    }

    map<string, vector<path_elem> > GFAKluge::get_seq_to_paths(){
      return seq_to_paths;
    }

    map<string, vector<contained_elem> > GFAKluge::get_seq_to_contained(){
        return seq_to_contained;
    }

    map<string, vector<alignment_elem> > GFAKluge::get_seq_to_alignment(){
        return seq_to_alignment;
    }

    map<string, string> GFAKluge::get_header(){
        return header;
    }

    string GFAKluge::opt_string(map<string, string>& opts){
        string ret = "";
        map<string, string>::iterator it;
        for (it = opts.begin(); it != opts.end(); it++){
            // TODO needs to add a ';' between fields I think?
            // Using a tab while building without the spec.
            ret += it->first + ":" + it->second + "\t";
        }
        // There must be a better way to remove the last tab. TODO
        // A join method would be nice.
        return ret.substr(0, ret.size() - 1);
    }

    std::string GFAKluge::to_string(){
        string ret = "";
        int i;
        //First print header lines.
        if (header.size() > 0){
            map<std::string, std::string>::iterator it;
            for (it = header.begin(); it != header.end(); it++){
                ret += "H\t" + it->first + "\t" + it->second + "\n";
            }
        }
        if (name_to_seq.size() > 0){
            map<std::string, sequence_elem>::iterator st;
            for (st = name_to_seq.begin(); st != name_to_seq.end(); st++){
                ret += "S\t" + (st->second).name + "\t" + (st->second).sequence + "\n";
                //TODO iterate over links
                //L    segName1,segOri1,segName2,segOri2,CIGAR      Link
                if (seq_to_link[st->first].size() > 0){
                    for (i = 0; i < seq_to_link[st->first].size(); i++){
                        string link = "L\t" + seq_to_link[st->first][i].source_name + "\t";
                        link += seq_to_link[st->first][i].source_orientation_forward ? "+" : "-";
                        link += "\t";
                        link += seq_to_link[st->first][i].sink_name + "\t";
                        link += seq_to_link[st->first][i].sink_orientation_forward ? "+" : "-";
                        link += "\t";
                        link += seq_to_link[st->first][i].cigar + "\n";
                        ret += link;
                    }

                }

                if (seq_to_paths[st->first].size() > 0){
                    for (i = 0; i < seq_to_paths[st->first].size(); i++){
                        string pat = "P\t" + seq_to_paths[st->first][i].source_name + "\t";
                        pat += seq_to_paths[st->first][i].name + "\t";
                        pat += seq_to_paths[st->first][i].is_reverse ? "-" : "+";
                        pat+= "\t";
                        pat += seq_to_paths[st->first][i].cigar + "\n";
                        ret += pat;
                    }
                }

                //TODO iterate over contained segments
                if (seq_to_contained[st->first].size() > 0){
                    for (i = 0; i < seq_to_contained[st->first].size(); i++){
                        string cont = "C\t" + seq_to_contained[st->first][i].source_name + "\t";
                        cont += seq_to_contained[st->first][i].source_orientation_forward ? "+" : "-";
                        cont += "\t";
                        cont += seq_to_contained[st->first][i].sink_name + "\t";
                        cont += seq_to_contained[st->first][i].sink_orientation_forward ? "+" : "-";
                        cont += "\t";
                        cont += seq_to_contained[st->first][i].cigar + "\n";
                        ret += cont;
                    }
                }
            }


        }
        //TODO iterate over annotation lines.


        //Print sequences and links in order, then annotation lines.

        return ret;
    }

     std::ostream& operator<<(std::ostream& os, GFAKluge& g){
         os << g.to_string();
         return os;
     }
}
