#include "gfakluge.hpp"

using namespace std;
using namespace gfak;

GFAKluge::GFAKluge(){

}

GFAKluge::~GFAKluge(){

}

// Borrow from
//http://stackoverflow.com/questions/236129/split-a-string-in-c
// Thanks StackOverflow!
vector<string> split(string s, char delim){
    vector<string> ret;
    stringstream sstream;
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
        return false;
    }
    string line;
    vector<string> line_tokens;
    while (getline(gfi, line)){
        vector<string> tokens = split(line, '\t');
        //switch (tokens[0]){
        if (tokens[0] == "H"){
            header_elem h;
            line_tokens = split(tokens[1], ':');
            h.key = line_tokens[0];
            h.val = line_tokens[1];
            header[h.key] = h.val;
        }
        else if (tokens[0] ==  "S"){
            //TODO: we've got some tokens at the end of the line
            //that have not been handled yet.
            sequence_elem s;
            s.name = tokens[1];
            s.seq = tokens[2];
            s.id = atol(s.name.c_str());
            name_to_seq[s.name] = s;
        }
        else if (tokens[0] ==  "L"){
            // TODO: we need to deal with links where the link is given before
            // its corresponding sequence in the file. TODO this is probably
            // now fixed by using the string: sequence map.
            link_elem l;
            l.source_name = tokens[0];
            l.sink_name = tokens[0];
            //TODO: search the input strings for "-" and "+" and set using ternary operator
            l.source_orientation = tokens[0];
            l.sink_orientation = tokens[0];
            //l.pos = tokens[0];
            l.cigar = tokens[0];
            add_link(l.source_name, l);
        }
        else if (tokens[0] == "C"){
            contained_elem c;
            //TODO fix token indices here
            c.source_name = tokens[0];
            c.sink_name = tokens[0];
            c.source_orientation = tokens[0];
            c.sink_orientation = tokens[0];
            c.pos = atoi(tokens[0].c_str());
            c.cigar = tokens[0];
            add_contained(c.sink_name, c);
        }
        else if (tokens[0] == "x"){
            annotation_elem x;
            x.key = tokens[0];
            x.info = tokens[1];
        }
        else if (tokens[0] == "a"){
            annotation_elem a;
            a.key = tokens[0];
            a.info = tokens[0];
        }
        else{
        }

    }

    return true;

    }

    bool GFAKluge::parse_gfa_file(fstream fs){
        cerr << "Not implemented: parse_gfa_file(fstream)" << endl; exit(1);
        return true;
    }

    void GFAKluge::add_link(sequence_elem seq, link_elem link){
        seq_to_link[seq.name].push_back(link);
    }

    void GFAKluge::add_contained(sequence_elem seq, contained_elem con){
        seq_to_contained[seq.name].push_back(con);
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

    string GFAKluge::to_string(){
        string ret = "";
        //First print header lines.
        map<std::string, std::string>::iterator it;
        for (it = header.begin(); it != header.end(); it++){
            ret += "H " + it->first + " " + it->second + "\n";
        }

        //Print sequences and links in order, then annotation lines.
    }
