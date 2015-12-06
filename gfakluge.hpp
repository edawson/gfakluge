#ifndef GFAK_HPP
#define GFAK_HPP
#include <string>
#include <fstream>


using namespace std;
struct header_elem{
    string key;
    string val;
};

struct sequence_elem{
    // A list of links to neighbors

    // Node sequence
    string seq;
    string name;
    long id;

};

struct link_elem{

    // Source sequence
    string ssource_name;
    string sink_name;
    bool source_orientation;
    bool sink_orientation;
    string cigar;
    // Sink sequence
    //
    //
};

struct contained_elem{
    string ssource_name;
    string sink_name;
    bool source_orientation;
    bool sink_orientation;
    int pos;
    string cigar;

};

// TODO devise a clever data structure for gfa entries
class GFAK{
    public:
        bool set_gfa_file(string filename);
        bool parse_gfa_file(fstream gfa_stream);
        map<sequence_elem, vector<link_elem>> seq_to_links;
        map<sequence_elem, vector<contained_elem>> seq_to_contained;
    private:


};

#endif
