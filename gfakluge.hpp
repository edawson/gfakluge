#ifndef GFAK_HPP
#define GFAK_HPP
#include <string>
#include <fstream>
#include <map>
#include <vector>

using namespace std;
namespace gfak{
    struct header_elem{
        std::string key;
        std::string val;
    };

    struct annotation_elem{
        std::string key;
        std::string info;
    };

    struct sequence_elem{
        // A list of links to neighbors

        // Node sequence
        std::string seq;
        std::string name;
        long id;

    };

    struct link_elem{

        // Source sequence
        std::string source_name;
        // Sink sequence
        std::string sink_name;
        bool source_orientation;
        bool sink_orientation;
        std::string cigar;
    };

    struct contained_elem{
        std::string source_name;
        std::string sink_name;
        bool source_orientation;
        bool sink_orientation;
        int pos;
        std::string cigar;
    };

    class GFAKluge{
        public:
            bool parse_gfa_file(std::string filename);
            bool parse_gfa_file(std::fstream gfa_stream);
            void add_link(sequence_elem, link_elem);
            void add_contained(sequence_elem, contained_elem);
            vector<link_elem> get_links(sequence_elem seq);
            vector<contained_elem> get_contained(sequence_elem seq);
            std::string to_string();
        private:
            map<std::string, std::string> header;
            map<sequence_elem, vector<contained_elem> > seq_to_contained;
            map<sequence_elem, vector<link_elem> > seq_to_link;



    };
    std::ostream& operator<<(std::ostream& os, const GFAKluge g);
};
#endif
