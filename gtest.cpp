#include "gfakluge.hpp"

int main() {
    gfak::GFAKluge og;

    og.set_version(1);

    gfak::sequence_elem s;
    s.sequence = "ACCTT";
    s.name = "11";

    gfak::sequence_elem t;
    t.sequence = "TCAAGG";
    t.name = "12";

    gfak::sequence_elem u;
    u.sequence = "CTTGATT";
    u.name = "13";

    gfak::link_elem l;
    l.source_name = s.name;
    l.sink_name = t.name;
    l.source_orientation_forward = true;
    l.sink_orientation_forward = false;
    l.cigar = "4M";

    gfak::link_elem m;
    m.source_name = t.name;
    m.sink_name = u.name;
    m.source_orientation_forward = false;
    m.sink_orientation_forward = true;
    m.cigar = "5M";

    gfak::link_elem n;
    n.source_name = s.name;
    n.sink_name = u.name;
    n.source_orientation_forward = true;
    n.sink_orientation_forward = true;
    n.cigar = "3M";

    gfak::path_elem p;
    p.name = "14";
    p.segment_names = {s.name, t.name, u.name};
    p.orientations = {true, false, true};
    p.overlaps = {"4M", "5M"};

    og.add_sequence(s);
    og.add_sequence(t);
    og.add_sequence(u);
    og.add_link(s, l);
    og.add_link(t, m);
    og.add_link(s, n);
    og.add_path(p.name, p);

    og.gfa_2_ize();
    og.gfa_1_ize();
    og.set_version(2.0);
    cout << og << endl;
  
}
