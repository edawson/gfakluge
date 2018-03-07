gfakluge
--------------------

## What is it?  
GFAKluge is a C++ parser/writer for [GFA files](http://lh3.github.io/2014/07/19/a-proposal-of-the-grapical-fragment-assembly-format/). It parses
GFA to a set of data structures that represent the encoded graph.
You can use these components and their fields/members to build up your own
graph representation.  

A set of command line utilities for basic manipulations is also included.

## Command line utilities
When `make` is run, a bin directory is added with the following utilities:
[X] gfa extract : transform the GFA segment lines to a FASTA file.  
[X] gfa diff : check if two GFA files are different (not very sophisticated at the moment)  
[X] gfa sort : change the line order of a GFA file so that lines proceed in
Header -> Segment -> Link/Edge/Containment -> Path order.  
[X] gfa spec\_convert : convert between the different GFA specifications (e.g. GFA1 -> GFA2).  
[X] gfa stats : get the assembly stats of a GFA file (e.g. N50, L50)  
[X] gfa subset : extract a subgraph between two Segment IDs in a GFA file.  
[ ] gfa ids : manually coordinate / increment the ID spaces of two graphs, so that they can be concatenated.  
[ ] gfa merge : merge (i.e. concatenate) multiple GFA files.  

For CLI usage, run any of the above with no arguments or `-h`.

## Why gfaKluge?
Currently, most GFA parsers go directly from file to proprietary internal graph representations.
This library instead parses to standard C++ STL containers.
This means it is not beholden to any internal representation while
still providing easy I/O of GFA.
Since it only relies on the STL, it's easy to build and requires only what's included in the repo.

## How do I build it?  
You can build libgfakluge by typing ``make`` in the repo. 
To use GFAKluge in your program, you'll need to
add a few lines to your code. First, add the necessary include line to your C++ code:  
                #include "gfa_kluge.hpp"

Next, make sure that the library is on the proper system paths and compile line:

                g++ -o my_exe my_exe.cpp -L/path/to/gfakluge/ -lgfakluge


You should then be able to parse and manipulate gfa from your program:  

                    gg = GFAKluge();
                    gg.parse_gfa_file(my_gfa_file); 

                    cout << gg << endl;


## Internal Structures
Internally, lines of GFA are represented as structs with member variables that correspond to their defined fields.
Here's the definition for a sequence line, for example:

                struct sequence_elem{
                    std::string seq;
                    std::string name;
                    map<string, string> opt_fields;
                    long id;
                };

The structs for contained elements, link elements, and alignment elements are very similar. These individual structs
are then wrapped in a set of standard containers for easy access:

                map<std::string, std::string> header;
                map<string, sequence_elem> name_to_seq;
                map<std::string, vector<contained_elem> > seq_to_contained;
                map<std::string, vector<link_elem> > seq_to_link;
                map<string, vector<alignment_elem> > seq_to_alignment;

All of these structures can be accessed using the ``get_<Thing>`` method, where \<Thing\> is the name of the map you would like to retrieve.  

## GFA2
GFAKluge now supports GFA2! This brings with it four new structs: `edge_elem`, `gap_elem`, `fragment_elem`, and `group_elem`. They're contained in maps much like those for the GFA1 structs.  

A few caveats apply:  
    1. As GFA2 is a **superset** of GFA1, we support only support legal GFA2 -> GFA1 conversions. Information can be lost along the way (e.g. unordered groups won't be output).
    2. Our GFA2 testing is a bit limited but we've verified several times to be on-spec.

Tags we specifically do not (i.e. cannot) support in GFA2 -> GFA1 conversion: G - gap, U - unordered group, F - fragment.
Links and containments should get converted to edges correctly. Sequence elements should get converted, but watch out for the length field if you hit issues.

GFAKluge is fully compliant with reading GFA2 and GFA0.1 <-> GFA1.0 -> GFA2.0 conversion as of September 2017.

## Reading GFA
                GFAKluge gg;
                gg.parse_gfa_file("my_gfa.gfa");

You can then iterate over the aforementioned maps/structs and build out your own graph representation.

I'm working on a low-memory API for reading lines / emitting structs but it won't be this pretty.

## Writing GFA
                GFAKluge og;

                sequence_elem s;
                s.sequence = "GATTACA";
                s.name = "seq1";
                og.add_sequence(s);

                sequence_elem t;
                t.sequence = "AATTGN";
                t.name = "seq2";
                og.add_sequence(t);

                link_elem l;
                l.source = s.name;
                l.sink = s.name;
                l.source_orientation_forward = true;
                l.sink_orientation_forward = true;
                l.pos = 0;
                l.cigar = "";

                og.add_link(l.source, l);

                cout << og << endl;
                ofstream f = ofstream("my_file.gfa);
                f << og;

## Status
- GFAKluge is essentially a set of dumb containers - it does no error checking of your structure to detect if it is
valid GFA. This may change as the GFA spec becomes more formal.  
- Parses JSON structs in optional fields of sequence lines (just as strings though).  
- Full GFA1/GFA2 compatibility and interconversion is now implemented.  


## Getting Help 
Eric T Dawson  
github: [edawson](https://github.com/edawson/https://github.com/edawson/GFAKluge)  
Please post an issue for help.
