gfakluge
--------------------

[![DOI](https://zenodo.org/badge/47132080.svg)](https://zenodo.org/badge/latestdoi/47132080)
[![status](http://joss.theoj.org/papers/d731f6dfc6b77013caaccfd8333c684a/status.svg)](http://joss.theoj.org/papers/d731f6dfc6b77013caaccfd8333c684a)

## What is it?  
GFAKluge is a C++ parser/writer and
a set of command line utilities for manipulating [GFA files](http://lh3.github.io/2014/07/19/a-proposal-of-the-grapical-fragment-assembly-format).
It parses GFA to a set of data structures that represent the encoded graph.
You can use these components and their fields/members to build up your own
graph representation. You can also convert between GFA 0.1 <-> 1.0 <-> 2.0
to glue programs that use different GFA versions together.


**Homepage**: https://github.com/edawson/gfakluge  
**License**: MIT  

## Dependencies
A C++11 compliant compiler (we recommend GCC or clang)  
OpenMP (via GCC or clang)  
**NB**: GFAKluge cannot be compiled with Apple clang, as it does not include OpenMP.

## Command line utilities
When `make` is run, the `gfak` binary is built in the top level directory. It offers the following subcommands:  
+ gfak extract : transform the GFA segment lines to a FASTA file.  
+ gfak fillseq : fill in the sequence field of S lines with placeholders using sequences from a FASTA file.
+ gfak diff : check if two GFA files are different (not very sophisticated at the moment)  
+ gfak sort : change the line order of a GFA file so that lines proceed in
Header -> Segment -> Link/Edge/Containment -> Path order.  
+ gfak convert : convert between the different GFA specifications (e.g. GFA1 -> GFA2).  
+ gfak stats : get the assembly stats of a GFA file (e.g. N50, L50)  
+ gfak subset : extract a subgraph between two Segment IDs in a GFA file.  
+ gfak ids : manually coordinate / increment the ID spaces of two graphs, so that they can be concatenated.  
+ gfak merge : merge (i.e. concatenate) multiple GFA files. NB: Obliterates nodes with the same ID.  

For CLI usage, run any of the above (including `gfak` with no subcommand) with no arguments or `-h`. To change specification version, most commands take the `-S` flag and a single `double` argument.  

## Example CLI Usage
Examples of various commands are included in the [examples.md file](https://github.com/edawson/gfakluge/blob/master/examples.md).

## C++ API
Examples of the C++ API are included in the [interface.md file](https://github.com/edawson/gfakluge/blob/master/interface.md).  

## How do I build it?  
The `gfak` utilities are available via homebrew: `brew install brewsci/bio/gfakluge`  

Building GFAKluge from source requires OpenMP. This should be supported on Linux by default. On Apple Mac OS X, we recommend installing gcc:  

```
brew install gcc@8
make CXX=g++-8
```  
or  
```
sudo port install gcc8
make
```

You can then build libgfakluge and the command line `gfak` utilities by typing ``make`` in the repo.  
To use GFAKluge in your program, you'll need to
add a few lines to your code. First, add the necessary include line to your C++ code:  
                #include "gfakluge.hpp"

Next, make sure that the library is on the proper system paths and compile line:

                g++ -o my_exe my_exe.cpp -L/path/to/gfakluge/ -lgfakluge


You should then be able to parse and manipulate gfa from your program:  

                    gg = GFAKluge();
                    gg.parse_gfa_file(my_gfa_file); 

                    cout << gg << endl;


## Why gfak / gfakluge?
+ Simple command line utilities (no awk foo needed!)  
+ High level C++ API for many graph manipulations.  
+ Easy to build - no external dependencies; build with just a modern C++ compiler supporting C++11.
+ Easy to develop with - Backing library is mostly STL containers and a handful of structs.  
+ Performance - gfakluge is fast and relies on standard STL containers and basic structs.  


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
They reside in gfakluge.hpp.  

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
                // Write GFA1
                f << og;

                // To convert to GFA2:
                og.set_version(2.0);
                f << od;

## Status
- GFAKluge is essentially a set of dumb containers - it does no error checking of your structure to detect if it is
valid GFA. This may change as the GFA spec becomes more formal.  
- Diff is not a useful tool yet.
- Parses JSON structs in optional fields of sequence lines (just as strings though).  
- Full GFA1/GFA2 compatibility and interconversion is now implemented.  
- CLI has been refactored to a single executable
- Memory usage for to\_string is a bit high - be careful with large graphs.
- API for input / spec conversion / output is stable. API for merging graphs and coordinating ID namespaces may change slightly, but will strive for backwards compatibility.


## Getting Help 
Eric T Dawson   
github: [edawson](https://github.com/edawson/https://github.com/edawson/GFAKluge)   
Please post an issue for help.

## Contributing
GFAKluge is open-source and community contributions are welcome and appreciated! Please keep the following in mind when contributing to the repo:  

1. Please treat others with kindness and professionalism. Everyone is welcome and we will not tolerate harassment for any reason.
2. Please keep `gfakluge.hpp` header-only and update the build process if a modification alters it.  
3. Please update the dependency list if one is added.  
4. Please use semantic versioning. Minor changes bump the third versioning digit (e.g. 1.0.0 -> 1.0.1).  
Additional features, or changes that may or may not partially break backward compatibility 
but which do not require significant modifications to code depending on the library bump the second versioning digit (e.g. 1.0.0 -> 1.1.0).  
Changes which signficantly alter the API require a bump in the major version digit (e.g. 1.0.0 -> 2.0.0).
5. Please fully specify all namespace items (e.g. `std::stream` in place of just `stream`).  
6. To incorporate changes, please file a pull request on the Github page.
7. Bug reports or feature requests should be posted as "issues" on the Github page with the appropriate tag and referenced in any relevant pull requests.

