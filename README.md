gfakluge
--------------------

## What is it?  
GFAKluge is a C++ parser/writer for [GFA files](http://lh3.github.io/2014/07/19/a-proposal-of-the-grapical-fragment-assembly-format/). It parses
GFA to a set of data structures that represent the encoded graph.
You can use these components and their fields/members to build up your own
graph representation.



## How do I build it?  
You can build libgfakluge by simply typing ``make``. 
To use GFAKluge in your program, you'll need to
add a few lines to your code. First, add the necessary include line to your C++ code:  
				#include "gfa_kluge.hpp"

Nest, make sure that the library is on the proper system paths and compile line:
                g++ -o my_exe my_exe.cpp -L/path/to/gfakluge/ -lgfakluge
		

You should then be able to parse and manipulate gfa from your program:
                gg = GFAKluge();
                gg.parse_gfa_file(my_gfa_file); 

                cout << gg << endl;

## Status
Currently, parsing is implemented but the data structures are not publicly exposed. As of January 2016
I'm in the process of exposing them so GFA can be parsed into structures such as those used by [vg](https://github.com/ekg/vg.git).

## Who do I bother if it's broken?  
Eric T Dawson  
github: [edawson](https://github.com/edawson/https://github.com/edawson/GFAKluge)  
Please post an issue for help or file a bug report to report a bug.
