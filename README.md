GFAKluge
--------------------

## What is it?  
GFAKluge is a parser/manipulator for [GFA files](http://lh3.github.io/2014/07/19/a-proposal-of-the-grapical-fragment-assembly-format/). It parses
GFA to a set of graph components that represent the encoded graph.
You can use these components and their fields/members to build up your own
graph representation.


You can also output GFA by filling in these simple structures and then
converting them to GFA with GFAKluge.


## How do I build it?  
You can build a toy example by simply cloning the repo, entering it
and typing ```make; make test```. To use GFAKluge in your program, you'll need to
add a few lines to your code. First, add the necessary include line to your C++ code.
				#include "gfa_kluge.hpp"

Then, add a make target to you makefile and make sure that the object/include files
are on the proper system paths:
				gfa_kluge.o: gfakluge.cpp gfakluge.hpp
							$(CC) $(CFLAGS) -c -o gfakluge.o gfakluge.cpp

You should then be able to parse and manipulate gfa from your program.

## Who do I bother if it's broken?  
Eric T Dawson  
github: [edawson](https://github.com/edawson/https://github.com/edawson/GFAKluge)
Please post an issue for help or file a bug report to report a bug.
