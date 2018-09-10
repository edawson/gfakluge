The GFAKluge developer API
----------------------------

## Introduction
The GFAKluge.hpp header file defines an object-orientated API for managing graphs in GFA format. It also defines structs for each of the GFA line types.

## GFAKluge class and functions
The GFAKluge class defines a collection of structs representing a GFA file.

## Adding elements to a GFAKluge object
The GFAKluge object has a series of "add\_\*" methods which add a given element to the graph object.

- add\_sequence(sequence\_elem s) adds a sequence element to the graph
- add\_edge(string s, edge\_elem e) and add\_edge(sequence\_elem s, edge\_elem e) add an
  edge\_element with source *s* to the graph.
- add\_link(sequence\_elem s, link\_elem l) and add\_link(string s, link\_elem l) add an
  edge\_element that represents a link starting at source *s*.
- add\_contained(sequence\_elem s, contained\_elem c) and add\_contained(string s, contained\_elem c) add a contained\_elem that falls within *s*.
- add\_fragment(sequence\_elem s, fragment\_elem f) and add\_fragment(string s, fragment\_elem f) add
  a fragment\_element that overlaps *s*.
- add\_gap(gap\_elem g) adds a gap to the GFAKluge object by its source.
- add\_group(group\_elem g) adds a group\_elem to the GFAKluge object by its name.
- add\_path(string s, path\_elem p) adds a path to the GFAKluge object.
- add\_walk(std::string pathname, const int& rank, const string& segname, const bool& ori, const string& overlap, vector\<opt\_elem\> opts) adds a walk element (which is a single element on a path) to  a path in the GFAKluge object with name *pathname*.

## Retrieving elements from a GFAKluge object
The raw maps containing elements in a GFAKluge object can be retrieved using the following methods:

- get\_name\_to\_seq() : returns a map from the string name of a sequence element to the element itself.
- get\_groups() : returns a map from the string name of a group to the group\_elem itself.
- get\_name\_to\_path() : returns a map from the string name of a path to the corresponding path\_elem.

The other get\_\* methods return maps where the key is the string name of a source sequence element.  
- get\_seq\_to\_edges() : returns a map\<string, vector\<edge_elem\>\>, where all edge\_elems are stored relative to their source sequence element.
- get\_seq\_to\_fragments() : returns a map from a string identifier for a sequence element to fragments that fall within that element.
- get\_seq\_to\_gaps() : returns a map from a string identifier for a sequence element to gaps that fall within that element.
- get\_alignments(sequence\_elem s) : returns the alignments (GFA0.1) to a sequence element s.


We would encourage developers to look at the `std::string to_string_2();` and `std::string block_order_string();` methods in the `gfakluge.cpp` file for an example of how to iterate over elements in these maps.


## GFA2 vs. GFA1
There are notable and important differences between GFA1 and GFA2:
  - GFA1 L and C lines can be represented as GFA2 E lines. However, E lines may not always
  be valid C or L lines.
  - GFA2 has unordered groups (U lines), which GFA1 only has ordered groups (paths / P lines).
  Unordered groups are lost on conversion from GFA2 -> GFA2.
  - GFA2 has G lines (gaps) and F lines (fragments). These are lost on conversion from
  GFA2 -> GFA1.
  
GFAKluge represents all lines (except P lines) internally as GFA2-compatible structs. Paths
are stored as path\_elems but can be converted to ordered groups.

In addition, pre-GFA1 versions may use 'a' and 'x' lines, which are outside the official spec.
These lines cannot be expressed in GFA1 or GFA2 and are lost upon conversion.


## Data Structures
GFAKluge is a collection of maps that store graph structure information. The atomic data structures
of GFAKluge are structs with member variables representing the different fields in a GFA line. These
are defined in gfakluge.hpp with the suffix "\_elem."
The details of the backing "\*\_elem" structures are described below.

### sequence\_elem
The sequence\_elem struct defines a node in a GFA graph. It has the following members:
  1. *name* : a string identifier for the sequence.
  2. *sequence* : the base pairs the node contains. "\*" may also be used as a placeholder when sequences are held in an external FASTA.
  3. *length* : A uint64_t containing the length of *sequence*.
  4. *opt_fields* : a vector of opt\_elems (key:type:value triples) which describe annotations or tags on this sequence element.
  
  sequence\_elem structs have the following methods:
  1. to\_string\_1 : return a GFA1 string of the sequence elem.
  2. to\_string\_2 : return a GFA2 string of the sequence elem.
  3. as\_fasta\_record : return a valid FASTA record (name and sequence)for the sequence element.

### edge\_elem

### group\_elem

### fragment\_elem

### gap\_elem

### path\_elem

## gfa\_builder.hpp
The gfa\_builder.hpp header file defines a set of functions used to build variation graphs
from a FASTA file and a VCF file. These functions support the `gfak build` command, but we
discuss the use of this interface outside of GFAKluge.
