The GFAKluge developer API
----------------------------


## Introduction
The GFAKluge.hpp header file defines an object-orientated API for managing graphs in GFA format. It also defines structs for each of the GFA line types.

## GFAKluge class and functions
The GFAKluge class defines a collection of structs representing a GFA file.

## Adding elements to a GFAKluge object
The GFAKluge object has a series of "add\_\*" methods which add a given element to the graph object.

- `add_sequence(sequence_elem s)` adds a sequence element to the graph
- `add_edge(string s, edge_elem e)` and add\_edge(sequence\_elem s, edge\_elem e) add an
  edge\_element with source *s* to the graph.
- `add_link(sequence_elem s, link_elem l)` and add\_link(string s, link\_elem l) add an
  edge\_element that represents a link starting at source *s*.
- `add_contained(sequence_elem s, contained_elem c)` and add\_contained(string s, contained\_elem c) add a contained\_elem that falls within *s*.
- `add_fragment(sequence_elem s, fragment_elem f)` and add\_fragment(string s, fragment\_elem f) add
  a fragment\_element that overlaps *s*.
- `add_gap(gap_elem g)` adds a gap to the GFAKluge object by its source.
- `add_group(group_elem g)` adds a group\_elem to the GFAKluge object by its name.
- `add_path(string s, path_elem p)` adds a path to the GFAKluge object.
- `add_walk(std::string pathname, const int& rank, const string& segname, const bool& ori, const string& overlap, vector<opt_elem> opts)` adds a walk element (which is a single element on a path) at position *rank* on a path in the GFAKluge object with name *pathname*.

## Retrieving elements from a GFAKluge object
The raw maps containing elements in a GFAKluge object can be retrieved using the following methods:

- `get\_name\_to\_seq()` : returns a map from the string name of a sequence element to the element itself.
- `get\_groups() : returns a map from the string name of a group to the group\_elem itself.
- `get\_name\_to\_path()` : returns a map from the string name of a path to the corresponding path\_elem.

The other get\_\* methods return maps where the key is the string name of a source sequence element.  
- `get\_seq\_to\_edges()` : returns a map\<string, vector\<edge_elem\>\>, where all edge\_elems are stored relative to their source sequence element.
- `get\_seq\_to\_fragments()` : returns a map from a string identifier for a sequence element to fragments that fall within that element.
- `get\_seq\_to\_gaps()` : returns a map from a string identifier for a sequence element to gaps that fall within that element.
- `get\_alignments(sequence\_elem s)` : returns the alignments (GFA0.1) to a sequence element s.


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
  2. *sequence* : the base pairs the node contains (in a C++ string). "\*" may also be used as a placeholder when sequences are held in an external FASTA.
  3. *length* : A uint64_t containing the length of *sequence*.
  4. *opt_fields* : a vector of opt\_elems (key:type:value triples) which describe annotations or tags on this sequence element.
  
sequence\_elem structs have the following methods:
  1. *to\_string\_1* : return a GFA1 string of the sequence elem.
  2. *to\_string\_2* : return a GFA2 string of the sequence elem.
  3. as\_fasta\_record : return a valid FASTA record (name and sequence)for the sequence element.

### edge\_elem
The edge\_elem struct defines an edge (link or containment) in a GFA graph. It has the following members:
  1. *type* : an integer that describes what type of edge this elem represents. A link is given type = 1; a containment, type = 2; unknown, type = 3; unset, type = 0.
  2. *source\_name* : The string name of the source sequence\_elem (node) from which this edge originates.
  3. *sink\_name* : The string name of the sink sequence\_elem (node) to which this edge goes.
  4. *source\_orientation\_forward* : A boolean which describes if the edge starts at the beginning (false) or end (true) of its source node.
  5. *sink\_orientation\_forward* : A boolean which describes if the edge starts at the beginning (true) or end (false) of its sink node.
  6. *ends* :A bitset that holds four (boolean) values indicating whether an edge falls at the end(s) of its source/sink.
  7. *tags* : a vector opt\_elem structs that holds annotations on the edge.
  8. *alignment* : a string field that describes an overlaps.
  9. *source_begin* : A 64-bit int that describes the position in the edge's source node where it begins.
  10. *source_end* : A 64-bit int that describes the position in the edge's source node where it ends.
  11. *sink_begin* : A 64-bit int that describes the position in the edge's sink node where it begins.
  12. *sink_end* : A 64-bit int that describes the position in the edge's sink node where it ends.

edge\_elem structs have the following methods:
  1. *to\_string\_1()* : Return a string representing the edge as a GFA1 'L' or 'C' line, depending on its *type*. **NB** Edges where *type* is unset or unknown are lost upon conversion to GFA1.
  2. *to\_string\_2()* : Returns a string GFA2 'E' line representing the edge.
  3. *determine\_type()* : Attempts to set the *type* of the edge based on its *ends*. If the *type* is unclear, it is set to 3 (unknown).

### group\_elem
The group\_elem struct defines an ordered group ("path") or an unordered group ("collection") of sequence elements. It has the following members:
  1. *id* : A string identifier for the group.
  2. *ordered* : A boolean indicating whether this group is ordered (the group is a set, or path) or unordered (the group is a collection).
  3. *items* : A vector containing the string IDs of the sequence\_elems in the group.
  4. *orientations* : A vector containing boolean values, one for each item in *items*, that describes whether an item is in the forward or reverse direction. This vector may be empty if the group is not ordered.
  5. *tags* : A map from string to opt\_elem structs which holds annotations on the group.

group\_elem structs have the following methods:
  1. *to\_string\_1()* : Returns a string representation of the group as a GFA1 'P' (path) line. Unordered groups return an empty string.
  2. *to\_string\_2()* : Returns a string representation of the group as a GFA2 'U' or 'O' line.
  3. *to\_walk\_string()* : Returns a string representation of the group as GFA0.1 / GFA1 walks, which is a multi-line description of a path.

### fragment\_elem
fragment\_elem structs define fragments, which are containments or alignments of a sequence within another sequence in GFA2. They have the following members:
  1. *id* : A string ID for the fragment
  2. *ref* : A string name for the sequence\_elem which contains the fragment.
  3. *ref\_orientation* : A boolean describing whether the fragment is in the forward (true) or reverse (false) direction.
  4. *seg\_begin* : An unsigned 32-bit integer representing the base pair within the ref where the fragment's alignment begins.
  5. *seg\_end* : An unsigned 32-bit integer representing the base pair within the ref where the fragment's alignment ends.
  6. *frag\_begin* : An unsigned 32-bit integer representing the base pair within the fragment where its alignment to *ref* begins.
  7. *frag\_end* : An unsigned 32-bit integer representing the base pair within the fragment where its alignment to *ref* ends.
  8. *ends* : A bitset which marks whether the fragment's starts/ends align with the ends of the segment.
  9. *alignment* :
  10. *tags* :

fragment\_elem structs have the following methods:
  1. *to\_string\_2()* : returns a string of the GFA2 'F' line the struct represents.
  2. *to\_string()* : a wrapper for *to\_string\_2()*.


### gap\_elem
gap\_elem structs define a gap in a GFA2 graph. They have the following members:
  1. *id* : A string ID for the gap.
  2. *source\_name* : A string ID for the sequence\_elem at the gap's source.
  3. *sink\_name* : A string ID for the sequence\_elem at the gap's sink.
  4. *distance* : a 32-bit int which gives the gap's size. 
  5. *tags* : A map containing opt\_elems that provide annotations on the gap.

gap\_elem structs have the following methods:
  1. *to\_string\_2()* : returns a string of the GFA2 'G' line the struct represents.
  2. *to\_string()* : a wrapper for *to\_string\_2()*.

### path\_elem
path\_elem structs define a path (an ordered trace along sequences / nodes) in a GFA graph. They have the following members:
  1. *name* : A string name for the path.
  2. *segment\_names* : A vector containing the string names of the sequence\_elems on a path.
  3. *orientations* : A vector of booleans that holds whether a given element is oriented in the forward (true) or reverse (false) direction.
  4. *overlaps* : A vector of CIGAR string overlaps for each element in segment\_names.
  5. *opt\_fields* : A map containing opt\_elems providing annotations on the path.

They also have the following methods:
  1. *add\_ranked\_segment(rank, seg\_name, ori, overlap, opts)* : adds a single element to the path at the 1-based position *rank*.
  2. *to\_string\_1()* : returns a string representation of the path as a GFA1 'P' line.

## gfa\_builder.hpp
The gfa\_builder.hpp header file defines a set of functions used to build variation graphs
from a FASTA file and a VCF file. These functions support the `gfak build` command, but we
discourage the use of this interface outside of GFAKluge.
