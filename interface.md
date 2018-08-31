The GFAKluge developer API
----------------------------

## Introduction
The GFAKluge.hpp header file defines an object-orientated API for managing graphs in GFA format. It also defines structs for each of the GFA line types.

## GFAKluge class and functions
The GFAKluge class defines a collection of structs representing a GFA file.

## Adding elements to a GFAKluge object

## Retrieving elements from a GFAKluge object

## GFA2 vs. GFA1

## Data Structures

### sequence\_elem
The sequence\_elem struct defines a node in a GFA graph. It has the following members:
  1. *name* : a string identifier for the sequence.
  2. *sequence* : the base pairs the node contains. "*" may also be used as a placeholder when sequences are held in an external FASTA.
  3. *length* : A uint64_t containing the length of *sequence*.
  4. *opt_fields* : a vector of opt\_elems (key:type:value triples) which describe annotations or tags on this sequence element.
  
  sequence\_elem structs have the following methods:
  1. to_string_1 : 
  2. to_string_2 : 
  3. as_fasta_record : 

### edge\_elem

### group\_elem

### fragment\_elem

### gap\_elem

### path\_elem

## gfak\_build.hpp

## Consistent functions across structs
### to_string()