There are other utilities that might be nice for gfa manipulation:
    [X] Merge: merges GFA files, respecting graph structure. Probably requires an ID run beforehand.  
    Can be accomplished by a simple `cat` and sort, but might be nice to implement natively.  
    [X] IDs: coordinates the ids of a set of GFA files.  
    [] Diff: Disregard IDs and sequence and check if two graphs are isomorphic (should just be an  
    order of output lines?  
    [X] Stats: print various stats about a GFA file such as the number of nodes, number of links,  
    total length, alignments, containments, etc  
    [] Verify: make sure there are no dangling links or containments that lack seq_elems  
    [] subgraph: extract a subgraph within an interval of two seg_ids.  
