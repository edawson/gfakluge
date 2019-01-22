---
title: 'GFAKluge: A C++ library and command line utilities for the Graphical Fragment Assembly formats'
tags:
  - GFA
  - C++
  - genome assembly
  - bioinformatics
authors:
 - name: Eric T. Dawson
   orcid:  0000-0001-5448-1653
   affiliation: "1, 2, 3"
 - name: Richard Durbin
   orcid: 0000-0002-9130-1006
   affiliation: "2, 3"
affiliations:
 - name: Division of Cancer Epidemiology and Genetics, National Cancer Institute, Rockville, MD, USA
   index: 1
 - name: Department of Genetics, University of Cambridge, Cambridge, UK
   index: 2
 - name: Wellcome Sanger Institute, Hinxton, UK
   index: 3
date: 17 September 2018
bibliography: paper.bib
---

# Summary
GFAKluge is a set of command line utilities and a C++ library for parsing and
manipulating the Graphical Fragment Assembly (GFA) format. 
Genome assembly algorithms often use graph structures
to represent relationships between reads during the assembly process. Such information
is typically thrown away when assemblies are converted to FASTA files of contig sequences. 
Previous attempts to convey graph information did not gain widespread acceptance because there were no standard representations that were easily parsed and extensively used.
The Graphical Fragment Assembly
(GFA) format was proposed as a way to encode the graph structure of an assembly in a human-readable
text format [@GFAOriginal]. 
GFA aims to provide a single format for interchange between software for assembly, scaffolding, assessment and
visualization. Such programs are often written in high-performance
programming languages such as C or C++. GFAKluge facilitates interprogram exchange by providing
a high-level C++ API for developers and a set of command line tools for users. We hope the availability of an open-source,
easily extensible API will encourage software developers to consider adding support for GFA to their
bioinformatics programs.  
*Homepage:* https://github.com/edawson/gfakluge
*License:* MIT

# Command Line Utilities
GFAKluge also provides a command line interface for working with GFA. This includes support for
common tasks on assemblies such as calculating assembly N50 or graph statistics. There are also methods for merging
assemblies, reformating files for readability, and converting between the GFA 1.0 and GFA 2.0 specifications. A tool for constructing basic variation graphs
from a FASTA file and a VCF file is also included.
Many other tools exist for manipulating the GFA formats [@GFA-SPEC], though only RGFA [@RGFA], GfaPy [@GfaPy] and ABySS2.0 [@ABySS2.0] are known to produce and consume both versions.
By allowing interconversion
between the compatible subsets of the formats, the `gfak convert` tool allows programs that usually can't communicate to share data
without changes to their code. We have used GFAKluge to convert GFA from TwoPaCo [@TwoPaCo] for visualization in Bandage [@Bandage], to calculate assembly
statistics from the Falcon assembler [@Falcon], and to extract FASTA from a `vg msga` assembly [@vg].

```
   # Convert GFA 2.0 from TwoPaCo to GFA 1.0 for ingestion by Bandage.  
   gfak convert -S 1.0 data/gfa_2.gfa  

   # Calculate assembly statistics  
   gfak stats -a data/gfa_2.gfa  

   # Extract FASTA entries from a GFA file  
   gfak extract data/gfa_2.gfa  

```

The full list of `gfak` commands follows:  
```
   convert: Convert between GFA 0.1 <-> 1.0 <-> 2.0
   diff:    Determine whether two GFA files have identical graphs
   extract: Convert the S lines of a GFA file to FASTA format.
   fillseq: Add sequences from a FASTA file to S lines.
   ids:     Coordinate the ID spaces of multiple GFA graphs.
   concat:  Merge GFA graphs (without ID collisions).
   sort:    Print a GFA file in HSLP / HSEFGUO order.
   stats:   Get assembly statistics (e.g. N50) for a GFA file.
   subset:  Extract the subgraph between two IDs in a graph.
   trim:    Remove elements from a GFA graph.
```

Examples of most commands are included in the [examples.md file](https://github.com/edawson/gfakluge/blob/master/examples.md).

# Integrating GFAKluge into an existing program
As an example of how to use the GFAKluge API, we briefly summarize its use in the variation graph toolkit [vg](https://github.com/vgteam/vg) [@vg].
vg creates bidirected sequence graphs from assemblies and population variation that can then be used for read mapping and variant calling. We incorporated
GFAKluge into vg to support input and output of GFA. Reading in a GFA file requires one line of code and is agnostic to
the GFA version used. Converting from GFA to vg's internal structures and vice versa requires approximately forty lines of code. Changing output from
GFA v1.0 to GFA v2.0 requires a single API call. This allows vg to take assemblies in GFA format from TwoPaCo and many other assembly algorithms.
The gfak command line tools can be used to calculate assembly graph statistics on graphs produced by vg. A full description of the developer API is available in
the [interface.md file](https://github.com/edawson/gfakluge/blob/master/interface.md).

