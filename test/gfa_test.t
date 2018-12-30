#!/bin/bash


BASH_TAP_ROOT=./bash-tap
. ./bash-tap/bash-tap-bootstrap

plan tests 8

## Importing a GFA file in GFA1 and converting it to GFA2
is $(./gfak convert -S 2.0 data/v1.gfa | md5sum | awk '{ print $1 }') "268e075f19c7600304b51247b11e5f0f" "Importing a GFA file in GFA1 and converting it to GFA2"


## Importing a GFA file in GFA2 and converting it to GFA1
is $(./gfak convert -S 1.0 data/gfa_2.gfa | md5sum | awk '{ print $1 }' ) "d7bb881a8880850acb2977efa28c7979" "Importing a GFA file in GFA2 and converting it to GFA1"

## Sorting a GFA1 file
is $(./gfak sort data/test.gfa | md5sum | awk '{ print $1 }') "6dd44a9a0cc7308c7d6b92e8f0d9e648" "Sorting a GFA1 file"

## Sorting a GFA2 file
is $(./gfak sort data/gfa_2.gfa | md5sum | awk '{ print $1 }') "fa3b92296d3a23f9db99e611815788d4" "Sorting a GFA2 file"

## Extracting FASTA records from a GFA file
is $(./gfak extract data/gfa_2.gfa | sort | md5sum | awk '{ print $1 }') "43bbe8fee3f67fd90b90ee885ddb15e3" "Extracting FASTA records from a GFA file"

## Replacing sequence placeholders with FASTA records
is $(./gfak fillseq -f data/no_seqs.fa data/no_seqs.gfa | md5sum | awk '{ print $1 }') "caaf91eac390521d68d56bad57f7b3b3" "A GFA file can have sequence placeholders replaced with seqs from a FASTA file."

## Bumping the ID space of a graph.
is $(./gfak ids -s 9:9:9 data/test.gfa | sort | md5sum | awk '{ print $1 }') $(cat data/re_id.gfa | sort | md5sum | awk '{ print $1 }') "gfak ids can bump the ID space of the graph."

## Merging two GFA files
is $(./gfak merge -S 2.0 data/test.gfa data/gfa_2.gfa | md5sum | awk '{ print $1 }') "ca3b52673b63de931cd64a50669e7147" "Two graphs can be merged."
