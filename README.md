# trmap
Streaming a large GTF (or BED) across a reference transcript database 
(GTF or BED) stored in memory as an interval tree, 
checking for transcript overlaps or matching.

The streaming input GFF/GTF MUST be well-formed -- i.e. exons MUST 
be grouped together by transcript ID and immediately follow
their parent feature if present. (for BED-12 this is a given
due to exons being embedded in the same single line).

