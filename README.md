# trmap
This utility allows streaming of a large GFF (or BED) file in order to compare it to a
reference transcript database (GTF/BED) which is stored in memory as an interval tree.
This allows for fast checking for transcript overlaps and classification of their relationship with reference transcripts.
"Class codes" like those assigned by gffcompare (see http://ccb.jhu.edu/software/stringtie/gffcompare.shtml/) 
are also provided in the output.

The streaming input GFF/GTF/BED must be well-formed -- i.e. exons MUST 
be grouped together by transcript ID and immediately follow
their parent feature if present. (for BED-12 this is a given
due to exons being embedded in the same single line).

