=================================================================
MinHash Alignment Process (MHAP): a probabilistic sequence overlap algorithm.
============================================================

=================
Overview
=================
MHAP (pronounced MAP) is a reference implementation of a probabilistic              
sequence overlapping algorithm. Designed to efficiently detect all overlaps
between noisy long-read sequence data. It efficiently estimates Jaccard similarity
by compressing sequences to their representative fingerprints composed on min-mers (minimum k-mer).

=================
Overview
=================
Usage 1 (direct execution): MHAP -s<fasta/dat from/self file> [-q<fasta/dat to file>] [-f<kmer filter list, must be sorted>]
Usage 2 (generate precomputed binaries): MHAP -p<directory of fasta files> -q <output directory> [-f<kmer filter list, must be sorted>]
Options: 
   -k [int merSize], default: 16
    --memory [do not store kmers in memory]
    --num-hashes [int # hashes], default: 1024
    --min-store-length [int # of minimum sequence length that is hashed], default: 0
    --threshold [int threshold for % matching minimums], default: 0.05
    --max-shift [double fraction of the overlap size where shift in k-mer match is still considered valid], default: 0.2
    --num-min-matches [int # hashes that maches before performing local alignment], default: 3
    --num-threads [int # threads to use for computation], default (2 x #cores): 16
    --subsequence-size [depricated, int size of maximum minhashed sequence], default: 100000
    --no-self [do not compute results to self], default: false
    --store-full-id [use full sequence id rather than order in file], default: false
    --threshold [int threshold for % matching minimums], default: 0.05
    --max-shift [int # max sequence shift allowed for a valid kmer relative to median value], default: 0.2
    
Contents:

.. toctree::
   :maxdepth: 2

   content/installation
   content/quickstart
   content/utilities
