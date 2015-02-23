############
Quick Start
############

Running MHAP
-----------------

Running MHAP provides command-line documenation if you run it without parameters. Assuming you have followed the `installation instructions <installation.html>`_ instructions, you can run:
 
.. code-block:: bash

    $ java -jar mhap-1.5b1.jar

MHAP has two main usage modes, the main finds all overlaps between the input sequences. The second  only constructs an index which can be subsequently reused. 

Finding overlaps
-----------------

.. code-block:: bash

   $ java -Xmx32g -server -jar mhap-1.5b1.jar -s<fasta/dat from/self file> [-q<fasta/dat to file or directory>] [-f<kmer filter list, must be sorted>]

Both the -s and -q options can accept either FastA sequences or binary dat files (generated as described below). The -q option can accept either a file or a directory, in which case all FastA/dat files in the specified directory will be used. By default, only the sequences specified by -s are indexed and the sequences in -q are streamed against the constructed index. Since MHAP is written in Java, the memory usage can be high. Generally, 32GB of RAM is sufficient to index 20K sequences. If you have more sequences, you can partition your data and run MHAP on the partitions. You can also increase the memory MHAP is allowed to use by changing the Xmx parameter to a larger limit.

The optional -f flag provides a file of repetitive k-mers which should not be selected as min-mers. The file is a two-column tab-delimited input specifying the kmer and the fraction of total kmers the k-mer comprises. For example:

.. code-block:: bash

   $ head kmers.ignore
   GGGGGGGGGGGGG	0.0005

means the k-mer GGGGGGGGGGG represents 0.05% of the k-mers in the dataset (so if there are 100,000 total k-mers, it occurs 50 times).

Constructing binary index
-----------------

.. code-block:: bash

   $ java -Xmx32g -server -jar mhap-1.5v1.jar -p<directory of fasta files> -q <output directory> [-f<kmer filter list, must be sorted>]

In this use case, files in the -p directory will be converted to binary dat files in the -q directory. Subsequent runs using the dat files (instead of FastA files) will be faster as the sequences no longer need to be indexed, only loaded into memory.

Output
-----------------
MHAP outputs overlaps in a format similar to BLASR's M4 format. Example output::

   [A ID] [B ID] [Jaccard score] [# shared min-mers] [0=A fwd, 1=A rc] [A start] [A end] [A length] [0=B fwd, 1=B rc] [B start] [B end] [B length]

An example of output from a small dataset is below::

   155 11 87.83225 206 0 69 1693 1704 0 1208 2831 5871
   155 15 85.08692 163 0 16 1041 1704 1 67 1088 2935
   155 27 87.11507 159 0 455 1678 1704 0 0 1225 1862

In this case sequence 155 overlaps 11, 15, and 27.

Options
-----------------
The full list of options is available via command-line help (--help or -h). Below is a list of commonly used options.

   -k  [int]  K-mer size, default=16
   --num-hashes  [int]  Sketch size, higher=more sensitive but more memory usage and runtime, default=1024
   --num-min-matches  [int]  The number of hashes that maches before performing local alignment, default=3
   --pacbio_fast  [boolean]  Set all the parameters for the PacBio fast setting. This is the current best guidance, and could change at any time without warning, default = false
   --pacbio_sensitive  [boolean]  Set all the parameters for the PacBio sensitive settings. This is the current best guidance, and could change at any time without warning, default = false
   --min-store-length  [int length (in bp)]  The minimum sequence length to index. Sequences shorter than this are ignored in the index, default=0
   --threshold  [int]   The threshold for percentage of matching min-mers for a hit to be considered significant. Lowering will output more overlaps but increase false positives, higher will reduce overlaps but remove false positives, default=0.04
   --filter-threshold  [double]  The cutoff at which the k-mer in the k-mer filter file is considered repetitive. This value for a specific k-mer is specified in the second column in the filter file. If no filter file is provided, this option is ignored, default = 1.0E-5
   --max-shift  [double]  The fraction of the overlap size by which the overlap sizes in two sequences may differ, default=0.2
   --num-threads  [int]  The number of threads to use for computation, default (2 x #cores on system)
   --no-self  Do not compute self-matches for sequences in the -s file, default=false
   --store-full-id  Output full sequence ID from the input FastA file. Otherwise, the output is the position of the sequence in the file (i.e. first sequence gets ID=1, second gets ID=2, and so on), default=false

