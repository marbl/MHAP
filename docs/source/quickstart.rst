############
Quick Start
############

Running MHAP
-----------------

Running MHAP provides command-line documenation if you run it without parameters. Assuming you have followed the ?? instructions, you can run:
 
.. code-block:: bash

    $ java -jar mhap-0.1.jar

MHAP has two main usage modes, the main finds all overlaps between the input sequences. The second  only constructs an index which can be subsequently reused. 

Finding overlaps
-----------------

.. code-block:: bash

   $ java -Xmx32g -jar mhap-0.1.jar -s<fasta/dat from/self file> [-q<fasta/dat to file or directory>] [-f<kmer filter list, must be sorted>]

Both the -s and -q options can accept either FastA sequences or binary dat files (generated as described below). The -q option can accept either a file or a directory, in which case all FastA/dat files in the specified directory will be used. By default, only the sequences specified by -s are indexed and the sequences in -q are streamed against the constructed index. Since MHAP is written in Java, the memory usage can be high. Generally, 32GB of RAM is sufficient to index 20K sequences. If you have more sequences, you can partition your data and run MHAP on the partitions.

Constructing binary index
-----------------

.. code-block:: bash

   $ java -Xmx32g -jar mhap-0.1.jar -p<directory of fasta files> -q <output directory> [-f<kmer filter list, must be sorted>]

In this use case, files in the -p directory will be converted to binary dat files in the -q directory. Subsequent runs using the dat files (instead of FastA files) will be faster as the sequences no longer need to be indexed, only loaded into memory.

Options
-----------------
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

