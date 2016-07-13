############
Quick Start
############

Running MHAP
-----------------

Running MHAP provides command-line documenation if you run it without parameters. Assuming you have followed the `installation instructions <installation.html>`_ instructions, you can run:
 
.. code-block:: bash

    $ java -jar mhap-2.1.jar

MHAP has two main usage modes, the main finds all overlaps between the input sequences. The second  only constructs an index which can be subsequently reused. 

Finding overlaps
-----------------

.. code-block:: bash

   $ java -Xmx32g -server -jar mhap-2.1.jar -s<fasta/dat from/self file> [-q<fasta/dat to file or directory>] [-f<kmer filter list, must be sorted>]

Both the -s and -q options can accept either FastA sequences or binary dat files (generated as described below). The -q option can accept either a file or a directory, in which case all FastA/dat files in the specified directory will be used. By default, only the sequences specified by -s are indexed and the sequences in -q are streamed against the constructed index. Since MHAP is written in Java, the memory usage can be high. Generally, 32GB of RAM is sufficient to index 40K sequences. If you have more sequences, you can partition your data and run MHAP on the partitions. You can also increase the memory MHAP is allowed to use by changing the Xmx parameter to a larger limit.

The optional -f flag provides a file of repetitive k-mers which should be biased against selected as min-mers. The file is a two-column tab-delimited input specifying the kmer and the fraction of total kmers the k-mer comprises. For example:

.. code-block:: bash

   $ head kmers.ignore
   464
   GGGGGGGGGGGGG	0.0005

means the k-mer GGGGGGGGGGG represents 0.05% of the k-mers in the dataset (so if there are 100,000 total k-mers, it occurs 50 times). The first line specifies the total number of k-mer entries in the file.

It is also possible to use the k-mer list as a positive selection as was used in `Carvalho et. al. <http://biorxiv.org/content/biorxiv/early/2016/05/14/053256.full.pdf>`_. Specify the k-mer list as above and the flag:

.. code-block:: bash
   --supress-noise 2

which will not allow any k-mer not in in the input file to be a minmer. The k-mers above --filter-threshold will be ignored as repeats.

.. code-block:: bash
   --supress-noise 1

will downweight any k-mer not in the input file to bias against its selection as a minmer. The k-mers above --filter-threshold will be downeighted as repeats.

Constructing binary index
-----------------

.. code-block:: bash

   $ java -Xmx32g -server -jar mhap-2.1.jar -p<directory of fasta files> -q <output directory> [-f<kmer filter list, must be sorted>]

In this use case, files in the -p directory will be converted to binary sketch files in the -q directory. Subsequent runs using these files (instead of FastA files) will be faster as the sequences no longer need to be sketched, only loaded into memory.

Output
-----------------
MHAP outputs overlaps in a format similar to BLASR's M4 format. Example output::

   [A ID] [B ID] [% error] [# shared min-mers] [0=A fwd, 1=A rc] [A start] [A end] [A length] [0=B fwd, 1=B rc] [B start] [B end] [B length]

An example of output from a small dataset is below::

   155 11 0.164156 206 0 69 1693 1704 0 1208 2831 5871
   155 15 0.157788 163 0 16 1041 1704 1 67 1088 2935
   155 27 0.185483 159 0 455 1678 1704 0 0 1225 1862

In this case sequence 155 overlaps 11, 15, and 27. The error percent is computed from the Jaccard estimate using `mash distance <http://www.biorxiv.org/content/early/2015/10/26/029827.abstract>`_. 

Options
-----------------
The full list of options is available via command-line help (--help or -h). Below is a list of commonly used options.

   -h                  Displays the help menu.
   --version           Displays the version.
   --settings          Set all unset parameters for the default settings. Same defaults are applied to Nanopore and Pacbio reads. 0) None, 1) Default, 2) Fast, 3) Sensitive.
   -k                		[int], k-mer size used for MinHashing. The k-mer size for second stage filter is seperate, and cannot be modified, default = 16
   -q                		Usage 1: The FASTA file of reads, or a directory of files, that will be compared to the set of reads in the box (see -s). Usage 2: The output directory for the binary formatted dat files, default="".
   -s                		Usage 1 only. The FASTA or binary dat file (see Usage 2) of reads that will be stored in a box, and that all subsequent reads will be compared to, default = ""
   -p                		Usage 2 only. The directory containing FASTA files that should be converted to binary format for storage, default = ""
   --num-hashes        [int], number of min-mers to be used in MinHashing, default = 512.
   --num-min-matches   [int], minimum # min-mer that must be shared before computing second stage filter. Any sequences below that value are considered non-overlapping, default = 3.
   --threshold         [double], the threshold cutoff for the second stage sort-merge filter. This is based on the mash distance computed from the Jaccard distance of k-mers (size given by ordered-kmer-size) in the overlapping regions, default = 0.78
   --max-shift         [double], region size to the left and right of the estimated overlap, as derived from the median shift and sequence length, where a k-mer matches are still considered valid. Second stage filter only, default = 0.2.
   --ordered-kmer-size  [int], The size of k-mers used in the ordered second stage filter, , default = 12.
   --ordered-sketch-size  [int], The sketch size for second stage filter, default = 1536.
   --filter-threshold  [double], the cutoff at which the k-mer in the k-mer filter file is considered repetitive. This value for a specific k-mer is specified in the second column in the filter file. If no filter file is provided, this option is ignored, default = 1.0E-5.
   --repeat-weight        [double] Repeat suppression strength for tf-idf weighing. <0.0 do unweighted MinHash (version 1.0), >=1.0 do only the tf weighing. To perform no idf weighting, do no supply -f option, default = 0.0
			--supress-noise 			    [int] 0) Does nothing, 1) completely removes any k-mers not specified in the filter file, 2) supresses k-mers not specified in the filter file, similar to repeats, default = 0
   -f                     [string], k-mer filter file used for filtering out highly repetative k-mers. Must be sorted in descending order of frequency (second column), default = "".
   --no-self           Do not compute the overlaps between sequences inside a box. Should be used when the to and from sequences are coming from different files, default = false.
   --no-tf             Do not perform the tf weighing, of the tf-idf weighing, default = false
			--num-threads       [int], number of threads to use for computation. Typically set to #cores, , default = 4.
			--store-full-id        Store full IDs as seen in FASTA file, rather than storing just the sequence position in the file. Some FASTA files have long IDS, slowing output of results. This options is ignored when using compressed file format, default = false.
