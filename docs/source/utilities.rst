############
Utilities
############

Using MHAP extras
-----------------

In addition to the main overlapping algorithm, MHAP indcludes several utilities for validating overlaps and simulating data.

Validating overlaps
-----------------

Assuming you have a mapping of sequences to a truth (such as a reference genome) in BLASR's M4 format, you can validate overlaps using MHAP's EstimateROC utility which will compute PPV/Sensitivity/Specificity:

.. code-block:: bash

   $ java -cp mhap-2.1.1.jar edu.umd.marbl.mhap.main.EstimateROC <reference mapping M4> <overlaps M4/MHAP> <fasta of sequences> [minimum overlap length to evaluate] [number of random trials] [use dynamic programming] [verbose] [minimum identity of overlap] [maximum different between expected overlap and reported] [load all overlaps]

The default minimum overlap length is 2000 and default number of trials is 10000. This will estimate sensitivity/specificity to within 1%. It can be increased at the expense of runtime. Specifying 0 will examine all possible N^2 overlap pairs. 

The dynamic programming flag (true/false) will check overlaps not present in the reference mapping by running a Smith-Watermann alignment to identify the overlap specified. This step requires the `SSW Library <https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library>`_ to be separately compiled and installed:

.. code-block:: bash

   $ wget https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library/archive/master.zip
   $ unzip master.gip && cd  Complete-Striped-Smith-Waterman-Library-master/src
   $ make all
   $ cd /full/path/to/mhap/target/lib
   $ ln -s /full/path/to/Complete-Striped-Smith-Waterman-Library-master/src/libsswjni.so

The verbose flag (true/false) enables logging to report true overlaps missing from the result and false-positives where no alignments could be found matching the required thresholds.

The minimum identity of the overlap (0.7 by default) is the lower bound for the sensitivity of an overlapper to evaluate. It is used to select matches to the reference that could be found by the overlapper. It is also used to threshold the minimum identity found by the Smith-Waterman alignment above.

The load all overlaps flag (true/false) will evaluate the specificity and PPV on all overlaps reported by the overlapper if enabled, not only those for good reads (where both reads were mapped to the reference in the truth set).

Simulating Data
-----------------

MHAP includes a tool to simulate sequencing data with random error as well as estimate Jaccard similarity for the simulated data.

.. code-block:: bash

   $ java -cp mhap-2.1.1.jar edu.umd.marbl.mhap.main.KmerStatSimulator <# sequences> <sequence length (bp)> <insertion error rate> <deletion error rate> <substitution error rate> [reference genome]

The error rates must be between 0 and 1 and are additive. Specifying 10% insertion, 2% deletion, and 1% substitution will result in sequences with a 13% error rate. If no reference sequence is given, completely random sequences are generated and errors added. Otherwise, random sequences are drawn from the reference and errors added. Errors are added randomly with no bias.

.. code-block:: bash

   $  java -cp mhap-2.1.1.jar edu.umd.marbl.mhap.main.KmerStatSimulator <# trials> <kmer size> <sequence length> <overlap length> <insertion error rate> <deletion error rate> <substitution error rate> [one-sided error] [reference genome] [kmer filter]

This usage will output a distribution of Jaccard similarity between a pair of overlapping sequences with the specified error rate (when using the specified k-mer size) and two random sequences of the same length. If no reference sequence is given, completely random sequences are generated and errors added, otherwise sequences are drawn from the reference. When one-sided error is specified (by typing true for the parameter), only one of the two sequences will have error simulated, matching a mapping of a noisy sequence to a reference. If a set of k-mers for filtering is given, they are excluded when computing Jaccard similarity, both between random and overlapping sequences.
