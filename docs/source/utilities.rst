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

   $ java -cp mhap-1.6.jar edu.umd.marbl.mhap.main.EstimateROC <reference mapping M4> <overlaps M4/MHAP> <fasta of sequences> [minimum overlap length to evaluate] [number of random trials] [use dynamic programming] [verbose]

The default minimum overlap length is 2000 and default number of trials is 10000. This will estimate sensitivity/specificity to within 1%. It can be increased at the expense of runtime. Specifying 0 will examine all possible N^2 overlap pairs. If the dynamic programming is turned on (by typing true for the parameter), overlaps not present in the reference mapping will be confirmed if a Smith-Watermann alignment can identify the overlap specified. 

Simulating Data
-----------------

MHAP includes a tool to simulate sequencing data with random error as well as estimate Jaccard similarity for the simulated data.

.. code-block:: bash

   $ java -cp mhap-1.6.jar edu.umd.marbl.mhap.main.KmerStatSimulator <# sequences> <sequence length (bp)> <insertion error rate> <deletion error rate> <substitution error rate> [reference genome]

The error rates must be between 0 and 1 and are additive. Specifying 10% insertion, 2% deletion, and 1% substitution will result in sequences with a 13% error rate. If no reference sequence is given, completely random sequences are generated and errors added. Otherwise, random sequences are drawn from the reference and errors added. Errors are added randomly with no bias.

.. code-block:: bash

   $  java -cp mhap-1.6.jar edu.umd.marbl.mhap.main.KmerStatSimulator <# trials> <kmer size> <sequence length> <overlap length> <insertion error rate> <deletion error rate> <substitution error rate> [one-sided error] [reference genome] [kmer filter]

This usage will output a distribution of Jaccard similarity between a pair of overlapping sequences with the specified error rate (when using the specified k-mer size) and two random sequences of the same length. If no reference sequence is given, completely random sequences are generated and errors added, otherwise sequences are drawn from the reference. When one-sided error is specified (by typing true for the parameter), only one of the two sequences will have error simulated, matching a mapping of a noisy sequence to a reference. If a set of k-mers for filtering is given, they are excluded when computing Jaccard similarity, both between random and overlapping sequences.
