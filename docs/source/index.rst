=============================================================================
MinHash Alignment Process (MHAP): a probabilistic sequence overlap algorithm.
=============================================================================

=================
Overview
=================
MHAP (pronounced MAP) is a reference implementation of a probabilistic              
sequence overlapping algorithm. Designed to efficiently detect all overlaps
between noisy long-read sequence data. It efficiently estimates Jaccard similarity
by compressing sequences to their representative fingerprints composed on min-mers (minimum k-mer).

MHAP is included within the Celera Assembler `PBcR <http://wgs-assembler.sourceforge.net/wiki/index.php?title=PBcR>`_ pipeline. The Celera Assembler can be downloaded `here <https://sourceforge.net/projects/wgs-assembler/files/wgs-assembler/wgs-8.3/>`_.

Contents:

.. toctree::
   :maxdepth: 2

   installation
   quickstart
   utilities
   contact
