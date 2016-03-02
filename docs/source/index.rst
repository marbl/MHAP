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

MHAP is included within the `Canu <http://canu.readthedocs.org/>`_ assembler. Canu can be downloaded `here <https://github.com/marbl/canu/releases>`_.

Contents:

.. toctree::
   :maxdepth: 2

   installation
   quickstart
   utilities
   contact
