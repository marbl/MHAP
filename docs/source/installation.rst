############
Installation
############

Before your start
=================
MHAP requires a recent version of the `JVM <http://www.oracle.com/technetwork/java/javase/downloads/jre7-downloads-1880261.html>`_ (1.7u51+). JDK 1.6 or earlier will not work. If you would like to build the code from source, you need to have the `JDK <http://www.oracle.com/technetwork/java/javase/downloads/jdk7-downloads-1880260.html>`_ and the `ANT <http://ant.apache.org/>`_ build system available.

Prerequisites
==============
    * java (1.7u51+)
    * ant (1.8.2+)

Here is a list of currently supported Operating Systems:

1. Mac OSX (10.7 or newer)
2. Linux 64-bit (tested on CentOS, Fedora, RedHat, OpenSUSE and Ubuntu)

Installation
======================
Pre-compiled
-----------------

To download a pre-compiled tar run:

.. code-block:: bash

    $ wget https://github.com/marbl/MHAP/releases/download/v0.1/mhap-0.1.tar.gz

And if ``wget`` not available, you can use ``curl`` instead:

.. code-block:: bash

    $ curl -L https://github.com/marbl/MHAP/releases/download/v0.1/mhap-0.1.tar.gz > mhap-0.1.tar.gz

Then run

.. code-block:: bash

   $ tar xvzf mhap-0.1.tar.gz

Source
-----------------

To build the code from the release:

.. code-block:: bash

    $ wget https://github.com/marbl/MHAP/archive/v0.1.zip

If you see a certificate not trusted error, you can add the following option to wget:

.. code-block:: bash

    $ --no-check-certificate

And if ``wget`` not available, you can use ``curl`` instead:

.. code-block:: bash

    $ curl -L https://github.com/marbl/MHAP/archive/v0.1.zip > v0.1.zip

You can also browse the https://github.com/marbl/MHAP/tree/v0.1
and click on Downloads. 

Once downloaded, extract to unpack:

.. code-block:: bash

    $ unzip v0.1.zip

Change to MetAMOS directory:

.. code-block:: bash

    $ cd MHAP-0.1

Once inside the MetAMOS directory, run:

.. code-block:: bash

    $ ant

This will compile the program and create a target/mhap-0.1.jar file which you can use to run MHAP. 
