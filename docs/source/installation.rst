############
Installation
############

Before your start
=================
MHAP requires a recent version of the `JVM <http://www.oracle.com/technetwork/java/javase/downloads/jre7-downloads-1880261.html>`_ (1.8u6+). JDK 1.7 or earlier will not work. If you would like to build the code from source, you need to have the `JDK <http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html>`_ and the `Maven <https://maven.apache.org>`_ build system available.

Prerequisites
==============
    * java (1.8u6+)
    * maven (3.0+)

Here is a list of currently supported Operating Systems:

1. Mac OSX (10.7 or newer)
2. Linux 64-bit (tested on CentOS, Fedora, RedHat, OpenSUSE and Ubuntu)
3. Windows (XP or newer)

Installation
======================
Pre-compiled
-----------------

The pre-compiled version is recommended to users who want to run MHAP, without doing development. To download a pre-compiled tar run:

.. code-block:: bash

    $ wget https://github.com/marbl/MHAP/releases/download/v2.0/mhap-2.0.tar.gz

And if ``wget`` not available, you can use ``curl`` instead:

.. code-block:: bash

    $ curl -L https://github.com/marbl/MHAP/releases/download/v2.0/mhap-2.0.tar.gz > mhap-2.0.tar.gz

Then run

.. code-block:: bash

   $ tar xvzf mhap-2.0.tar.gz

Source
-----------------

To build the code from the release:

.. code-block:: bash

    $ wget https://github.com/marbl/MHAP/archive/v2.0.zip

If you see a certificate not trusted error, you can add the following option to wget:

.. code-block:: bash

    $ --no-check-certificate

And if ``wget`` not available, you can use ``curl`` instead:

.. code-block:: bash

    $ curl -L https://github.com/marbl/MHAP/archive/v2.0.zip > v2.0.zip

You can also browse the https://github.com/marbl/MHAP/tree/v2.0
and click on Downloads. 

Once downloaded, extract to unpack:

.. code-block:: bash

    $ unzip v2.0.zip

Change to MetAMOS directory:

.. code-block:: bash

    $ cd MHAP-2.0

Once inside the MetAMOS directory, run:

.. code-block:: bash

    $ maven

This will compile the program and create a target/mhap-2.0.jar file which you can use to run MHAP. The quick-start instructions assume you are in the target directory when running the program. You can also use the target/mhap-2.0.tar file to copy MHAP to a different system or directory. If you would like to run the `validation utilties <utilities.html>`_ you must also download and build the `SSW Library <https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library>`_. Follow the instructions on the `utilities <utilities.html>`_ page.
