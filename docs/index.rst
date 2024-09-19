.. nucleoseeker documentation master file, created by
   sphinx-quickstart on Wed Jun 26 09:57:59 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to NucleoSeeker's documentation!

NucleoSeeker - A tool for precision filtering of RNA structures to enhance Deep learning predictions

.. toctree::
   :maxdepth: 2
   :caption: Contents:

Getting Started
===============

To get started with NucleoSeeker, follow the following steps:

Currently we only support Unix based systems including MacOS.

**Setup Steps**
**Get Clustal Omega ready**

* Instructions to setup clustal-omega can be found http://www.clustal.org/omega/INSTALL
* Clustal omega version supported `1.2.4`

.. code-block:: bash

      wget http://www.clustal.org/omega/clustal-omega-1.2.4.tar.gz
      tar zxf clustal-omega-1.2.4.tar.gz
      cd clustal-omega-1.2.4
      ./configure --prefix /your/install/path
      make
      make check                 # optional: run automated tests
      make install               # optional: install Infernal programs, man pages

      # or use this

      sudo apt-get install clustalo



**Get Emboss ready (Optional)**

* **NOTE** - *Emboss is very slow, unless you are experimenting we don't recommend using it. Clustal Omega should be sufficient for most use cases.*
* For setting up Emboss, please read http://emboss.open-bio.org/html/adm/ch01s01.html
* Emboss version supported `6.6.0`

**Get Infernal ready**

* For infernal follow instructions http://eddylab.org/infernal
* Infernal version supported `1.1.5`

.. code-block:: bash

      wget http://eddylab.org/software/infernal/infernal.tar.gz
      tar zxf infernal.tar.gz
      cd infernal-1.1.5
      ./configure --prefix /your/install/path
      make
      make check                 # optional: run automated tests
      make install               # optional: install Infernal programs, man pages

      # or use this

      sudo apt-get install infernal infernal-doc

**Get Rfam.cm file ready**

* To use this tool, you need to provide Rfam covariance model which is available for download at https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz. It also needs to be modified using `cmpress` command from `Infernal` tool (mentioned above). If you don't have it then use the code below - 

.. code-block:: bash

      cd nucleoseeker
      mkdir -p rfam
      wget https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz -O Rfam.cm.gz
      gunzip Rfam.cm.gz
      cmpress Rfam.cm

Installation
============

To install NucleoSeeker, you can use the following steps:

.. code-block:: bash

      git clone __repo_url__
      cd nucleoseeker
      pip install -r requirements.txt

Usage
=====
To generate a new dataset using NucleoSeeker, the following command can be used:

.. code-block:: bash

      export DATA_PATH=/path/to/data # the dataset will be saved here
      python3 src/dataset_creator.py --dataset_name test_dataset --rfam_cm_path your/rfam/path --exptl_method "X-RAY DIFFRACTION" --resolution 3.6 --year_range 2024 --save 1 --dend 500


After using this command a directory with the name `test` will be created in the `DATA_PATH` directory with the following subdirectories:

Directories:

.. code-block:: bash
   
      DATA_PATH
      ├── test_dataset
      │   ├── files
      │   ├── sequences
      ├──clean_tblout.tblout
      ├──cmscan.out
      ├──combined.fasta
      ├──fam_pdb_chain.csv
      ├──final.fasta
      ├──raw_experimental_RNA_0_500.csv
      ├──sequence_identity_mat_clustal.csv
      ├──tblout.tblout


* This tool generates various files mostly at each level of filter. The first file that is generated is the `raw` csv file which contains the raw data from the PDB database.

* Then the `combined.fasta` file is generated which contains sequences used in sequence identity calculation by `Clustal Omega` and `Emboss`. This file is obtained after applying `StructureLevelFilter` and `PDBFilter` on the raw data.

* The `sequence_identity_mat_clustal(emboss).csv` file contains the sequence identity matrix obtained from `Clustal Omega` and `Emboss` tools.

* The `final.fasta` file contains the final sequences in fasta format, these are the final sequences and if you don't want to analyse families then this is the final output.

* After this `cmscan.out`, `tblout.tblout`, `clean_tblout.tblout` files are generated which are the output of `Infernal` tool. The `fam_pdb_chain.csv` file contains the mapping of the family and the PDB chain.

* The `fam_pdb_chain.csv` is obtained after family search by `Infernal` tool. This is your final output if you want to analyse families.

* The `files` directory contains the dataframe and list for structures at each level of filter.

* The `sequences` directory sequences for all the final structures in individual fasta files.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
