Installation
============

This section provides instructions on how to install JWSToolKit.

The first step before installing JWSToolKit on your machine is to make sure you have created a python environment 
separate from your Python system. By installing JWSToolKit, dependencies will be installed that can break and modify your existing system dependencies. 


Requirements
^^^^^^^^^^^^^

Before installing JWSToolKit, you should have the following dependencies installed on your machine: 

* numpy - For mathematical processing
* matplotlib - For data visualisation
* scipy - For mathematical processing
* astropy.io - To manipulate data in .fits format
* tqdm - Progress bar in the terminal
* photutils - To handle telescope observations

.. warning:: 
    In particular, pay close attention to the photutils package. 
    To run JWSToolKit routines correctly, you need to have the latest version of photutils 
    installed on your machine, as well as a version of Python later than 3.11.


Installation with *pip*
^^^^^^^^^^^^^^^^^^^^^^^^^

To install the package with pip, use the following command:  

.. code-block:: console

    pip install JWSToolKit

If you want to install a specific version of the package, use the command line: 

.. code-block:: console

    pip install JWSToolKit==1.0.4

Finally, if the package is already installed on your machine but you wish to update it, use one of the commands: 

.. code-block:: console

    pip install --upgrade JWSToolKit
    pip install --upgrade JWSToolKit==1.0.4


Installation with *conda*
^^^^^^^^^^^^^^^^^^^^^^^^^^^

When using a conda environment, you should install the package via the *conda* command: 

.. code-block:: console

    conda install delabrov::jwstoolkit

To update the package: 

.. code-block:: console

    conda update jwstoolkit
    conda update jwstoolkit=1.0.4       # to download a specific version 

In this situation, it is advisable to create a Python environment separate from the default system environment. To do this, enter the following command:

.. code-block:: console

    conda create -n name_of_the_environment 

Once created, you can launch the environment with the command: 

.. code-block:: console

    conda activate name_of_the_environment 


