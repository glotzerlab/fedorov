########################################
Fedorov
########################################

A python package to initialize different crystal structures. Named after the Russian mathematician, crystallographer and mineralogist: `Evgraf Fedorov <https://en.wikipedia.org/wiki/Evgraf_Fedorov/>`_. This package provides methods to generate any user defined crystal structures in both 2D and 3D systems, and includes access to 590 known crystal structures defined in `Aflow <http://aflowlib.org/CrystalDatabase/>`_. Users can use this package to obtain the necessary data in form of numpy array to use for softwares such as `HOOMD-Blue <https://hoomd-blue.readthedocs.io/en/stable/index.html/>`_ to construct any systems with specific crystal structure for simulation and analysis. In addition, the package also provides access to all the 3D space group, 2D plane group as well as point group operations that allows user to apply different symmetry operations.

****************************************
Installation
****************************************

Install with pip
----------------------------------------

To install the package with the package manager pip_, execute

.. code:: bash

    $ pip install fedorov --user

To upgrade the package, simply execute the same command with the ``--upgrade`` option.

.. code:: bash

    $ pip install fedorov --user --upgrade

.. _pip: https://pip.pypa.io/en/stable/

Install from source
----------------------------------------

Alternatively you can clone the `git repository <https://github.com/glotzerlab/fedorov>`_ and execute the ``setup.py`` script to install the package.

.. code:: bash

  git clone https://github.com/glotzerlab/fedorov.git
  cd fedorov
  python setup.py install --user

****************************************
Documentation
****************************************

- Documentation link: Website link to be added

You can also build the documentation yourself with sphinx:

You can install sphinx using conda:

.. code-block:: bash

    conda install sphinx

or pip:

.. code-block:: bash

    pip install sphinx

To build the documentation, run the following commands in the source directory:

.. code-block:: bash

    cd doc
    make html
    # Then open _build/html/index.html

To build a PDF of the documentation (requires LaTeX and/or PDFLaTeX):

.. code-block:: bash

    cd doc
    make latexpdf
    # Then open _build/latex/fedorov.pdf

****************************************
Quickstart
****************************************

This package provides methods to generate known or user defined crystal structures. For example, to generate any known crystal structures from this `list <https://github.com/glotzerlab/fedorov/blob/master/fedorov/crystal_data/Aflow_processed_data.csv>`_:

.. code-block:: python

    import numpy as np
    from fedorov import SpaceGroup, Prototype, AflowPrototype
    from fedorov import convert_to_box
    # generate the exact prototype provided by Aflow, use prototype_index [0, 589]
    prototype_index = 5
    new_structure = AflowPrototype(prototype_index=prototype_index, print_info=True,
                                   set_type=True)
    basis_vectors, type_list = new_structure.get_basis_vectors()
    lattice_vectors = new_structure.get_lattice_vectors()
    Lx, Ly, Lz, xy, xz, yz = convert_to_box(lattice_vectors)

More example can be found `here <https://github.com/glotzerlab/fedorov/tree/master/demo>`_.

****************************************
Testing
****************************************

You can test this package by executing:

.. code-block:: bash

    python -m pytest tests/

within the repository root directory.

****************************************
Credits
****************************************
See `Credits <https://github.com/glotzerlab/fedorov/blob/master/Credits.rst>`_.
