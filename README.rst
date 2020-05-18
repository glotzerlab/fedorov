Fedorov
==================================================
A python package to initialize different crystal structures. Named after the Russian mathematician, crystallographer and mineralogist: `Evgraf Fedorov <https://en.wikipedia.org/wiki/Evgraf_Fedorov/>`_. This package provides simple ways to generate 590 known crystal structures defined in `Aflow <http://aflowlib.org/CrystalDatabase/>`_ or any user defined arbitrary crystal structures with proper space group and Wyckoff position information. The output of the crystal defined is provided through basis position: N by 3 numpy array with each row containing one particle position in the unit cell with N particles, and lattice vectors: 3 by 3 numpy array with each row vector describing the unit cell one unit cell dimension.

Installation
--------------------------------------------------

Install with pip
================

To install the package with the package manager pip_, execute

.. code:: bash

    $ pip install fedorov --user

To upgrade the package, simply execute the same command with the ``--upgrade`` option.

.. code:: bash

    $ pip install fedorov --user --upgrade

.. _pip: https://pip.pypa.io/en/stable/

Install from source
========================

Alternatively you can clone the `git repository <https://github.com/glotzerlab/fedorov>`_ and execute the ``setup.py`` script to install the package.

.. code:: bash

  git clone https://github.com/glotzerlab/fedorov.git
  cd fedorov
  python setup.py install --user

Documentation
--------------------------------------------------
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

Quickstart
--------------------------------------------------

This package provides simple ways to generate known or user defined crystal structures. For example, to generate any known crystal structures from this `list <https://github.com/glotzerlab/fedorov/blob/master/fedorov/crystal_data/Aflow_processed_data.csv>`_:

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

Testing
--------------------------------------------------

You can test this package by executing:

.. code-block:: bash

    python -m pytest tests/

within the repository root directory.

Authors
--------------------------------------------------
- Pengji Zhou (zhoupj@umich.edu)

Citation
--------------------------------------------------

The pre-defined crystal structures are obtained by `Aflow <http://aflowlib.org/CrystalDatabase/>`_

- \M. J. Mehl, D. Hicks, C. Toher, O. Levy, R. M. Hanson, G. L. W. Hart, and S. Curtarolo. The AFLOW Library of Crystallographic Prototypes: Part 1, Comp. Mat. Sci. 136, S1-S828 (2017). (`doi=10.1016/j.commatsci.2017.01.017 <http://doi.org/10.1016/j.commatsci.2017.01.017>`_)

- \D. Hicks, M. J. Mehl, E. Gossett, C. Toher, O. Levy, R. M. Hanson, G. L. W. Hart, and S. Curtarolo. The AFLOW Library of Crystallographic Prototypes: Part 2, Comp. Mat. Sci. 161, S1-S1011 (2019). (`doi=10.1016/j.commatsci.2018.10.043 <http://doi.org/10.1016/j.commatsci.2018.10.043/>`_)

The space group information are obtained from the `Bilbao Crystallographic Server <https://www.cryst.ehu.es/>`_ :

- \M. I. Aroyo, J. M. Perez-Mato, D. Orobengoa, E. Tasci, G. de la Flor, A. Kirov.
  "Crystallography online: Bilbao Crystallographic Server"
  Bulg. Chem. Commun. 43(2) 183-197 (2011).
  (`<http://bcc.bas.bg/BCC_Volumes/Volume_43_Number_2_2011/Volume_43_Number_2_2011_PDF/2011_43_2_1.pdf/>`_)

- \M. I. Aroyo, J. M. Perez-Mato, C. Capillas, E. Kroumova, S. Ivantchev, G. Madariaga, A. Kirov & H. Wondratschek.
  "Bilbao Crystallographic Server I: Databases and crystallographic computing programs"
  Z. Krist. 221, 1, 15-27 (2006). (`doi:10.1524/zkri.2006.221.1.15 <http://dx.doi.org/10.1524/zkri.2006.221.1.15/>`_)

- \M. I. Aroyo, A. Kirov, C. Capillas, J. M. Perez-Mato & H. Wondratschek.
  "Bilbao Crystallographic Server II: Representations of crystallographic point groups and space groups"
  Acta Cryst. A62, 115-128 (2006). (`doi:10.1107/S0108767305040286 <http://dx.doi.org/10.1107/S0108767305040286/>`_)
