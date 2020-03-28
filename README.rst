Fedorov
==================================================
A python package to initialize different crystal structures

Resources
--------------------------------------------------
- `Installation guide <installation.rst>`_
- Documentation: Website link to be added

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

Testing
--------------------------------------------------

You can test this package by executing:

.. code-block:: bash

    python -m pytest tests/

within the repository root directory.

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
