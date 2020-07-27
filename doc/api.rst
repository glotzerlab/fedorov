.. _api:

API Reference
=============

This is the API for the **Fedorov** package.

Classes for 3D crystal initialization
-------------------------------------------------
This section contains three classes that allow user to initialize a 3D crystal structure with
different user input.

.. currentmodule:: fedorov

.. autoclass:: SpaceGroup
    :members:

.. autoclass:: Prototype
    :members:

.. autoclass:: AflowPrototype
    :show-inheritance:
    :members:

    The list of crystal structures provided from Aflow are summarized below:

    .. csv-table::
       :file: ../fedorov/crystal_data/Aflow_processed_data.csv
       :widths: 20, 20, 20, 20, 20, 20, 50, 20, 50
       :header-rows: 1

Classes for 3D unit cell
-------------------------------------------------

.. currentmodule:: fedorov

.. autoclass:: Triclinic
    :members:

.. autoclass:: Monoclinic
    :members:

.. autoclass:: Orthorhombic
    :members:

.. autoclass:: Tetragonal
    :members:

.. autoclass:: Hexagonal
    :members:

.. autoclass:: Rhombohedral
    :members:

.. autoclass:: Cubic
    :members:

Class for 2D crystal initialization
-------------------------------------------------
This section contains one class that allows user to initialize a 2D crystal structure with different
user input.

.. currentmodule:: fedorov

.. autoclass:: PlaneGroup
   :members:

Classes for 2D unit cell
-------------------------------------------------

.. currentmodule:: fedorov

.. autoclass:: Oblique2D
   :members:

.. autoclass:: Rectangular2D
   :members:

.. autoclass:: Hexagonal2D
   :members:

.. autoclass:: Square2D
   :members:

Class for Point group symmetry operations
-------------------------------------------------
This section contains one class that allows user to obtain all point group symmetry operations.

.. currentmodule:: fedorov

.. autoclass:: PointGroup
  :members:

Some methods for crystal initialization
-------------------------------------------------

.. currentmodule:: fedorov

.. autofunction:: wrap

.. autofunction:: convert_to_box

.. autofunction:: convert_to_vectors

.. autofunction:: get_volumn

.. autofunction:: fractional_to_cartesian

.. autofunction:: translate_to_vector

.. autofunction:: translate_to_vector_2D
