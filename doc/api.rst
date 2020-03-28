.. _api:

API Reference
=============

This is the API for the **Fedorov** package.

Class for crystal initialization
--------------------------------

.. currentmodule:: fedorov

.. autoclass:: SpaceGroup
    :members:

.. autoclass:: Prototype
    :members:

.. autoclass:: AflowPrototype
    :show-inheritance:
    :members:

Class for unit cell
-------------------

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

Some methods for crystal initialization
---------------------------------------

.. currentmodule:: fedorov

.. autofunction:: wrap

.. autofunction:: convert_to_box

.. autofunction:: convert_to_vectors

.. autofunction:: get_volumn

.. autofunction:: fractional_to_cartesian

.. autofunction:: translate_to_vector
