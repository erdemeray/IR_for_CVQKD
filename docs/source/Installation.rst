.. _installation:

.. note::
   Python users should prefer installing via pip. Building from source is only required for development or C++ integration.

Installation Guide
==========================

The recommended way to install the library is via **pip** from PyPI:

.. code-block:: bash

    pip install information-reconciliation

After installation, you can directly import the library in Python:

.. code-block:: python

    import information_reconciliation

--------------------------------------------------
Alternative Installation Methods (Advanced Users)
--------------------------------------------------

If you need to build the library from source (e.g., for development or C++ usage), you can clone the repository.

**Clone the repository:**

.. code-block:: bash

    git clone https://github.com/erdemeray/IR_for_CVQKD.git

or

.. code-block:: bash

    git clone git@github.com:erdemeray/IR_for_CVQKD.git

++++++++++++++++++++++++++++++++++++++++++++
Using the Library in C++
++++++++++++++++++++++++++++++++++++++++++++

Prerequisites:
- CMake
- A C++ compiler supporting C++11 and OpenMP (GCC, Clang, MSVC)

Build steps:

.. code-block:: bash

    cd IR_lib/cpp
    mkdir build
    cd build
    cmake ..
    make

The compiled library (`information_reconciliation`) can then be linked in your C++ projects.

++++++++++++++++++++++++++++++++++++++++++++
Using the Library in Python (from source)
++++++++++++++++++++++++++++++++++++++++++++

If you prefer building from source instead of using PyPI:

.. code-block:: bash

    pip install .

Alternatively, you can use a Conda-based setup for full dependency isolation.

Using Conda
-----------------------------

1. Create environment:

.. code-block:: bash

    conda env create -f environment.yml

2. Activate environment:

.. code-block:: bash

    conda activate IR_lib

3. Build the project:

.. code-block:: bash

    mkdir build
    cd build
    cmake ..
    cmake --build . --target install

++++++++++++++++++++++++++++++++++++++++++++
Using Docker
++++++++++++++++++++++++++++++++++++++++++++

Build and run the Docker container:

.. code-block:: bash

    docker build -t IR_lib -f docker/Dockerfile .
    docker run --rm -it --mount "type=bind,src=${PWD},target=/developer" IR_lib

Then build the project inside the container as described above.

+++++++++++++++++++++++++++++++++++++++++++
Building the Documentation
+++++++++++++++++++++++++++++++++++++++++++

To build the documentation locally:

.. code-block:: bash

    pip install .[docs]
    cmake -S . -B build -DBUILD_DOCS=ON
    cmake --build build --target Sphinx
