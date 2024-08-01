.. _installation:

Installation Guide
==========================

To obtain the library, clone the repository to your local machine by executing the following command in your terminal:

**Clone the repository** to your local machine by running the following command in your terminal:

    .. code-block:: bash

        git clone https://github.com/erdemeray/IR_for_CVQKD.git

    or

    .. code-block:: bash

        git clone git@github.com:erdemeray/IR_for_CVQKD.git

There are multiple methods to utilize the library, depending on the programming language you intend to use. The installation steps may vary accordingly.

* If you plan to use the library in C++, you can compile it using CMake and integrate it into your C++ projects. Ensure that CMake and a C++ compiler are installed on your system. Python is not required for using the library in C++.
* If you plan to use the library in Python and MATLAB, configure your system so that CMake can locate the Python installation. The recommended method is to use a Conda environment or a Docker container to install all dependencies in an isolated environment.

++++++++++++++++++++++++++++++++++++++++++++
Using the Library in C++
++++++++++++++++++++++++++++++++++++++++++++

Prerequisites include CMake and a C++ compiler. You can download CMake from [here](https://cmake.org/download/). Use any C++ compiler that supports the C++11 standard and the OpenMP library, such as GCC, Clang, or MSVC.

After installing the prerequisites, you can compile the library and integrate it into your C++ projects. To compile the library, navigate to the `IR_lib/cpp` directory within the cloned repository. Then, create a directory named `build` and change your directory to this folder. Execute the following commands:

.. code-block:: bash

    cd IR_lib/cpp
    mkdir build
    cd build

Now, execute `CMake` to configure and generate the build files, and compile the library by running the following commands:

.. code-block:: bash

    cmake ..
    make

The compiled library can be utilized in C++ projects by linking the `information_reconciliation` library found in the `build` directory. Example usage of the library in C++ is provided in the `IR_lib/cpp/tests` directory.

++++++++++++++++++++++++++++++++++++++++++++
Using the Library in Python and MATLAB
++++++++++++++++++++++++++++++++++++++++++++

Using Conda
-----------------------------

The simplest method to install the library is to create an environment in Conda and compile it. For this, ensure that Conda is installed on your system. If Conda is not installed, you can download it from [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

First, create a Conda environment and then compile the C++ code with Python bindings, allowing it to be used in Python by importing the `information_reconciliation` module.

Steps:

1. **Create a Conda Environment** Change directory to the `conda` folder within the cloned repository. Open the `environment.yml` file in a text editor and comment/uncomment the lines based on your operating system. Save the file.

Then, run the following command to create a Conda environment with the required dependencies:

    .. code-block:: bash

        conda env create -f environment.yml

Now, you have all the dependencies installed in the Conda environment named `IR_lib`.

2. **Activate the Conda Environment** by running the following command:

    .. code-block:: bash

        conda activate IR_lib

3. **Navigate to the `IR_lib` directory** and create a directory named `build`, then change your directory to this folder. Execute the following commands:

    .. code-block:: bash

        cd IR_lib
        mkdir build
        cd build

4. **Compile the Library** by running the following commands:

    .. code-block:: bash

        cmake ..
        make

The compiled library can be used in Python/MATLAB by importing the `information_reconciliation` module found in the `build` directory. Examples in Python are provided in the `IR_lib/tests` directory.

Using Docker
-----------------------------

Another method to use the library is to utilize the Docker container. Use the provided Dockerfile in the `docker` folder to build the Docker image and run the container. Ensure that Docker is installed on your system. If Docker is not installed, you can download it from [here](https://docs.docker.com/get-docker/).

To use the Docker container, build the image by running the following command in the root directory of the cloned repository:

.. code-block:: bash

    docker build -t IR_lib -f docker/Dockerfile .

After building the image, run the container by executing the following command:

.. code-block:: bash

    docker run --rm -it --mount "type=bind,src=${PWD},target=/developer" IR_lib

Now, you can compile the library and use it in the running container. First, connect to the running container. Then, navigate to the `IR_lib` directory and open the `CMakeLists.txt` file. Comment out the lines between 29 and 34, which are related to the Conda environment. The lines should look like this:

.. code-block:: cmake

    # Detect conda environment - COMMENT IF YOU DON'T USE CONDA
    # if(DEFINED ENV{CONDA_PREFIX})
    #   set(CMAKE_PREFIX_PATH $ENV{CONDA_PREFIX})
    #   message(STATUS "Conda environment detected: $ENV{CONDA_PREFIX}")
    # else()
    #   message(WARNING "Conda environment not detected. You may need to activate your environment.")
    # endif()

Then, follow steps 3 and 4 in the previous section to compile the library. The compiled library can be used in Python/MATLAB by importing the `information_reconciliation` module found in the `build` directory.
