Installation Guide
==========================

++++++++++++++++++++++++++++++++ 
Prerequisites
++++++++++++++++++++++++++++++++    

Before installing the library, ensure you have the following dependencies installed:

- **g++ compiler** (version 9.3.0 or higher)
- **CMake** (version 3.18 or higher)
- **OpenMP** (version 5.0 or higher)

For Python integration:

- **Python** (version 3.9 or higher)
- **Pybind11** (version 2.6.2 or higher)

++++++++++++++++++++++++++++++++
Installing the Library
++++++++++++++++++++++++++++++++

For C++ Usage:
----------------

1. **Clone the repository** to your local machine.
2. Navigate to the `cpp` directory within the cloned repository.
3. Execute `CMake` to configure and generate the build files, then compile the library.

For Python Integration:
------------------------

1. **Clone the repository** just like for C++.
2. Change directory to the `python` folder.
3. Use `CMake` to set up the build environment and compile the library.
4. The compiled library can be used in Python by importing the `information_reconciliation` module found in the build directory.

++++++++++++++++++++++++++++++++
Windows-Specific Installation Tips
++++++++++++++++++++++++++++++++

- You can install `mingw-w64 <https://www.mingw-w64.org/>`_ to get the `g++` compiler and OpenMP library.
- When building on **Windows**, if you face issues with locating `PyBind11`, specify the path to `PyBind11`'s CMake file to aid in its discovery.
- If you get a `dll not found` error when trying to import the library in Python, add the following directory to the DLL search path as follows:
  
    .. code-block:: python
    
        import os 
        os.add_dll_directory('C:\\Windows\\SysWOW64')
