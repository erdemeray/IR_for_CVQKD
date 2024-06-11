# CV-QKD Information Reconciliation Library

This library, written in C++, offers tools for information reconciliation tailored for continuous-variable quantum key distribution (CV-QKD) systems. It also includes Python bindings for use within Python environments.

## Article

The link to the article will be added soon.

## Quick Start

Follow these steps to set up the project on your local system.

### Requirements

Before you begin, ensure you have the following installed:
- CMake (version 3.18 or newer)
- g++ compiler
- Python (only if using the Python bindings)
- pybind11 (only if using the Python bindings)

### Compilation

Execute these commands from the root directory of the project to compile the library:

```sh
mkdir build
cd build
cmake ..
make
```

### Documentation

The complete documentation for the library is available at the following link: [Library Documentation](https://information-reconciliation-for-cv-qkd.readthedocs.io/).

## License
The project is open-source, under the GPL-3.0 License.

## Acknowledgements

- Funding for this project was provided by the German Federal Ministry of Education and Research (BMBF), grant number 16KISQ056 (DE-QOR).

- The C++ test suite utilizes the [cxx_argp](https://github.com/pboettch/cxx_argp) library for command-line input parsing, which is available under the [LGPL3](https://www.gnu.org/licenses/lgpl-3.0.html) license. 


## Cite this work as:
E. E. Cil and L. Schmalen, "Rate-adaptive protograph-based raptor-like LDPC code for continuous-variable quantum key distribution," *Proc. Advanced Photonic Congress: Signal Processing in Photonic Communications (SPPCom)*, Québec City, Canada, Jul. 2024
