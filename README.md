[![build](https://github.com/erdemeray/IR_for_CVQKD/actions/workflows/main.yml/badge.svg)](https://github.com/erdemeray/IR_for_CVQKD/actions/workflows/main.yml)
[![Documentation Status](https://readthedocs.org/projects/information-reconciliation-for-cv-qkd/badge/?version=latest)](https://information-reconciliation-for-cv-qkd.readthedocs.io/en/latest/?badge=latest)
[![PyPI version](https://img.shields.io/pypi/v/information-reconciliation?logo=pypi&logoColor=white)](https://pypi.org/project/information-reconciliation/)
[![Python versions](https://img.shields.io/pypi/pyversions/information-reconciliation?logo=python&logoColor=white)](https://pypi.org/project/information-reconciliation/)
[![Install with pip](https://img.shields.io/badge/install-pip%20install%20information--reconciliation-blue?logo=pypi&logoColor=white)](https://pypi.org/project/information-reconciliation/)



# Information Reconciliation Library for CV-QKD Systems

This C++ library provides a comprehensive set of tools for information reconciliation, specifically designed for continuous-variable quantum key distribution (CV-QKD) systems. The library also includes Python bindings, enabling seamless integration with Python and MATLAB environments.


## Installation

### Python (recommended)

Install the package directly from PyPI:

```bash
pip install information-reconciliation
```

After installation, you can import it in Python as:

```python
import information_reconciliation
```

### From source / C++

If you want to build the project from source, use it in C++, or set up a development environment, please refer to the full documentation.

## Documentation

For detailed information on installation, library usage, and implementation details, please refer to our comprehensive documentation: [Library Documentation](https://information-reconciliation-for-cv-qkd.readthedocs.io/).

## Research Article

This library is built upon the principles and methodologies detailed in our associated research article. For a deeper understanding of the library's foundations, we encourage you to review the publication:
[Research Article](https://arxiv.org/abs/2408.00569). 

## Licensing

This project is open-source and distributed under the terms of the GPL-3.0 License. For more details, please see the LICENSE file in the repository.

## Acknowledgements

We gratefully acknowledge the following contributions and support:

- Funding support from the German Federal Ministry of Education and Research (BMBF), grant number 16KISQ056 (DE-QOR).

- The C++ example provided in the library utilizes the [cxx_argp](https://github.com/pboettch/cxx_argp) library for command-line input parsing, which is available under the [LGPL3](https://www.gnu.org/licenses/lgpl-3.0.html) license. 


## Citation:

If you use this library in your research, please cite our work as follows:
```
E. E. Cil and L. Schmalen, "An open-source library for information reconciliation in continuous-variable QKD," *Proc. International Conference on Quantum Cryptography (QCRYPT)*, Vigo, Spain, Sep. 2024
```

For any questions or feedback, please open an issue in the repository or contact the maintainers directly.
