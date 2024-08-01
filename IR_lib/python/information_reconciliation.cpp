/********************************************************************
 * Project Name: Information Reconciliation Library for CV-QKD
 * File Name: infromation_reconciliation.cpp
 * Description:
 *    This file contains the implementation of the Python bindings for the information reconciliation library.
 * Author: Erdem Eray Cil
 * License: GPL-3.0 License
 * Revision History:
 *    11/06/2024 - v0.1 - First pre-release version
 *    26/06/2024 - v0.2 - Added template definitions for Python bindings
 ********************************************************************/

#include <pybind11/pybind11.h>

namespace py = pybind11;

void init_statistics(py::module &);
void init_reconciliation(py::module &);

namespace mcl{
    PYBIND11_MODULE(information_reconciliation, m) {
    m.doc() = "C++ information_reconciliation wrapped for Python use.";

    init_statistics(m);
    init_reconciliation(m);

}}