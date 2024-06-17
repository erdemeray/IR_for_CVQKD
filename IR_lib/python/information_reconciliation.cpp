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