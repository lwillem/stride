#include <pybind11/stl.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

void init_age_group(py::module &);
void init_mdp(py::module &);
void init_vaccine_types(py::module &);
//void init_health(py::module &);

// Export the library pylibstride as a module
PYBIND11_MODULE(pylibstride, m) {
    m.doc() = "STRIDE module";

    init_age_group(m);
    init_mdp(m);
    init_vaccine_types(m);
//    init_health(m);
}
