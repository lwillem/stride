// Reference to header file(s) needed
#include "../../../cpp/mdp/Vaccines.h"
// Some extra (external) header files of the project required for successfully building the pybind library
//

#include <pybind11/stl.h>
#include <pybind11/pybind11.h>
namespace py = pybind11;


void init_vaccine_types(py::module &m) {
    // VaccineType enumeration
    py::enum_<stride::VaccineType>(m, "VaccineType")
            .value("noVaccine", stride::VaccineType::noVaccine)
            .value("mRNA", stride::VaccineType::mRNA)
            .value("adeno", stride::VaccineType::adeno);
    // All vaccine types
    m.attr("AllVaccineTypes") = stride::AllVaccineTypes;
    // Vaccine Properties
    py::class_<stride::VaccineProperties, std::shared_ptr<stride::VaccineProperties>>(m, "VaccineProperties", "The properties of a vaccine.");
    // Constant Vaccine Properties
    py::class_<stride::ConstantVaccineProperties, std::shared_ptr<stride::ConstantVaccineProperties>, stride::VaccineProperties>(m, "ConstantVaccineProperties", "The properties of a constant vaccine.")
            .def(py::init<std::string&, double, double, double>(),
                    py::arg("id"), py::arg("ve_susceptible"), py::arg("ve_infectiousness"), py::arg("ve_severe"));
    // Linear Vaccine Properties
    py::class_<stride::LinearVaccineProperties, std::shared_ptr<stride::LinearVaccineProperties>, stride::VaccineProperties>(m, "LinearVaccineProperties", "The properties of a linear vaccine.")
            .def(py::init<std::string&, double, double, double, unsigned short>(),
                 py::arg("id"), py::arg("ve_susceptible"), py::arg("ve_infectiousness"), py::arg("ve_severe"),
                 py::arg("max_ve_day"));
}
