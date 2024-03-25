// Reference to header file(s) needed
#include "../../../cpp/mdp/AgeGroup.h"
// Some extra (external) header files of the project required for successfully building the pybind library
//

#include <pybind11/stl.h>
#include <pybind11/pybind11.h>
namespace py = pybind11;


void init_age_group(py::module &m) {
    // AgeGroup enumeration
    py::enum_<stride::AgeGroup>(m, "AgeGroup")
            .value("children", stride::AgeGroup::children)
            .value("youngsters", stride::AgeGroup::youngsters)
            .value("young_adults", stride::AgeGroup::young_adults)
            .value("adults", stride::AgeGroup::adults)
            .value("elderly", stride::AgeGroup::elderly);
    // All age groups
    m.attr("AllAgeGroups") = stride::AllAgeGroups;
    // Get AgeGroup given an age
    m.def("GetAgeGroup", &stride::GetAgeGroup, py::arg("age"), "Get the age group for the given age");

    // ChildlessAgeGroup enumeration
    py::enum_<stride::ChildlessAgeGroup>(m, "ChildlessAgeGroup")
            .value("young_adults", stride::ChildlessAgeGroup::young_adults_c)
            .value("adults", stride::ChildlessAgeGroup::adults_c)
            .value("elderly", stride::ChildlessAgeGroup::elderly_c);
    // All age groups
    m.attr("AllChildlessAgeGroups") = stride::AllChildlessAgeGroups;
    // Get AgeGroup given an age
    m.def("GetChildlessAgeGroup", &stride::GetChildlessAgeGroup, py::arg("age"), "Get the age group for the given age");
}
