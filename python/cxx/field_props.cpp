#include <string>

#include <opm/input/eclipse/EclipseState/Grid/FieldPropsManager.hpp>

#include <pybind11/stl.h>

#include "export.hpp"
#include "converters.hpp"

#include <python/cxx/OpmCommonPythonDoc.hpp>

namespace {

    bool contains( const FieldPropsManager& manager, const std::string& kw) {
        if (manager.has_int(kw))
            return true;
        if (manager.has_double(kw))
            return true;

        return false;
    }

    py::array_t<double> get_double_array(const FieldPropsManager& m, const std::string& kw) {
        if (m.has_double(kw))
            return convert::numpy_array( m.get_double(kw) );
        else
            throw std::invalid_argument("Keyword '" + kw + "'is not of type double.");
    }

    py::array_t<int> get_int_array(const FieldPropsManager& m, const std::string& kw) {
        if (m.has_int(kw))
            return convert::numpy_array( m.get_int(kw) );
        else
            throw std::invalid_argument("Keyword '" + kw + "'is not of type int.");
    }


    py::array get_array(const FieldPropsManager& m, const std::string& kw) {
        if (m.has_double(kw))
            return convert::numpy_array(m.get_double(kw));

        if (m.has_int(kw))
            return convert::numpy_array(m.get_int(kw));

        throw std::invalid_argument("No such keyword: " + kw);
    }


}


void python::common::export_FieldProperties(py::module& module) {

    using namespace Opm::Common::DocStrings;

    py::class_< FieldPropsManager >( module, "FieldProperties", FieldProperties_docstring)
    .def( "__contains__", &contains, py::arg("keyword"), FieldProperties_contains_docstring)
    .def("__getitem__", &get_array, py::arg("keyword"), FieldProperties_getitem_docstring)
    .def( "get_double_array",  &get_double_array, py::arg("keyword"), FieldProperties_get_double_array_docstring)
    .def( "get_int_array",  &get_int_array, py::arg("keyword"), FieldProperties_get_int_array_docstring)
    ;

}
