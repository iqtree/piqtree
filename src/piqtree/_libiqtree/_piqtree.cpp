#include "_piqtree.h"
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

namespace py = pybind11;

namespace PYBIND11_NAMESPACE {
namespace detail {
template <>
struct type_caster<StringArray> {
 public:
  PYBIND11_TYPE_CASTER(StringArray, const_name("StringArray"));

  // Conversion from Python to C++
  bool load(handle src, bool) {
    /* Extract PyObject from handle */
    PyObject* source = src.ptr();
    if (!py::isinstance<py::sequence>(source)) {
      return false;
    }

    auto seq = reinterpret_borrow<py::sequence>(src);
    value.length = seq.size();

    tmpStrings.reserve(value.length);
    tmpCStrs.reserve(value.length);

    for (size_t i = 0; i < seq.size(); ++i) {
      auto item = seq[i];
      if (!py::isinstance<py::str>(item)) {
        return false;
      }

      tmpStrings.push_back(item.cast<std::string>());
      tmpCStrs.push_back(tmpStrings[i].c_str());
    }

    value.strings = tmpCStrs.data();

    return true;
  }

  // Conversion from C++ to Python
  static handle cast(StringArray src, return_value_policy, handle) {
    throw std::runtime_error("Unsupported operation");
  }

 private:
  vector<string> tmpStrings;
  vector<const char*> tmpCStrs;
};

template <>
struct type_caster<DoubleArray> {
 public:
  PYBIND11_TYPE_CASTER(DoubleArray, _("DoubleArray"));

  // Conversion from Python to C++
  bool load(handle src, bool) {
    if (!py::isinstance<py::array_t<double>>(src)) {
      return false;  // Only accept numpy arrays of float64
    }

    auto arr = py::cast<py::array_t<double>>(src);
    if (arr.ndim() != 1) {
      return false;  // Only accept 1D arrays
    }

    value.length = arr.size();
    value.doubles = new double[value.length];
    std::memcpy(value.doubles, arr.data(), value.length * sizeof(double));
    return true;
  }

  // Conversion from C++ to Python
  static handle cast(DoubleArray src, return_value_policy, handle) {
    auto result = py::array_t<double>(src.length);
    std::memcpy(result.mutable_data(), src.doubles,
                src.length * sizeof(double));
    return result.release();
  }
};
}  // namespace detail
}  // namespace PYBIND11_NAMESPACE

int mine() {
  return 42;
}

PYBIND11_MODULE(_piqtree, m) {
  m.doc() = "piqtree - Unlock the Power of IQ-TREE 2 with Python!";

  m.attr("__iqtree_version__") = version();

  m.def("iq_robinson_fould", &robinson_fould,
        "Calculates the robinson fould distance between two trees");
  m.def("iq_random_tree", &random_tree,
        "Generates a set of random phylogenetic trees. tree_gen_mode "
        "allows:\"YULE_HARDING\", \"UNIFORM\", \"CATERPILLAR\", \"BALANCED\", "
        "\"BIRTH_DEATH\", \"STAR_TREE\".");
  m.def("iq_build_tree", &build_tree,
        "Perform phylogenetic analysis on the input alignment (in string "
        "format). With estimation of the best topology.");
  m.def("iq_fit_tree", &fit_tree,
        "Perform phylogenetic analysis on the input alignment (in string "
        "format). With restriction to the input toplogy.");
  m.def("iq_model_finder", &modelfinder,
        "Find optimal model for an alignment.");
  m.def("iq_jc_distances", &build_distmatrix,
        "Construct pairwise distance matrix for alignment.");
  m.def("iq_nj_tree", &build_njtree,
        "Build neighbour-joining tree from distance matrix.");
  m.def("mine", &mine, "The meaning of life, the universe (and everything)!");
}
