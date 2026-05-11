#include "_piqtree.h"
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <iostream>
#include <string>
#include <vector>
#include "_piqtree.h"

namespace py = pybind11;

void checkError(char* errorStr) {
  if (errorStr && std::strlen(errorStr) > 0) {
    string msg(errorStr);
    iqtree_free(errorStr);
    throw std::runtime_error(msg);
  }
  if (errorStr)
    iqtree_free(errorStr);
}

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

      tmpStrings.push_back(item.cast<string>());
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

    tmpDoubles.assign(arr.data(), arr.data() + value.length);
    value.doubles = tmpDoubles.data();

    return true;
  }

  // Conversion from C++ to Python
  static handle cast(DoubleArray src, return_value_policy, handle) {
    throw std::runtime_error("Unsupported operation");
  }

 private:
  vector<double> tmpDoubles;
};

template <>
struct type_caster<IntegerResult> {
 public:
  PYBIND11_TYPE_CASTER(IntegerResult, _("IntegerResult"));

  // Conversion from Python to C++
  bool load(handle /* src */, bool /* convert */) {
    throw std::runtime_error("Unsupported operation");
  }

  // Conversion from C++ to Python
  static handle cast(const IntegerResult& src, return_value_policy, handle) {
    checkError(src.errorStr);

    return py::int_(src.value).release();
  }
};

template <>
struct type_caster<StringResult> {
 public:
  // Indicate that this caster only supports conversion from C++ to Python
  PYBIND11_TYPE_CASTER(StringResult, _("StringResult"));

  // Reject Python to C++ conversion
  bool load(handle src, bool) {
    throw std::runtime_error("Unsupported operation");
  }

  // Conversion from C++ to Python
  static handle cast(const StringResult& src, return_value_policy, handle) {
    checkError(src.errorStr);

    PyObject* py_str = PyUnicode_FromString(src.value);
    if (!py_str)
      throw error_already_set();

    iqtree_free(src.value);

    return handle(py_str);
  }
};

template <>
struct type_caster<DoubleArrayResult> {
 public:
  PYBIND11_TYPE_CASTER(DoubleArrayResult, _("DoubleArrayResult"));

  // Conversion from Python to C++
  bool load(handle src, bool) {
    throw std::runtime_error("Unsupported operation");
  }

  // Conversion from C++ to Python
  static handle cast(DoubleArrayResult src, return_value_policy, handle) {
    checkError(src.errorStr);

    auto result = py::array_t<double>(src.length);

    std::memcpy(result.mutable_data(), src.value, src.length * sizeof(double));
    iqtree_free(src.value);

    return result.release();
  }
};

template <>
struct type_caster<HessianResult> {
 public:
  PYBIND11_TYPE_CASTER(HessianResult, _("HessianResult"));

  bool load(handle src, bool) {
    throw std::runtime_error("Unsupported operation");
  }

  // Returns a dict { 'tree': str, 'branch_lengths': ndarray(n,),
  //                  'gradient': ndarray(n,), 'hessian': ndarray(n, n) }.
  static handle cast(HessianResult src, return_value_policy, handle) {
    // checkError raises if errorStr is non-empty; in that case we still need
    // to release any payload pointers IQ-TREE may have allocated before the
    // failure, to avoid leaking. checkError frees src.errorStr itself.
    if (src.errorStr && std::strlen(src.errorStr) > 0) {
      if (src.tree) iqtree_free(src.tree);
      if (src.branch_lengths) iqtree_free(src.branch_lengths);
      if (src.gradient) iqtree_free(src.gradient);
      if (src.hessian) iqtree_free(src.hessian);
      checkError(src.errorStr);  // throws
    }
    checkError(src.errorStr);

    py::ssize_t n = static_cast<py::ssize_t>(src.n_branches);

    auto branch_lengths = py::array_t<double>(n);
    std::memcpy(branch_lengths.mutable_data(), src.branch_lengths,
                static_cast<size_t>(n) * sizeof(double));
    iqtree_free(src.branch_lengths);

    auto gradient = py::array_t<double>(n);
    std::memcpy(gradient.mutable_data(), src.gradient,
                static_cast<size_t>(n) * sizeof(double));
    iqtree_free(src.gradient);

    auto hessian = py::array_t<double>({n, n});
    std::memcpy(hessian.mutable_data(), src.hessian,
                static_cast<size_t>(n) * static_cast<size_t>(n) * sizeof(double));
    iqtree_free(src.hessian);

    PyObject* py_tree = PyUnicode_FromString(src.tree ? src.tree : "");
    if (!py_tree)
      throw error_already_set();
    iqtree_free(src.tree);

    py::dict result;
    result["tree"] = py::reinterpret_steal<py::object>(py_tree);
    result["branch_lengths"] = branch_lengths;
    result["gradient"] = gradient;
    result["hessian"] = hessian;
    return result.release();
  }
};
}  // namespace detail
}  // namespace PYBIND11_NAMESPACE

int mine() {
  return 42;
}

PYBIND11_MODULE(_piqtree, m) {
  m.doc() = "_piqtree - Linking IQ-TREE to Python!";

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
  m.def("iq_fit_tree_hessian", &fit_tree_hessian,
        "Fit a tree and additionally return the log-likelihood gradient "
        "(first derivatives w.r.t. branch lengths) and the full Hessian "
        "matrix at the fitted point. Returns a dict with keys 'tree', "
        "'branch_lengths', 'gradient', 'hessian'. Non-partitioned models only.");
  m.def("iq_model_finder", &modelfinder,
        "Find optimal model for an alignment.");
  m.def("iq_jc_distances", &build_distmatrix,
        "Construct pairwise distance matrix for alignment.");
  m.def("iq_nj_tree", &build_njtree,
        "Build neighbour-joining tree from distance matrix.");
  m.def("iq_consensus_tree", &consensus_tree,
        "Compute a consensus tree from a sequence of trees.");
  m.def("iq_simulate_alignment", &simulate_alignment,
        "Simulate an alignment with AliSim.");
  m.def("mine", &mine, "The meaning of life, the universe (and everything)!");
}
