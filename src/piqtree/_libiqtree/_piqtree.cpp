#include "_piqtree.h"
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/string_view.h>
#include <nanobind/stl/vector.h>
#include <cstring>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace nb = nanobind;
using namespace nb::literals;

using StringList = std::vector<std::string_view>;
using DoubleArrayIn =
    nb::ndarray<const double, nb::ndim<1>, nb::c_contig, nb::device::cpu>;
using DoubleArrayOut = nb::ndarray<nb::numpy, double, nb::ndim<1>>;

void checkError(char* errorStr) {
  if (errorStr && std::strlen(errorStr) > 0) {
    string msg(errorStr);
    iqtree_free(errorStr);
    throw std::runtime_error(msg);
  }
  if (errorStr)
    iqtree_free(errorStr);
}

/*
 * Arguments passed to IQ-TREE. These point into the Python objects, which the
 * converted arguments keep alive for the duration of the call.
 */
class StringArrayArg {
 public:
  explicit StringArrayArg(const StringList& strings) {
    cStrs.reserve(strings.size());
    for (const auto& str : strings) {
      cStrs.push_back(str.data());
    }
    array.strings = cStrs.data();
    array.length = cStrs.size();
  }

  StringArray& get() { return array; }

 private:
  std::vector<const char*> cStrs;
  StringArray array;
};

class DoubleArrayArg {
 public:
  explicit DoubleArrayArg(const DoubleArrayIn& arr) {
    array.doubles = const_cast<double*>(arr.data());
    array.length = arr.size();
  }

  DoubleArray& get() { return array; }

 private:
  DoubleArray array;
};

/*
 * Results returned by IQ-TREE. These raise a RuntimeError if IQ-TREE reported
 * an error, and take ownership of memory allocated by IQ-TREE.
 */
int unwrap(const IntegerResult& result) {
  checkError(result.errorStr);

  return result.value;
}

std::string unwrap(const StringResult& result) {
  checkError(result.errorStr);

  std::string value(result.value);
  iqtree_free(result.value);

  return value;
}

DoubleArrayOut unwrap(const DoubleArrayResult& result) {
  checkError(result.errorStr);

  // The numpy array uses IQ-TREE's memory, which is freed with the array.
  nb::capsule owner(result.value, [](void* p) noexcept { iqtree_free(p); });

  return DoubleArrayOut(result.value, {result.length}, owner);
}

/*
 * Call an IQ-TREE function and convert its result for Python. IQ-TREE can
 * throw exceptions not derived from std::exception, which nanobind would raise
 * as a SystemError. These are raised as a RuntimeError like all other errors.
 */
template <typename Func>
auto callIqtree(Func&& func) {
  try {
    return unwrap(func());
  } catch (const std::exception&) {
    throw;
  } catch (...) {
    throw std::runtime_error("Caught an unknown exception!");
  }
}

int mine() {
  return 42;
}

NB_MODULE(_piqtree, m) {
  m.doc() = "_piqtree - Linking IQ-TREE to Python!";

  m.attr("__iqtree_version__") = callIqtree(version);

  m.def(
      "iq_robinson_fould",
      [](const char* tree1, const char* tree2) {
        return callIqtree([&] { return robinson_fould(tree1, tree2); });
      },
      "tree1"_a, "tree2"_a,
      "Calculates the robinson fould distance between two trees");
  m.def(
      "iq_random_tree",
      [](int num_taxa, const char* tree_gen_mode, int num_trees,
         int rand_seed) {
        return callIqtree([&] {
          return random_tree(num_taxa, tree_gen_mode, num_trees, rand_seed);
        });
      },
      "num_taxa"_a, "tree_gen_mode"_a, "num_trees"_a, "rand_seed"_a,
      "Generates a set of random phylogenetic trees. tree_gen_mode "
      "allows:\"YULE_HARDING\", \"UNIFORM\", \"CATERPILLAR\", \"BALANCED\", "
      "\"BIRTH_DEATH\", \"STAR_TREE\".");
  m.def(
      "iq_build_tree",
      [](const StringList& names, const StringList& seqs, const char* model,
         int rand_seed, int bootstrap_replicates, int num_threads,
         const char* other_options) {
        StringArrayArg names_arg(names), seqs_arg(seqs);
        return callIqtree([&] {
          return build_tree(names_arg.get(), seqs_arg.get(), model, rand_seed,
                            bootstrap_replicates, num_threads, other_options);
        });
      },
      "names"_a, "seqs"_a, "model"_a, "rand_seed"_a, "bootstrap_replicates"_a,
      "num_threads"_a, "other_options"_a,
      "Perform phylogenetic analysis on the input alignment (in string "
      "format). With estimation of the best topology.");
  m.def(
      "iq_fit_tree",
      [](const StringList& names, const StringList& seqs, const char* model,
         const char* tree, bool bl_fixed, int rand_seed, int num_threads,
         const char* other_options) {
        StringArrayArg names_arg(names), seqs_arg(seqs);
        return callIqtree([&] {
          return fit_tree(names_arg.get(), seqs_arg.get(), model, tree,
                          bl_fixed, rand_seed, num_threads, other_options);
        });
      },
      "names"_a, "seqs"_a, "model"_a, "tree"_a, "bl_fixed"_a, "rand_seed"_a,
      "num_threads"_a, "other_options"_a,
      "Perform phylogenetic analysis on the input alignment (in string "
      "format). With restriction to the input toplogy.");
  m.def(
      "iq_model_finder",
      [](const StringList& names, const StringList& seqs, int rand_seed,
         const char* model_set, const char* freq_set, const char* rate_set,
         int num_threads, const char* other_options) {
        StringArrayArg names_arg(names), seqs_arg(seqs);
        return callIqtree([&] {
          return modelfinder(names_arg.get(), seqs_arg.get(), rand_seed,
                             model_set, freq_set, rate_set, num_threads,
                             other_options);
        });
      },
      "names"_a, "seqs"_a, "rand_seed"_a, "model_set"_a, "freq_set"_a,
      "rate_set"_a, "num_threads"_a, "other_options"_a,
      "Find optimal model for an alignment.");
  m.def(
      "iq_jc_distances",
      [](const StringList& names, const StringList& seqs, int num_threads) {
        StringArrayArg names_arg(names), seqs_arg(seqs);
        return callIqtree([&] {
          return build_distmatrix(names_arg.get(), seqs_arg.get(), num_threads);
        });
      },
      "names"_a, "seqs"_a, "num_threads"_a,
      "Construct pairwise distance matrix for alignment.");
  m.def(
      "iq_nj_tree",
      [](const StringList& names, const DoubleArrayIn& distances) {
        StringArrayArg names_arg(names);
        DoubleArrayArg distances_arg(distances);
        return callIqtree(
            [&] { return build_njtree(names_arg.get(), distances_arg.get()); });
      },
      "names"_a, "distances"_a,
      "Build neighbour-joining tree from distance matrix.");
  m.def(
      "iq_consensus_tree",
      [](const StringList& trees, double min_support) {
        StringArrayArg trees_arg(trees);
        return callIqtree(
            [&] { return consensus_tree(trees_arg.get(), min_support); });
      },
      "trees"_a, "min_support"_a,
      "Compute a consensus tree from a sequence of trees.");
  m.def(
      "iq_simulate_alignment",
      [](const char* tree, const char* subst_model, int rand_seed,
         const char* partition_info, const char* partition_type, int seq_length,
         double insertion_rate, double deletion_rate, const char* root_seq,
         int num_threads, const char* insertion_size_distribution,
         const char* deletion_size_distribution, int population_size) {
        return callIqtree([&] {
          return simulate_alignment(
              tree, subst_model, rand_seed, partition_info, partition_type,
              seq_length, insertion_rate, deletion_rate, root_seq, num_threads,
              insertion_size_distribution, deletion_size_distribution,
              population_size);
        });
      },
      "tree"_a, "subst_model"_a, "rand_seed"_a, "partition_info"_a,
      "partition_type"_a, "seq_length"_a, "insertion_rate"_a, "deletion_rate"_a,
      "root_seq"_a, "num_threads"_a, "insertion_size_distribution"_a,
      "deletion_size_distribution"_a, "population_size"_a,
      "Simulate an alignment with AliSim.");
  m.def("mine", &mine, "The meaning of life, the universe (and everything)!");
}
