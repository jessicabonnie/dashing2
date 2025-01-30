// First, include the sketch headers
#include "sketch/sketch.h"
#include "sketch/dist.h"
#include "sketch/hll.h"

// Now include the rest
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "fastxsketch.h"
#include "counter.h"
#include <cstring>

namespace py = pybind11;

// Define the hasher type before use
using Hasher = sketch::hash::CEHasher;

// Wrapper for FinalRMinHash class
class PyRMinHash {
    sketch::mh::RangeMinHash<uint64_t> rmh_;

public:
    PyRMinHash(size_t sketch_size = 1024, uint64_t seed = 0) 
        : rmh_(sketch_size) {}

    void add_sequence(const std::string& seq, size_t kmer_size = 31) {
        if (seq.empty()) {
            throw std::invalid_argument("Sequence cannot be empty");
        }
        if (kmer_size < 1 || kmer_size > seq.length()) {
            throw std::invalid_argument("Invalid k-mer size");
        }
        bns::Encoder<uint64_t> enc(kmer_size);
        const char* data = seq.data();
        for(size_t i = 0; i + kmer_size <= seq.length(); ++i) {
            // Use the rolling hash function
            uint64_t hash = 0;
            for(size_t j = 0; j < kmer_size; ++j) {
                hash = hash * 31 + (uint8_t)data[i + j];
            }
            rmh_.add(hash);
        }
    }

    double jaccard_index(const PyRMinHash& other) const {
        auto s1 = rmh_.finalize();
        auto s2 = other.rmh_.finalize();
        return s1.jaccard_index(s2);
    }

    size_t get_size() const {
        return rmh_.size();
    }

    py::array_t<uint64_t> get_mins() const {
        auto finalized = rmh_.finalize();
        auto result = py::array_t<uint64_t>(finalized.first.size());
        auto r = result.mutable_unchecked<1>();
        for(size_t i = 0; i < finalized.first.size(); ++i) {
            r(i) = finalized.first[i];
        }
        return result;
    }

    void merge(const PyRMinHash& other) {
        for(const auto& x : other.rmh_.finalize().first) {
            rmh_.add(x);
        }
    }

    void clear() {
        rmh_.clear();
    }
};

// Fix the sketch_fastx function
py::dict sketch_fastx(const std::string& path, 
                     size_t sketch_size = 1024,
                     size_t kmer_size = 31,
                     bool canonicalize = true,
                     size_t window_size = 0) {
    
    // Create a vector of paths
    std::vector<std::string> paths{path};
    
    // Create a result object
    dashing2::FastxSketchingResult result;
    
    // Create options with a default value and set member variables
    dashing2::Dashing2Options options(kmer_size, window_size);  // Pass k and w to constructor
    options.sketchsize() = sketch_size;
    options.canonicalize(canonicalize);
    
    // Call fastx2sketch with result reference, options, paths and output path
    dashing2::fastx2sketch(result, options, paths, path);

    // Convert result to Python dictionary
    py::dict ret;
    ret["names"] = py::cast(result.names_);
    
    // Convert sketches to numpy array
    std::vector<ssize_t> shape;
    shape.push_back(static_cast<ssize_t>(result.signatures_.size() / sketch_size));  // Changed to signatures_
    shape.push_back(static_cast<ssize_t>(sketch_size));
    
    auto sketches = py::array_t<uint64_t>(shape);
    auto r = sketches.mutable_unchecked<2>();
    
    for(size_t i = 0; i < result.signatures_.size() / sketch_size; ++i) {  // Changed to signatures_
        for(size_t j = 0; j < sketch_size; ++j) {
            r(i, j) = result.signatures_[i * sketch_size + j];  // Changed to signatures_
        }
    }
    ret["sketches"] = sketches;
    
    return ret;
}

// Fix the compute_distances function
py::array_t<double> compute_distances(py::array_t<uint64_t> sketches) {
    auto buf = sketches.request();
    if (buf.ndim != 2) 
        throw std::runtime_error("Number of dimensions must be 2");
    
    size_t n_sketches = buf.shape[0];
    size_t sketch_size = buf.shape[1];
    
    std::vector<ssize_t> shape = {static_cast<ssize_t>(n_sketches), 
                                 static_cast<ssize_t>(n_sketches)};
    auto distances = py::array_t<double>(shape);
    auto dist = distances.mutable_unchecked<2>();
    
    uint64_t* ptr = static_cast<uint64_t*>(buf.ptr);
    
    for(size_t i = 0; i < n_sketches; ++i) {
        dist(i, i) = 0.0;  // Distance to self is 0
        for(size_t j = i + 1; j < n_sketches; ++j) {
            // Create temporary RangeMinHash objects
            sketch::mh::RangeMinHash<uint64_t> rmh1(sketch_size), rmh2(sketch_size);
            
            // Add values to RangeMinHash
            for(size_t k = 0; k < sketch_size; ++k) {
                rmh1.add(ptr[i * sketch_size + k]);
                rmh2.add(ptr[j * sketch_size + k]);
            }
            
            // Convert to FinalRMinHash
            auto s1 = rmh1.finalize();
            auto s2 = rmh2.finalize();
            
            double jac = s1.jaccard_index(s2);
            double dist_val = -std::log(jac);
            
            dist(i, j) = dist_val;
            dist(j, i) = dist_val;
        }
    }
    
    return distances;
}

// Add this new wrapper class after PyRMinHash class
class PyHLL {
    sketch::hll::hll_t hll_;

public:
    PyHLL(int p = 14) : hll_(p) {
        if (p < 4 || p > 30) {
            throw std::invalid_argument("Precision parameter p must be between 4 and 30");
        }
    }

    void add_sequence(const std::string& seq, size_t kmer_size = 31) {
        if (seq.empty()) {
            throw std::invalid_argument("Sequence cannot be empty");
        }
        if (kmer_size < 1 || kmer_size > seq.length()) {
            throw std::invalid_argument("Invalid k-mer size");
        }
        bns::Encoder<uint64_t> enc(kmer_size);
        const char* data = seq.data();
        for(size_t i = 0; i + kmer_size <= seq.length(); ++i) {
            // Use the same rolling hash function as in PyRMinHash
            uint64_t hash = 0;
            for(size_t j = 0; j < kmer_size; ++j) {
                hash = hash * 31 + (uint8_t)data[i + j];
            }
            hll_.add(hash);
        }
    }

    double get_estimate() const {
        return hll_.report();
    }

    void merge(const PyHLL& other) {
        hll_ += other.hll_;
    }

    void clear() {
        hll_.clear();
    }

    size_t get_precision() const {
        return hll_.p();
    }

    std::pair<double, double> get_error_bounds() const {
        double est = hll_.report();
        double rel_err = hll_.est_err();
        return std::make_pair(est * (1.0 - rel_err), est * (1.0 + rel_err));
    }
};

PYBIND11_MODULE(core, m) {
    m.doc() = R"pbdoc(
        Python bindings for Dashing2
        ---------------------------

        This module provides Python bindings for the Dashing2 C++ library,
        implementing efficient MinHash sketching for biological sequences.

        Key Classes:
            RMinHash: A MinHash sketch implementation for sequences
        
        Key Functions:
            sketch_fastx: Create sketches from FASTA/FASTQ files
            compute_distances: Compute pairwise distances between sketches
    )pbdoc";
    
    // Bind RMinHash class
    py::class_<PyRMinHash>(m, "RMinHash")
        .def(py::init<size_t, uint64_t>(),
             py::arg("sketch_size") = 1024,
             py::arg("seed") = 0,
             R"pbdoc(
                Create a new RMinHash sketch.

                Args:
                    sketch_size: Size of the sketch (default: 1024)
                    seed: Random seed for hashing (default: 0)
             )pbdoc")
        .def("add_sequence", &PyRMinHash::add_sequence,
             py::arg("sequence"),
             py::arg("kmer_size") = 31,
             "Add a sequence to the sketch using k-mers")
        .def("jaccard_index", &PyRMinHash::jaccard_index,
             py::arg("other"),
             "Compute Jaccard index between this sketch and another")
        .def("get_size", &PyRMinHash::get_size,
             "Get the size of the sketch")
        .def("get_mins", &PyRMinHash::get_mins,
             "Get the sketch values as a numpy array")
        .def("merge", &PyRMinHash::merge,
             py::arg("other"),
             "Merge another sketch into this one")
        .def("clear", &PyRMinHash::clear,
             "Clear all values from the sketch");
    
    // Bind HLL class
    py::class_<PyHLL>(m, "HLL")
        .def(py::init<int>(),
             py::arg("precision") = 14,
             R"pbdoc(
                Create a new HyperLogLog sketch.

                Args:
                    precision: Precision parameter p (default: 14).
                              Controls the accuracy and memory usage.
                              Memory used is approximately 2^p bytes.
                              Error is approximately 1.04/sqrt(2^p).
             )pbdoc")
        .def("add_sequence", &PyHLL::add_sequence,
             py::arg("sequence"),
             py::arg("kmer_size") = 31,
             "Add a sequence to the HLL sketch using k-mers")
        .def("get_estimate", &PyHLL::get_estimate,
             "Get the estimated number of distinct elements")
        .def("merge", &PyHLL::merge,
             py::arg("other"),
             "Merge another HLL sketch into this one")
        .def("clear", &PyHLL::clear,
             "Clear all values from the sketch")
        .def("get_precision", &PyHLL::get_precision,
             "Get the precision parameter p used by this sketch")
        .def("get_error_bounds", &PyHLL::get_error_bounds,
             "Get the estimated lower and upper bounds (95% confidence interval)");
    
    // Bind standalone functions
    m.def("sketch_fastx", &sketch_fastx,
          "Sketch sequences from a FASTA/FASTQ file",
          py::arg("path"),
          py::arg("sketch_size") = 1024,
          py::arg("kmer_size") = 31,
          py::arg("canonicalize") = true,
          py::arg("window_size") = 0);
    
    m.def("compute_distances", &compute_distances,
          "Compute pairwise distances between sketches",
          py::arg("sketches"));
}
