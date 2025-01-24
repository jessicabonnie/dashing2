#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "sketch/sketch.h"  // Include main dashing2 header
#include "sketch/dist.h"
#include "fastxsketch.h"
#include "counter.h"

namespace py = pybind11;

// Wrapper for FinalRMinHash class
class PyRMinHash {
    sketch::mh::FinalRMinHash sketch_;

public:
    PyRMinHash(size_t sketch_size = 1024, uint64_t seed = 0) 
        : sketch_(sketch_size, seed) {}

    void add_sequence(const std::string& seq, size_t kmer_size = 31) {
        sketch::Encoder<> enc(kmer_size);
        sketch_.add_string(seq.data(), seq.size(), enc);
    }

    double jaccard_index(const PyRMinHash& other) const {
        return sketch_.jaccard_index(other.sketch_);
    }

    size_t get_size() const {
        return sketch_.size();
    }

    // Add method to get the sketch values as numpy array
    py::array_t<uint64_t> get_mins() const {
        auto result = py::array_t<uint64_t>(sketch_.size());
        auto r = result.mutable_unchecked<1>();
        for(size_t i = 0; i < sketch_.size(); ++i) {
            r(i) = sketch_.at(i);
        }
        return result;
    }
};

// Wrapper for FastxSketchingResult
py::dict sketch_fastx(const std::string& path, 
                     size_t sketch_size = 1024,
                     size_t kmer_size = 31,
                     bool canonicalize = true,
                     size_t window_size = 0) {
    
    sketch::FastxSketchingResult result = sketch::fastx_sketch(path, 
        sketch::FastxSketchingOptions()
            .sketch_size(sketch_size)
            .kmer_length(kmer_size)
            .canonicalize(canonicalize)
            .window_size(window_size));

    // Convert result to Python dictionary
    py::dict ret;
    ret["names"] = py::cast(result.names_);
    
    // Convert sketches to numpy array
    auto sketches = py::array_t<uint64_t>({result.sketches_.size(), sketch_size});
    auto r = sketches.mutable_unchecked<2>();
    for(size_t i = 0; i < result.sketches_.size(); ++i) {
        for(size_t j = 0; j < sketch_size; ++j) {
            r(i, j) = result.sketches_[i].at(j);
        }
    }
    ret["sketches"] = sketches;
    
    return ret;
}

// Compute pairwise distances between sketches
py::array_t<double> compute_distances(py::array_t<uint64_t> sketches) {
    auto buf = sketches.request();
    if (buf.ndim != 2) 
        throw std::runtime_error("Number of dimensions must be 2");
    
    size_t n_sketches = buf.shape[0];
    size_t sketch_size = buf.shape[1];
    
    // Create result matrix
    auto distances = py::array_t<double>({n_sketches, n_sketches});
    auto dist = distances.mutable_unchecked<2>();
    
    // Get pointer to sketch data
    uint64_t* ptr = static_cast<uint64_t*>(buf.ptr);
    
    // Compute pairwise distances
    for(size_t i = 0; i < n_sketches; ++i) {
        dist(i, i) = 0.0;  // Distance to self is 0
        for(size_t j = i + 1; j < n_sketches; ++j) {
            // Create temporary sketches
            sketch::mh::FinalRMinHash s1(sketch_size), s2(sketch_size);
            std::copy(ptr + i * sketch_size, 
                     ptr + (i + 1) * sketch_size, 
                     s1.begin());
            std::copy(ptr + j * sketch_size, 
                     ptr + (j + 1) * sketch_size, 
                     s2.begin());
            
            // Compute Jaccard distance
            double jac = s1.jaccard_index(s2);
            double dist_val = -std::log(jac);  // Convert to distance
            
            // Fill both triangles of the symmetric matrix
            dist(i, j) = dist_val;
            dist(j, i) = dist_val;
        }
    }
    
    return distances;
}

PYBIND11_MODULE(core, m) {
    m.doc() = "Python bindings for Dashing2"; // optional module docstring
    
    // Bind RMinHash class
    py::class_<PyRMinHash>(m, "RMinHash")
        .def(py::init<size_t, uint64_t>(),
             py::arg("sketch_size") = 1024,
             py::arg("seed") = 0)
        .def("add_sequence", &PyRMinHash::add_sequence,
             py::arg("sequence"),
             py::arg("kmer_size") = 31)
        .def("jaccard_index", &PyRMinHash::jaccard_index)
        .def("get_size", &PyRMinHash::get_size)
        .def("get_mins", &PyRMinHash::get_mins);
    
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
