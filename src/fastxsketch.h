#pragma once
#ifndef DASHING2_FASTX_SKETCH_H__
#define DASHING2_FASTX_SKETCH_H__
#include "d2.h"
#include "mmvec.h"
#include "tmpseqs.h"
#include "cmp_main.h"
#include <variant>

namespace dashing2 {
using std::to_string;

/**
 * Helper function to convert pointer to string representation
 * @param ptr Pointer to convert
 * @return String representation of pointer address
 */
template<typename T>
static inline std::string to_string(const T *ptr) {
    std::ostringstream oss;
    oss << static_cast<const void *>(ptr);
    return oss.str();
}

// Flag indicating if sequences should be kept in memory
static bool seqs_in_memory = false;

struct Dashing2DistOptions;

/**
 * Class to store results from sketching sequences
 * Contains sketches, k-mers, cardinalities and other metadata
 */
struct SketchingResult {
    // Constructor initializes sequences with memory flag
    SketchingResult(): sequences_(seqs_in_memory) {}
    
    // Move constructor and assignment
    SketchingResult(SketchingResult &&o) = default;
    SketchingResult &operator=(SketchingResult &&o) = default;
    
    // Delete copy constructor and assignment
    SketchingResult(const SketchingResult &o) = delete;
    SketchingResult(SketchingResult &o) = delete;
    SketchingResult &operator=(const SketchingResult &o) = delete;

    std::vector<std::string> names_;              // List of files, potentially multiple per line
    std::vector<std::string> destination_files_;  // Contains sketches/kmer-sets,kmer-sequences etc.
    std::vector<std::string> kmerfiles_;          // Contains file-paths for k-mers, if saved
    std::vector<std::string> kmercountfiles_;     // Contains k-mer counts, if saved
    
    // nperfile_ is either empty (each filename/row has its own sketch)
    // or contains number of sketches per file/line
    std::vector<uint32_t> nperfile_;             
    
    std::vector<double> cardinalities_;           // Cardinality estimates
    tmpseq::MemoryOrRAMSequences sequences_;      // Sequence storage
    
    // Signatures matrix - only used for LSH pre-filtering with edit distance
    mm::vector<RegT> signatures_;    
    
    // TODO: mmap these matrices to reduce peak memory footprint
    mm::vector<uint64_t> kmers_;                 // K-mer storage
    std::vector<float> kmercounts_;              // K-mer count storage
    // This contains the k-mers corresponding to signatures, if asked for 128-bit k-mers, these are stored in chunks of 2 64-bit integers.
    size_t nq = 0;                               // Number of queries

    /**
     * Get total number of sequences
     * @return Sum of nperfile if nonempty, otherwise number of names
     */
    size_t total_seqs() const {
        // Sum of nperfile if nonempty
        // otherwise, just one sequence/bag of k-mers per "name"
        return nperfile_.size() ? std::accumulate(nperfile_.begin(), nperfile_.end(), size_t(0)): names_.size();
    }

    std::string str() const;
    static SketchingResult merge(SketchingResult *start, size_t n, const std::vector<std::string> &);
    void print();
    
    /**
     * Get number of queries
     * @return Number of queries
     */
    size_t nqueries() const {return nq;}

    /**
     * Set number of queries
     * @param nqnew New number of queries
     */
    void nqueries(size_t nqnew) {nq = nqnew;}
};

using FastxSketchingResult = SketchingResult;

// Function declarations for sequence resizing
void seq_resize(std::vector<std::string>& seqs, const size_t num_seqs);
void seq_resize(tmpseq::Seqs&seqs, const size_t num_seqs) noexcept;
void seq_resize(tmpseq::MemoryOrRAMSequences& , const size_t num_seqs);
int32_t num_threads();

// Core sketching function declarations
FastxSketchingResult &fastx2sketch(FastxSketchingResult &res, Dashing2Options &opts, const std::vector<std::string> &paths, std::string path);
FastxSketchingResult &fastx2sketch_byseq(FastxSketchingResult &res, Dashing2DistOptions &opts, const std::string &path, kseq_t *kseqs, std::string outpath, bool parallel=false, const size_t seqs_per_batch = 8192);
std::string makedest(Dashing2Options &opts, const std::string &path, bool iskmer=false);

namespace variation {
// Type definitions for different set sketch variants
using ByteSetS = sketch::setsketch::CFByteSetS;
using NibbleSetS = sketch::setsketch::CFNibbleSetS;
using ShortSetS = sketch::setsketch::CFShortSetS;
using UintSetS = sketch::setsketch::CFUintSetS;
using VSetSketch = std::variant<NibbleSetS, ByteSetS, ShortSetS, UintSetS>;

/**
 * Get raw data pointer from variant set sketch
 * @param o Variant set sketch
 * @return Pointer to underlying data
 */
INLINE const RegT *getdata(VSetSketch &o) {
    const RegT *ret;
    std::visit([&ret](auto &x) {ret = (const RegT *)x.data();}, o);
    return ret;
}

/**
 * Get cardinality from variant set sketch
 * @param o Variant set sketch
 * @return Cardinality value
 */
INLINE double getcard(VSetSketch &o) {
    double ret;
    std::visit([&ret](auto &x) {ret = x.getcard();}, o);
    return ret;
}
} // namespace variation
using variation::VSetSketch;
} // namespace dashing2

#endif
