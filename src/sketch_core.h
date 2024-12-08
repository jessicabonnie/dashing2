#ifndef SKETCH_CORE_H__
#define SKETCH_CORE_H__
#include "fastxsketch.h"
#include "lfsketch.h"
#include "bwsketch.h"
#include "bedsketch.h"
namespace dashing2 {

/**
 * Gets file sizes for a list of paths
 * @param paths Vector of file paths to get sizes for
 * @return Vector of pairs containing file index and size in bytes
 */
std::vector<std::pair<size_t, uint64_t>> get_filesizes(const std::vector<std::string> &paths);

/**
 * Core sketching function that processes input files and generates sketches
 * @param result SketchingResult object to store the sketches and metadata
 * @param opts Options controlling the sketching behavior
 * @param paths Vector of input file paths to sketch
 * @param outfile Output file path to write results
 * @return Reference to the populated SketchingResult
 */
SketchingResult &sketch_core(SketchingResult &, Dashing2DistOptions &opts, const std::vector<std::string> &paths, std::string &outfile);

}

#endif
