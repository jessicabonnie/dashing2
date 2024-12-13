#include "fastxsketch.h"

namespace dashing2 {

/**
 * Merges multiple SketchingResult objects into a single result
 *
 * Takes an array of SketchingResult objects and combines their contents, including names,
 * signatures, kmers, cardinalities, and sequences if present. Optionally prepends provided
 * names to sequence identifiers.
 *
 * @param start Pointer to array of SketchingResult objects to merge
 * @param n Number of SketchingResult objects to merge
 * @param names Optional vector of names to prepend to sequence identifiers
 * @return Merged SketchingResult containing combined data from all inputs
 */
SketchingResult SketchingResult::merge(SketchingResult *start, size_t n, const std::vector<std::string> &names=std::vector<std::string>()) {
    DBG_ONLY(std::fprintf(stderr, "About to merge from %p of size %zu, names has size %zu\n", (void *)start, n, names.size());)
    SketchingResult ret;
    
    // Handle empty or single input cases
    if(n == 0) return ret;
    else if(n == 1) {
        ret = std::move(*start);
        std::transform(ret.names_.begin(), ret.names_.end(), ret.names_.begin(), [&names](const auto &x) {return names.front() + ":" + x;});
        return ret;
    }
    //ret.nperfile_.resize(total_seq);
    // Combine nperfile_ vectors from all inputs
    for(size_t i = 0; i < n; ++i) {
        ret.nperfile_.insert(ret.nperfile_.end(), start[i].nperfile_.begin(), start[i].nperfile_.end());
    }

    // Calculate total sizes and offsets
    size_t total_seqs = 0, total_sig_size = 0;
    std::vector<size_t> offsets(n + 1);
    std::vector<size_t> sig_offsets(n + 1);
    for(size_t i = 0; i < n; ++i) {
        const size_t nseqsi = start[i].names_.size();
        const size_t nregs = start[i].signatures_.size();
        assert(start[i].names_.size() == start[i].cardinalities_.size());
        total_seqs += nseqsi;
        total_sig_size += nregs;
        offsets[i + 1] = total_seqs;
        sig_offsets[i + 1] = total_sig_size;
    }

    // Resize result containers
    ret.names_.resize(total_seqs);
    if(std::any_of(start, start + n, [](auto &x) {return x.sequences_.size();})) {
        seq_resize(ret.sequences_, total_seqs);
    }
    const size_t sketchsz = start->signatures_.size() / start->names_.size();
    if(total_sig_size > 0) {
        ret.signatures_.resize(total_sig_size);
    }
    if(start->kmers_.size()) {
        ret.kmers_.resize(total_seqs * sketchsz);
    }
    ret.cardinalities_.resize(total_seqs);
    if(start->kmercounts_.size()) {
        ret.kmercounts_.resize(total_sig_size);
    }

    // Track what data types are present
    const bool seqsz = total_seqs,
               kmercountsz = !start->kmercounts_.empty();

    // Merge data from each input
    for(size_t i = 0; i < n; ++i) {
        auto &src = start[i];
        assert(src.names_.size() == offsets[i + 1] - offsets[i]);
        const auto ofs = offsets[i], sofs = sig_offsets[i];
        
        // Extract filename for prepending to sequence names
        std::string fname;
        if(names.size() > i) fname = names[i].substr(0, names[i].find_first_of(' '));
        // Append filename to sequence names to ensure seq names
        // Copy and transform names
        std::transform(src.names_.begin(), src.names_.end(), &ret.names_[ofs], [&fname](const auto &x) {
            return x + ':' + fname;
        });

        // Copy remaining data fields
        std::copy(src.cardinalities_.begin(), src.cardinalities_.end(), &ret.cardinalities_.at(ofs));
        if(seqsz) {
            ret.sequences_.add_set(src.sequences_.begin(), src.sequences_.end());
        }
        if(!start[i].signatures_.empty())
            std::copy(src.signatures_.begin(), src.signatures_.end(), &ret.signatures_.at(sofs));
        if(!start[i].kmers_.empty())
            std::copy(src.kmers_.begin(), src.kmers_.end(), &ret.kmers_[sofs]);
        if(kmercountsz)
            std::copy(src.kmercounts_.begin(), src.kmercounts_.end(), &ret.kmercounts_[sofs]);
    }
    return ret;
}

/**
 * Generates output destination path based on input parameters and options
 *
 * Constructs a filename incorporating various sketch parameters and options including:
 * - Seed value
 * - Canonicalization
 * - Spaced seed pattern
 * - Sketch size
 * - k-mer size
 * - Window size
 * - Count threshold
 * - Space type
 * - Compression parameters
 *
 * @param opts Dashing2Options containing sketch parameters
 * @param path Input path to base output name on
 * @param iskmer Whether input is k-mer based
 * @return Constructed output path incorporating relevant parameters
 */
std::string makedest(Dashing2Options &opts, const std::string &path, bool iskmer) {
    std::string ret(path);
    ret = ret.substr(0, ret.find_first_of(' '));
    
    // Handle path trimming and output prefix
    if(opts.trim_folder_paths() || opts.outprefix_.size()) {
        ret = trim_folder(path);
        if(opts.outprefix_.size())
            ret = opts.outprefix_ + '/' + ret;
    }

    // Add parameter suffixes
    if(opts.seedseed_ != 0)
        ret += ".seed" + std::to_string(opts.seedseed_);
    if(opts.canonicalize())
        ret += ".rc_canon";
    if(!opts.sp_.unspaced()) {
        ret += opts.sp_.to_string();
    }
    if(opts.kmer_result_ <= FULL_SETSKETCH)
        ret = ret + std::string(".sketchsize") + std::to_string(opts.sketchsize_);
    ret = ret + std::string(".k") + std::to_string(opts.k_);
    if(opts.w_ > opts.k_) {
        ret = ret + std::string(".w") + std::to_string(opts.w_);
    }
    if(opts.count_threshold_ > 0) {
        ret = ret + ".ct_threshold";
        if(std::fmod(opts.count_threshold_, 1.)) ret = ret + std::to_string(opts.count_threshold_);
        else ret = ret + std::to_string(int(opts.count_threshold_));
    }
    if(opts.sspace_ != SPACE_SET && opts.sspace_ != SPACE_EDIT_DISTANCE) {
        ret += '.';
        ret += to_string(opts.ct());
        if(opts.ct() != EXACT_COUNTING)
            ret += std::to_string(opts.cssize_);
    }
    if(opts.sspace_ == SPACE_SET && opts.sketch_compressed()) {
        char buf[256];
        auto l = std::sprintf(buf, ".a=%0.16Lg.b=%0.16Lg.fd=%0.16Lg", opts.compressed_a_, opts.compressed_b_, static_cast<long double>(opts.fd_level_));
        ret += std::string(buf, l);
    }
    ret += ".";
    if(opts.kmer_result_ <= FULL_SETSKETCH)
        ret += to_string(opts.sspace_);
    else {
        auto ks = opts.kmer_result_;
        if(iskmer && ks == FULL_MMER_COUNTDICT) {
            ks = FULL_MMER_SET;
        }
        ret += to_string(ks);
    }
    ret = ret + "." + bns::to_string(opts.rht_) + to_suffix(opts);
    DBG_ONLY(std::fprintf(stderr, "Source %s->%s\n", path.data(), ret.data());)
    return ret;
}


}
