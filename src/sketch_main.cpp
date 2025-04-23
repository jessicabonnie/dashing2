#include "sketch_core.h"
#include "options.h"
#include "cmp_main.h"


// Define sketch options structure with shared options
#define SKETCH_OPTS \
static option_struct sketch_long_options[] = {\
    SHARED_OPTS\
    {"whole-line", no_argument, 0, OPTARG_WHOLE_LINE},\
};

namespace dashing2 {

/**
 * Print usage information for sketch subcommand
 * 
 * Displays help text explaining the sketch command usage and options
 * to stderr.
 */
void sketch_usage() {
    std::fprintf(stderr, "dashing2 sketch <opts> [fastas... (optional)]\n"
                         "We use only m-mers; if w <= k, however, this reduces to k-mers if the -w/--window-size is unspecified.\n"
                         SHARED_DOC_LINES
    );
}

/**
 * Main function for sketch subcommand
 * 
 * Handles command line argument parsing and orchestrates the sketching workflow:
 * - Parses and validates command line options
 * - Sets up sketching parameters and options
 * - Processes input files to generate sketches
 * - Optionally performs distance comparisons between sketches
 * 
 * @param argc Number of command line arguments
 * @param argv Array of command line argument strings
 * @return 0 on successful execution, 1 on error or if help is displayed
 */
int sketch_main(int argc, char **argv) {
    // Initialize parameters with default values
    int c;
    int k = -1;                     // k-mer size
    int w = -1;                     // window size  
    int nt = -1;                    // number of threads
    SketchSpace sketch_space = SPACE_SET;  // sketching space type
    KmerSketchResultType res = ONE_PERM;   // sketching result type
    
    // Flags for various options
    bool save_kmers = false;        // save k-mers
    bool save_kmercounts = false;   // save k-mer counts
    bool cache = false;             // cache sketches
    bool use128 = false;            // use 128-bit hashing
    bool canon = true;              // use canonical k-mers
    bool exact_kmer_dist = false;   // compute exact k-mer distances
    bool hpcompress = false;        // homopolymer compression
    bool refine_exact = false;      // refine distances exactly
    bool fasta_dedup = false;       // deduplicate FASTA entries
    bool whole_line_sketch = false;  // enable whole line sketch
    
    // Parameters for compressed sketching
    long double compressed_a = -1.L;
    long double compressed_b = -1.L;
    
    // Distance calculation parameters
    double similarity_threshold = -1.;
    unsigned int count_threshold = 0.;
    size_t cssize = 0;             // compressed sketch size
    size_t sketchsize = 1024;      // sketch size
    
    // Input/output parameters
    std::string ffile;             // file containing paths
    std::string outfile;           // output file
    std::string qfile;             // query file
    std::string outprefix;         // output prefix
    std::string cmpout;            // comparison output file
    
    // Other parameters
    int option_index = 0;
    bns::RollingHashingType rht = bns::DNA;
    DataType dt = DataType::FASTX;
    OutputKind ok = SYMMETRIC_ALL_PAIRS;
    int topk_threshold = -1;
    int truncate_mode = 0;
    double nbytes_for_fastdists = sizeof(RegT);
    bool parse_by_seq = false;
    int by_chrom = false;
    double downsample_frac = 1.;
    uint64_t seedseed = 0;
    size_t batch_size = 0;
    int nLSH = 2;
    
    // Vector for storing pairs of genome IDs to compare
    std::vector<std::pair<uint32_t, uint32_t>> compareids;
    
    Measure measure = SIMILARITY;
    std::ios_base::sync_with_stdio(false);
    std::string fsarg;
    // By default, use full hash values, but allow people to enable smaller
    bool normalize_bed = false;
    OutputFormat of = OutputFormat::HUMAN_READABLE;
    std::string spacing;
    std::vector<std::string> paths;

    // Validate command line options
    validate_options(argv);
    
    // Define sketch options
    SKETCH_OPTS
    
    // Parse command line options
    for(;(c = getopt_long(argc, argv, "m:p:k:w:c:f:S:F:Q:o:L:CNs2BPWh?ZJGHv", sketch_long_options, &option_index)) >= 0;) {
        switch(c) {
            SHARED_FIELDS
            case OPTARG_WHOLE_LINE: whole_line_sketch = true; break;
            case OPTARG_HELP: case '?': case 'h': sketch_usage(); return 1;
        }
        //std::fprintf(stderr, "After getopt argument %d, of is %s\n",c , to_string(of).data());
    }

    // Set default k if not specified
    if(k < 0) k = nregperitem(rht, use128);
    
    // Set number of threads from environment if not specified
    if(nt < 0) {
        char *s = std::getenv("OMP_NUM_THREADS");
        if(s) nt = std::max(std::atoi(s), 1);
    }
    OMP_ONLY(omp_set_num_threads(nt));

    // Handle input paths
    if(compareids.empty()) {
        paths.insert(paths.end(), argv + optind, argv + argc);
    } else if(optind != argc) throw std::runtime_error("CLI paths must be empty to use pairlist mode.");
    std::unique_ptr<std::vector<std::string>> qup;
    if(ffile.size()) {
        if(!bns::isfile(ffile)) THROW_EXCEPTION(std::runtime_error("No path found at "s + ffile));
        std::ifstream ifs(ffile);
        static constexpr size_t bufsize = 1<<18;
        std::unique_ptr<char []> buf(new char[bufsize]);
        ifs.rdbuf()->pubsetbuf(buf.get(), bufsize);
        for(std::string l;std::getline(ifs, l);) {
            paths.push_back(l);
        }
        if(paths.empty()) {
            THROW_EXCEPTION(std::runtime_error("No paths read from "s + ffile));
        }
    }

    // Track number of reference sequences
    size_t nref = paths.size();
    
    // Read query paths if specified
    if(qfile.size()) {
        std::ifstream ifs(qfile);
        for(std::string l;std::getline(ifs, l);)
            paths.push_back(l);
    }
    size_t nq = paths.size() - nref;

    // Configure sketching options
    Dashing2Options opts(k, w, rht, sketch_space, dt, nt, use128, spacing, canon, res);
    opts
        .cache_sketches(cache)
        .cssize(cssize)
        .sketchsize(sketchsize)
        .save_kmers(save_kmers)
        .outprefix(outprefix)
        .save_kmercounts(save_kmercounts)
        .save_kmers(save_kmers)
        .parse_by_seq(parse_by_seq)
        .count_threshold(count_threshold)
        .homopolymer_compress_minimizers(hpcompress)
        .seedseed(seedseed)
        .fasta_dedup(fasta_dedup)
        .whole_line_sketch(whole_line_sketch);
    opts.by_chrom_ = by_chrom;
    opts.downsample(downsample_frac);
    opts.compressed_a_ = compressed_a;
    opts.compressed_b_ = compressed_b;
    opts.fd_level_ = nbytes_for_fastdists;
    opts.set_sketch_compressed();

    // Validate homopolymer compression setting
    if(hpcompress) {
        if(!opts.homopolymer_compress_minimizers_) THROW_EXCEPTION(std::runtime_error("Failed to hpcompress minimizers"));
    }

    // Set filterset if specified
    opts.filterset(fsarg);

    // Adjust sketch result type based on space
    if((opts.sspace_ == SPACE_PSET || opts.sspace_ == SPACE_MULTISET || opts.sspace_ == SPACE_EDIT_DISTANCE)
            && opts.kmer_result_ == ONE_PERM) {
        opts.kmer_result_ = FULL_SETSKETCH;
    }

    opts.bed_parse_normalize_intervals_ = normalize_bed;

    // Configure distance calculation options
    Dashing2DistOptions distopts(opts, ok, of, nbytes_for_fastdists, truncate_mode, topk_threshold, similarity_threshold, cmpout, exact_kmer_dist, refine_exact, nLSH);

    // Validate input paths
    if(paths.empty()) {
        std::fprintf(stderr, "No paths provided. See usage.\n");
        sketch_usage();
        return 1;
    }

    // Perform sketching
    SketchingResult result;
    if(verbosity >= EXTREME) {
        std::fprintf(stderr, "About to sketch\n");
    }
    sketch_core(result, distopts, paths, outfile);
    if(verbosity >= EXTREME) {
        std::fprintf(stderr, "Finished sketching\n");
    }

    // Set number of queries
    result.nqueries(nq);

    // Perform distance comparisons if requested
    if(cmpout.size()) {
        distopts.measure_ = measure;
        distopts.cmp_batch_size_ = default_batchsize(batch_size, distopts);
        cmp_core(distopts, result);
    }

    return 0;
}

} // namespace dashing2
