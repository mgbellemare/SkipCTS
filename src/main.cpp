/******************************
      Author: Joel Veness
        Date: 2011-2013
******************************/

#include <iostream>
#include <fstream>
#include <bitset>
#include <cmath>
#include <sstream>

#include "common.hpp"
#include "ctw.hpp"
#include "facmodels.hpp"
#include "cts.hpp"
#include "ac.hpp"
#include "skipcts.hpp"

// boost includes
#include <boost/program_options.hpp>
#include <boost/progress.hpp>
#include <boost/filesystem/operations.hpp>


namespace po = boost::program_options;

static po::variables_map options;
static po::options_description options_desc("Usage", 100);


// bit history data structure
history_t history;


/* display the program usage options and quit. */
static void showHelp(void) {

    std::cout << std::endl << options_desc << std::endl;
    exit(0);
}


/* load program options into a variable map */
static void initOptions(int argc, char *argv[], po::variables_map &vm) {

    std::string file_msg;

    file_msg += "File to compress/decompress. A compressed file with the\n"
                "same name plus an extension (e.g. foo.cts) will be\n"
                "produced. Decompression is chosen automatically if\n"
                "the file extension matches a method (i.e. .ctw for\n"
                "ctw.)";

    options_desc.add_options()

        ("help", "")

        ("depth",
         po::value<unsigned int>()->default_value(48),
         "Maximum depth to look back in bits. Higher values use\nmore RAM.")

        ("skips",
         po::value<unsigned int>()->default_value(0),
         "Number of context skips, if a skipping model is used.")

        ("slots",
         po::value<unsigned int>()->default_value(18),
         "The base 2 logarithm of the number of slots in the \n"
          "hashtable (SkipCTS only).")

        ("method",
         po::value<std::string>()->default_value("faccts"),
         "Compression method. (ctw, cts, facctw, faccts, skipcts)")

        ("file",
         po::value<std::string>(),
         file_msg.c_str())
    ;

    po::store(po::parse_command_line(argc, argv, options_desc), vm);
    po::notify(vm);

    if (argc == 1) showHelp();
}


/* checks for illegal combinations of configuration parameters,
   throwing an exception if an illegal combination is detected. */
static void processOptions(const po::variables_map &vm) {

    if (vm.count("help") > 0) showHelp();
}


/* We initialize the context with zeros for convenience. */ 
void zeroFill(history_t &h, size_t n) {

    for (size_t i=0; i < n; ++i) h.push_back(0);
}


/* create the file compressor, or die if creation fails */
static Compressor *buildCompressor() {

    unsigned int d     = options["depth"].as<unsigned int>();
    unsigned int skips = options["skips"].as<unsigned int>();
    unsigned int slots = options["slots"].as<unsigned int>();

    if (slots < 10) throw std::invalid_argument("slots must be at least 10.");

    // ensure history always has enough context
    history.clear();
    zeroFill(history, d);

    Compressor *c = NULL;
    if (options["method"].as<std::string>() == "cts") {
        c = new SwitchingTree(history, d);
    } else if (options["method"].as<std::string>() == "ctw") {
        c = new ContextTree(history, d);
    } else if (options["method"].as<std::string>() == "facctw") {
        c = new FactoredContextTree(history, d);
    } else if (options["method"].as<std::string>() == "faccts") {
        c = new FactoredSwitchingTree(history, d);
    } else if (options["method"].as<std::string>() == "skipcts") {
        c = new SkipCTS(history, d, skips, slots);
    } else if (options["method"].as<std::string>() == "facskipcts") {
        c = new FactoredSkipCTS(history, d);
    } 

    if (c == NULL) {
        std::cout << "Unknown compression method." << std::endl;
        std::exit(1);
    }

    return c;
}


/* compress the input file */
static void compress() {

    std::auto_ptr<Compressor> c(buildCompressor());
    std::ifstream::pos_type size;

    // open the file to be compressed
    std::string fn = options["file"].as<std::string>();
    std::ifstream file(fn.c_str(), std::ios::in | std::ios::binary | std::ios::ate);
    if (!file.is_open()) {
        throw std::runtime_error("cannot load or open input file");
    }
    size = file.tellg(); file.seekg(0, std::ios::beg);
    std::cerr << "Compressing " << fn << " (" << size << " bytes)" << std::endl;

    // determine the destination filename
    std::string ext = std::string(".") + c->fileExtension();
    std::string ofn = fn + ext;

    // now compress the file, bit by bit
    BinaryArithmeticEncoder bae(ofn);
    std::ostringstream oss; boost::progress_timer ti(oss);
    boost::progress_display show_progress(static_cast<unsigned long>(size));

    for (int i=0; i < size; ++i) {
        
        char buffer = file.get();

        for (size_t j=0; j < 8; j++) {
        
            int bit = buffer & 1;
            double p = c->prob(1);
            double before = c->logBlockProbability();
            bae.encode(bit, p);
            c->update(bit);
            double after = c->logBlockProbability();
           
            assert(after < before);
            assert(p > 0.0 && p < 1.0);

            buffer >>= 1;
        }
        ++show_progress;
    }
    bae.close();
    file.close();
    std::cerr << std::endl;

    // display compression statistics
    boost::uintmax_t compr_size = boost::filesystem::file_size(ofn) * 8;
    double bpc = double(compr_size) / double(size);
    double ratio = double(compr_size) / double(size*8) * 100.0;
    double kbsec = double(size / 1024) / ti.elapsed();
    std::cout << "Compression statistics: " << std::endl;
    std::cout << "    Ratio:          " <<  ratio << "%" << std::endl;
    std::cout << "    Size:           " <<  (compr_size/8) << " bytes" << std::endl;
    std::cout << "    Bits per Byte:  " <<  bpc << std::endl;
    std::cout << "    Time:           " << ti.elapsed() << "s" << std::endl;
    std::cout << "    KB/s:           " << kbsec << std::endl;
}


/* decompress the input file */
static void decompress() {

    std::string fn  = options["file"].as<std::string>();
    BinaryArithmeticDecoder dec(fn.c_str());
    std::auto_ptr<Compressor> dc(buildCompressor());

    std::cerr << "Decompressing " << fn << " (" << (dec.size()/8) << " bytes)" << std::endl;

    // open destination filename
    std::string decfn = fn.substr(0, fn.find_last_of("."));
    std::ofstream decoded(decfn.c_str(), std::ios::out | std::ios::binary);
    if (!decoded.is_open()) {
        throw std::runtime_error("cannot open/write to ");
    }

    // now decode the file, one bit at a time
    char buf = 0, upto = 0;
    std::ostringstream oss; boost::progress_timer ti(oss);
    boost::progress_display dec_progress(dec.size());
    for (uint32_t i=0; i < dec.size(); ++i) {

        double p = dc->prob(1);

        bit_t bit = dec.decode(p);
        dc->update(bit);

        buf |= char(bit) << upto;
        upto++;
        if (upto == 8) {
            decoded.put(buf);
            buf = 0; upto = 0;
        }

        ++dec_progress;
    }

    dec.close();
    decoded.close();

    // display decompression statistics
    double kbsec = double(dec.size() / 8 / 1024) / ti.elapsed();

    std::cerr << std::endl;
    std::cout << "Decompression statistics: " << std::endl;
    std::cout << "    Time:           " << ti.elapsed() << "s" << std::endl;
    std::cout << "    KB/s:           " << kbsec << std::endl;
}


/* determine whether we are compressing or decompressing */
static bool compressionMode() {

    // if file ends in ".cts" or ".ctw" then we are in decompression mode
    std::string fn  = options["file"].as<std::string>();
    size_t pos = fn.find_last_of(".");
    if (pos == std::string::npos) return true;

    std::string ext = fn.substr(pos);
    if (ext == ".ctw"    || ext == ".cts" || 
        ext == ".faccts" || ext == ".facctw" || 
        ext == ".skipcts" || ext == ".facskipcts") return false;

    return true;
}


/* application entry point */
int main(int argc, char* argv[]) {

    // load program configuration
    initOptions(argc, argv, options);

    // check and apply program configuration
    processOptions(options);

    // determine whether we are compressing or decompressing
    try {
        if (compressionMode()) {
            compress();
        } else {
            decompress();
        }
    } catch (std::exception &e) {
        std::cout << "error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}



