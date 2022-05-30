#include <iostream>
#include <csignal>
#include <chrono>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include "DirectoryUtils.h"
#include <omp.h>
using namespace std::chrono;
using namespace std;

std::string temp_dir_global;  // for interrupt handling
bool temp_dir_flag_global = false;

void signalHandler(int signum) {
    std::cout << "Interrupt signal (" << signum << ") received.\n";
    std::cout << "Program terminated unexpectedly\n";
    if (temp_dir_flag_global) {
        std::cout << "Deleting temporary directory:" << temp_dir_global << "\n";
        boost::filesystem::remove_all( temp_dir_global);
    }
    exit(signum);
}
int main(int argc, char **argv) {
    // register signal SIGINT and signal handler
    signal(SIGINT, signalHandler);
    auto start = high_resolution_clock::now();
    omp_set_nested(1);
    std::srand(unsigned(std::time(0)));
    //program options
    namespace po = boost::program_options;
    bool help_flag = false, compress_flag = false, decompress_flag = false;
    std::string infile, outfile;
    int num_thr;
    std::string working_dir;
    size_t k, n,l, m_k, m_w, max_chain_iter, edge_threshold;
    po::options_description desc("Allowed options");
    desc.add_options()("help,h", po::bool_switch(&help_flag), "produce help message")(
            "compress,c", po::bool_switch(&compress_flag), "compress")(
            "decompress,d", po::bool_switch(&decompress_flag), "decompress")(
            "input-file,i",
            po::value<std::string>(&infile),
            "input file name")(
            "output-file,o",po::value<std::string>(&outfile),
            "output file name")(
            "num-threads,t", po::value<int>(&num_thr)->default_value(20),
            "number of threads (default 20)")(
            "kmer,k", po::value<size_t>(&k)->default_value(23))(
            "num_kmer,l",po::value<size_t>(&l)->default_value(2),"orderhash number kmer for per hashval")(
            "num-hash,n", po::value<size_t>(&n)->default_value(60),
            "number of hash functions for minhash (default 60)")(
            "minimap-k", po::value<size_t>(&m_k)->default_value(20),
            "kmer size for the minimap2 (default 20)")(
            "minimap-w, w", po::value<size_t>(&m_w)->default_value(50),
            "window size for the minimap2 (default 50)")(
            "max-chain-iter", po::value<size_t>(&max_chain_iter)->default_value(400),
            "the max number of partial chains during chaining for minimap2 (default 400)")(
            "working-dir,w", po::value<std::string>(&working_dir)->default_value("."),
            "directory to create temporary files (default current directory)");
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    if (help_flag) {
        std::cout << desc << "\n";
        return 0;
    }
    if(infile.empty()){
        std::cout<<"No input file specified"<<"\n";
        std::cout << desc << "\n";
        return 0;
    }
    if(outfile.empty()){
        std::cout<<"No output file specified"<<"\n";
        std::cout << desc << "\n";
        return 0;
    }
    if ((!compress_flag && !decompress_flag) ||
        (compress_flag && decompress_flag)) {
        std::cout
                << "Exactly one of compress or decompress needs to be specified \n";
        std::cout << desc << "\n";
        return 1;
    }
    // generate randomly named temporary directory in the working directory
    std::string temp_dir;
    while (true) {
        std::string random_str = "tmp." + DirectoryUtils::random_string(10);
        temp_dir = working_dir + "/" + random_str + "/";
        if (!boost::filesystem::exists(temp_dir)) break;
    }
    if (!boost::filesystem::create_directory(temp_dir)) {
        throw std::runtime_error("Cannot create temporary directory.");
    }
    std::cout << "Temporary directory: " << temp_dir << "\n";
    temp_dir_global = temp_dir;
    temp_dir_flag_global = true;
    return 0;
}
