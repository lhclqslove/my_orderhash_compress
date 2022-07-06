//
// Created by USER on 2022/5/31.
//
#include "Compressor.h"
#include <boost/filesystem.hpp>
#include <chrono>
#include <fstream>
#include <iostream>
#include <ios>
#include <unistd.h>
#include <string>
#include <omp.h>
#include "bsc_helper.h"
#include "lzma2_helper.h"
#include "DirectoryUtils.h"
#include "ReadFilter.h"
#include "Consensus.h"
#ifdef USE_MALLOC_TRIM
    #include <malloc.h>
#endif
void mem_usage(double& vm_usage, double& resident_set) {
    // from https://www.tutorialspoint.com/how-to-get-memory-usage-at-runtime-using-cplusplus
    using namespace std;
    vm_usage = 0.0;
    resident_set = 0.0;
    ifstream stat_stream("/proc/self/stat",ios_base::in); //get info from proc
    // directory
    //create some variables to get info
    string pid, comm, state, ppid, pgrp, session, tty_nr;
    string tpgid, flags, minflt, cminflt, majflt, cmajflt;
    string utime, stime, cutime, cstime, priority, nice;
    string O, itrealvalue, starttime;
    unsigned long vsize;
    long rss;
    stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
                >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
                >> utime >> stime >> cutime >> cstime >> priority >> nice
                >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care
    // about the rest
    stat_stream.close();
    long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // for x86-64 is configured
    // to use 2MB pages
    vm_usage = vsize / 1024.0;
    resident_set = rss * page_size_kb;
    cout << "resident_set: " << resident_set << " KB\n";
}
void Compressor::compress(const char *inputFileName, const int numThr) const {
    std::cout << "number of threads: " << numThr << std::endl;
    double vm_usage, resident_set;
    std::cout << "Entering compress():\n";
    size_t loop_index=0;
    mem_usage(vm_usage, resident_set);
    read_t mergeCnt=0;
    omp_set_num_threads(numThr);
    {
        //先读入文件，加载到数据类中
        ReadData rD;
        rD.tempDir = tempDir;
        ReadData::recover_cnt=recover_cnt;
        auto read_start = std::chrono::high_resolution_clock::now();
        rD.loadFromFile(inputFileName, filetype);

        auto read_end = std::chrono::high_resolution_clock::now();
        auto duration =
                std::chrono::duration_cast<std::chrono::milliseconds>(read_end - read_start);
        std::cout << "Reading file took " << duration.count()
                  << " milliseconds" << std::endl;
        std::cout << "After rD.loadFromFile():\n";
        mem_usage(vm_usage, resident_set);
        size_t lcnt=0;
//        size_t occ_cnt=0;
        while(rD.getNumReads()>ReadData::sequence_number_threshold&&loop_index<1)
//       while(rD.getNumReads()>1)
        {
                OrderHashReadFilter rF;
                rF.k = k;
                rF.n = n;
//                if(l-lcnt>1&&loop_index>0&&loop_index%4==0)lcnt++;
//                rF.l = l-lcnt;
                rF.l = l;
                rF.que_cnt=que_cnt;
//                if(loop_index>0&&loop_index%5==0)occ_cnt+=5;
                rF.max_occ = max_occ;
                rF.tempDir = tempDir;
                {
                    auto start = std::chrono::high_resolution_clock::now();
                    rF.initialize(rD);
                    auto end = std::chrono::high_resolution_clock::now();
                    auto duration =
                            std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
                    std::cout<<"loop_index"<<loop_index << "Calculating orderHashes took " << duration.count()
                              << " millis"
                                 ""
                                 "econds" << std::endl;
                }
                std::cout <<"loop_index"<<loop_index<< "After orderhash computation:\n";
                mem_usage(vm_usage, resident_set);
//                malloc_trim(0);
    #ifdef USE_MALLOC_TRIM
                malloc_trim(0); // clear memory
    #endif

                Consensus consensus;
                {
                    consensus.rD = &rD;
                    consensus.rF = &rF;
                    consensus.k= k;
                    consensus.tempDir = tempDir;
                    consensus.tempFileName = tempFileName;
                    consensus.numThr = numThr;
                    //parameters for minimap2
                    consensus.m_k = m_k;
                    consensus.m_w = m_w;
                    consensus.max_chain_iter = max_chain_iter;
                    consensus.edge_threshold = edge_threshold;
#ifdef  LOG
                    std::cout <<"loop_index"<<loop_index<<std::endl;
#endif
                    consensus.generateAndWriteConsensus(loop_index++,mergeCnt);
                }
        }
        /*
         * 把剩下下的序列全部写入lone文件中
         *
         * */
        Consensus consensus;
        {
            consensus.rD = &rD;
            consensus.tempDir = tempDir;
            consensus.tempFileName = tempFileName;
            consensus.writeLoneReadtofile(loop_index);
        }
        /*
         * 写入说明信息
         * */
        std::ofstream metaData;
        metaData.open(tempDir + "metaData");
        metaData<<"loop_index="<<loop_index<<'\n';
        std::cout << "numreads = " << ReadData::index<< '\n';
        metaData<< "numreads = " << ReadData::index<< '\n';
        metaData << "numThr=" << numThr << '\n';
        metaData.close();

    }
#ifdef  LOG
    std::cout <<"mergeCnt"<<mergeCnt<<std::endl;
#endif
#ifdef USE_MALLOC_TRIM
    malloc_trim(0); // clear memory
#endif

    {
// We first compress all the files in the temp directory (we compress
// only files with extensions)
        size_t total_compressed_size = 0;
        size_t total_uncompressed_size=0;
        std::cout << "bsc/lzma2 compression starts" << std::endl;
        std::set<std::string> extensions;
        DirectoryUtils::getAllExtensions(tempDir, std::inserter(extensions, extensions.end()));
        for (const std::string &ext: extensions) {
            std::vector<size_t> uncompressedSizes(numThr);
            std::vector<size_t> compressedSizes(numThr);
#pragma omp parallel num_threads(numThr)
#pragma omp for
            for (int i = 0; i < numThr; i++) {
                for(int j=0;j<=loop_index;j++)
                {
                    if(j==loop_index&&i!=0)continue;
                    std::string fullPath = tempDir + tempFileName +".loop."+ std::to_string(j) + ".tid." + std::to_string(i) + ext;
                    std::string outPath = fullPath + "Compressed";
                    // use lzma2 for .base to get better compression
//                    if(ext==std::string(".loneid"))continue;
                    if (ext == std::string(".base"))
                        lzma2::lzma2_compress(fullPath.c_str(), outPath.c_str());
                    else
                        bsc::BSC_compress(fullPath.c_str(), outPath.c_str());
                    uncompressedSizes[i] += boost::filesystem::file_size(fullPath);
//



                    compressedSizes[i] += boost::filesystem::file_size(outPath);

                    boost::filesystem::remove(fullPath);
                }

            }
            size_t totalUncompressed = 0, totalCompressed = 0;
            for (int i = 0; i < numThr; i++) {
                totalUncompressed += uncompressedSizes[i];
                totalCompressed += compressedSizes[i];
            }
            total_uncompressed_size+=totalUncompressed;
            total_compressed_size += totalCompressed;
            std::cout << "Extension " << ext << ": Compressed " << totalUncompressed << " bytes to " << totalCompressed
                      << " bytes\n";
        }
        std::cout <<"Total uncompressed size of all streams is"<<total_uncompressed_size<<" bytes\n"<< "Total compressed size of all streams is " << total_compressed_size << " bytes\n";
    }

    std::cout << "Creating tar archive ..." << std::endl;
    std::string tar_command =
            "tar -cf " + outputFileName + " -C " + tempDir + " . ";
    int tar_status = std::system(tar_command.c_str());
    if (tar_status)
        throw std::runtime_error(
                "Error occurred during tar archive generation.");
    std::cout << "Tar archive done!\n";

    {
        std::string lsCommand = "ls -l " + outputFileName;
        int ls_status = std::system(lsCommand.c_str());
        if (ls_status)
            throw std::runtime_error("Error occurred during ls command.");
    }

}