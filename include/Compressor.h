//
// Created by USER on 2022/5/31.
//

#ifndef MY_ORDERHASH_COMPRESS_COMPRESSOR_H
#define MY_ORDERHASH_COMPRESS_COMPRESSOR_H

#include <iostream>
#include "ReadData.h"
#include "ThreadPool.h"

class Compressor {
public:
    /** Parameters for filtering **/
    size_t k, n,l, m_k, max_occ,m_w, max_chain_iter, edge_threshold,que_cnt,recover_cnt;
    ReadData::Filetype filetype = ReadData::Filetype::READ;
//    ReadAligner *rA;
//    size_t loop_index;
    std::string tempDir;
    /** The temp directories **/
    // std::string tempDir = "tempRaw/";
    /** The output filenames to use in temp directories **/
    std::string tempFileName = "Stream";

    std::string outputFileName = "compressedFile";

    bool low_mem;
    // std::string tarFileName = "originalFile";

    /**
     * @brief Compresses the read data file inputFileName and stores the result
     * in this->outputFileName
     *
     * @param inputFileName
     */
    void compress(const char *inputFileName, const int numThr) const;

};
#endif //MY_ORDERHASH_COMPRESS_COMPRESSOR_H
