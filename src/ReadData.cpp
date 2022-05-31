//
// Created by USER on 2022/5/31.
//
#include "ReadData.h"
#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <limits> // std::numeric_limits
#include <stdexcept>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
void ReadData::loadFromFile(const char *filename, enum Filetype filetype) {
    switch(filetype){
        case FASTQ:
            loadFromFastqFile(filename, false);
        case GZIP:
            loadFromFastqFile(filename, true);
    }
}
void ReadData::loadFromFastqFile(const char *fileName, bool gzip_flag) {
    numReads=0;
    readData.clear();
    std::ifstream infile;
    boost::iostreams::filtering_streambuf<boost::iostreams::input> *inbuf;
    std::istream *fin = &infile;
    if (gzip_flag) {
        infile.open(fileName, std::ios_base::binary);
        inbuf =
                new boost::iostreams::filtering_streambuf<boost::iostreams::input>;
        inbuf->push(boost::iostreams::gzip_decompressor());
        inbuf->push(infile);
        fin = new std::istream(inbuf);
    } else {
        infile.open(fileName);
    }
    std::string line;
    size_t totalNumBases = 0;
    maxReadLen = 0;
    size_t numReadsPerBlock = 5000;
    std::vector<std::string> lines(numReadsPerBlock);
    size_t numReadsCurrBlock = 0;
    size_t numReadsInserted = 0;
    while (true) {
        while (std::getline(*fin, line)) {
            std::getline(*fin, lines[numReadsCurrBlock++]);
            auto readLen = lines[numReadsCurrBlock-1].size();
            totalNumBases += readLen;
            if (readLen > maxReadLen)
                maxReadLen = readLen;
            numReads++;
            if (numReads == std::numeric_limits<read_t>::max())
                throw std::runtime_error(
                        "Too many reads for read_t type to handle.");
//            readPos.push_back(0);
            std::getline(*fin, line);
            std::getline(*fin, line);
            if (numReadsCurrBlock == numReadsPerBlock)
                break;
        }
        readData.resize(numReadsInserted + numReadsCurrBlock);
#pragma omp parallel for
        for (size_t i = 0; i < numReadsCurrBlock; i++) {
            std::unique_ptr<Read> ptr(new Read(
                    lines[i],numReadsInserted+i,que_cnt));
            readData[numReadsInserted+i] = std::move(ptr);
        }
        numReadsInserted += numReadsCurrBlock;
        if (numReadsCurrBlock < numReadsPerBlock)
            break;
        numReadsCurrBlock = 0;
    }
    assert(numReads != 0);
    avgReadLen = totalNumBases / numReads;
    std::cout << "numReads " << numReads << std::endl;
    std::cout << "avgReadLen " << avgReadLen << std::endl;
    std::cout << "maxReadLen " << maxReadLen << std::endl;
    if (gzip_flag) {
        delete fin;
        delete inbuf;
    }
    // close files
    infile.close();
}