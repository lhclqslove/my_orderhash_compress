//
// Created by USER on 2022/6/1.
//
#include <random>
#include "ReadFilter.h"
ReadFilter::~ReadFilter() {}
void OrderHashReadFilter::initialize(ReadData &rD) {
    this->rD = &rD;
    numReads = rD.getNumReads();
    std::vector<kMer_t> sketches(n * numReads);//每序列得到n个hashval,不满足条件的max_occ序列会hashval为0
    generateRandomNumbers(n);
    size_t maxNumkMers;
    if (rD.maxReadLen < k-1)
        maxNumkMers = 0;
    else
        maxNumkMers = rD.maxReadLen - k + 1;
#pragma omp parallel
    {
        // We define these vectors here rather than allocate inside string2Sketch
        // to avoid thread contention during repeated allocation and deallocation.
        // Note that memory allocation typically leads to waits when multiple threads
        // do it at the same time.
        std::vector<kMer_t> kMersVec(maxNumkMers), hashesVec(n);
        std::string readStr;
#pragma omp for
        for (read_t i = 0; i < numReads; ++i){
            rD.getRead(i, readStr);
            string2Sketch(readStr, sketches.data() + i * n, kMersVec, hashesVec);
        }
    } // pragma omp parallel

}
void OrderHashReadFilter::generateRandomNumbers(size_t n) {
    if (randNumbers)
        delete[] randNumbers;
    randNumbers = new kMer_t[n];

    std::random_device rd;
    std::mt19937_64 gen(rd());

    /* This is where you define the number generator for unsigned long long: */
    std::uniform_int_distribution<unsigned long long> dis;

    for (size_t i = 0; i < n; ++i) {
        randNumbers[i] = dis(gen);
    }
}
