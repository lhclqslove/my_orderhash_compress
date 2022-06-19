//
// Created by USER on 2022/6/7.
//

#ifndef MY_ORDERHASH_COMPRESS_CONSENSUS_H
#define MY_ORDERHASH_COMPRESS_CONSENSUS_H
#include "ConsensusGraph.h"
#include "Edits.h"
#include "OmpMutex.h"
#include "ReadData.h"
#include "omp.h"
#include "ReadFilter.h"
#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>


// structure for stats
struct CountStats {
    size_t countMinHash; // number of reads passing minhash filter
    size_t countMinHashNotInGraph; // number of reads passing minhash filter that are not already in graph
    size_t countMergeSort; // number of reads passing merge sort
    size_t countAligner; // number of reads passing aligner
    CountStats() {
        countMinHash = countMinHashNotInGraph = countMergeSort = countAligner = 0;
    }
    CountStats operator+(const CountStats &c) {
        CountStats cs;
        cs.countMinHash = countMinHash + c.countMinHash;
        cs.countMinHashNotInGraph = countMinHashNotInGraph + c.countMinHashNotInGraph;
        cs.countMergeSort = countMergeSort + c.countMergeSort;
        cs.countAligner = countAligner + c.countAligner;
        return cs;
    }
};

struct match_info
{
    read_t vid;//在rD类中下标id
    read_t id;//序列的唯一id
    std::vector<Edit> editScript;
    ssize_t beginOffset, endOffset;
    ssize_t pos;
    bool reverseComplement;
    size_t match_length;
};

class Consensus {
public:
    ReadData *rD;

    ReadFilter *rF;


    // Directory for storing the temp files (.genome, .pos, .type, .id metaData)
    std::string tempDir = "tempRaw/";

    std::string tempFileName = "Contig";

    int numThr;

    //parameters for minimap2
    size_t  k,m_k, m_w, max_chain_iter, edge_threshold;

    /**
     * @brief Generates consensus, calls writeReads and writeMainPath on each
     * of the consensus graphs, and combines their output
     */
    void generateAndWriteConsensus(size_t loopindex);

    /**
     * @brief Combine files from threads and write metadata file
     *
     * @param numReadsInContig Vector with number of reads per contig in each thread
     */
    void finishWriteConsensus(const std::vector<std::vector<read_t>>& numReadsInContig);

    /**
     * @brief Checks that the read read in cG is equal to the read in ReadData
     * rD
     *
     * Assumes read is contained inside cG
     * @param cG
     * @param read
     * @return true
     * @return false
     */
    bool checkRead(ConsensusGraph *cG, read_t read);

    Consensus();

private:
    /**
     * @brief vector for checking repetitives
     *
     */
    std::vector<uint8_t> isRepetitive;
    /**
     * @brief Whether the reads have been added to a graph
     *
     */
    std::vector<uint8_t> inGraph;

    /**
     * @brief Number of locks to use (read i protected by lock i%numLocks)
     * Use large enough number of locks to avoid lock contention
     */
    static const uint32_t numLocks = (1<<24); // 16777216 locks

    /**
     * @brief Protects inGraph
     *
     * Note this is an array of locks, with read i protected by lock i%numLocks
     */
    std::vector<OmpMutex> readStatusLock;

    read_t numReads;

    /**
     * @brief Initializes numReads, readStatusLock
     *
     */
    void initialize();

    /**
     * @brief Gets an unadded read and updates its status to added
     *
     * @param read - starting point for searching unadded reads
     * @return true
     * @return false
     */
    bool getRead(read_t &read);

    /**
     * @brief Add reads overlapping with the mainPath of cG at position curPos
     *
     * @param cG
     * @param curPos
     * @param len The length of the mainPath used to get filtered reads
     * @param cs for collecting stats
     * @param logfile for writing pass/fail info to log
     * @param contigId current contig id for logging purposes
     */
    void addRelatedReads(ConsensusGraph *cG, ssize_t curPos, int len, CountStats &cs, std::ofstream &logfile, int contigId);

    /***
     * @brief 给cG的主链，找寻相似序列进行合并。
     * @param cG
     * @return 如果能找到就是true,否则FALSE
     */
    bool addRelatedReads(ConsensusGraph *cG);
    /**
     * @brief Create consensus graph by picking previously unadded read
     *
     * @param firstUnaddedRead starting point for searching unadded reads
     */
    ConsensusGraph *createGraph(read_t &firstUnaddedRead);

    /**
     * @brief check if a read string is reptitive or not
     *不同kmer的个数
     * @param readID the read id
     */
    bool checkRepetitive(read_t readID);
};

template<typename T>
std::shared_ptr<std::vector<T>> CombineVectors(std::shared_ptr<std::vector<T>> &&a, std::shared_ptr<std::vector<T>> &&b) {
    std::shared_ptr<std::vector<T>> ab = move(a);
    (*ab).insert(
            ab->end(), std::make_move_iterator(b->begin()), std::make_move_iterator(b->end()));
    return ab;
}
template<typename T>
std::vector<T> CombineVectors(const std::vector<T> &a, const std::vector<T> &b) {
    std::vector<T> ab = a;
    ab.insert(ab.end(), b.begin(), b.end());
    return ab;
}
#endif //MY_ORDERHASH_COMPRESS_CONSENSUS_H
