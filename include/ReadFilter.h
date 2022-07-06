//
// Created by USER on 2022/6/1.
//

#ifndef MY_ORDERHASH_COMPRESS_READFILTER_H
#define MY_ORDERHASH_COMPRESS_READFILTER_H
#include "ReadData.h"
#include "Types.h"
#include "BBHashMap.h"
#include <xxhash.hpp>


struct mer_info {
    size_t pos;
    uint64_t kmer;
    uint64_t hash;
    unsigned occ;
//    mer_info(size_t p, unsigned o, uint64_t h,uint64_t mer=0)
//            : pos(p)
//            , hash(mer)
//            , occ(o)
//            , kmer(h)
//    { };
    void init(size_t p, unsigned o, uint64_t mer){
        pos=p;
        occ=o;
        kmer=mer;
    };
    mer_info(){}
};
/**
 * @brief Filters reads that are likely to overlap
 *
 */
class ReadFilter {
public:
    /**
     * @brief Rechieves reads that are likely to overlap with string s and
     * stores them in results
     *
     * @param s string to find overlapping reads with
     * @param results Vector where we store the results
     */
    virtual void getFilteredReads(const std::string &s,
                                  std::vector<read_t> &results,size_t cnt) = 0;//纯虚函数，必须由派生类提供定义,之后才可以实例化

    virtual void initialize(ReadData &rD) = 0;

    virtual ~ReadFilter();
};

class OrderHashReadFilter : public ReadFilter {
public:

    /** [k]-mer **/
    size_t k;
    /** size of sketch **/
    size_t n;
    size_t l;
    static size_t max_occ;
    static size_t que_cnt;
    std::string tempDir;

    /**
     * @brief Builds the hash tables from the data in rD
     *
     * @param rD
     */
    void initialize(ReadData &rD) override;

    /**
     * @brief Generates a sequence of n kMer_t random numbers
     *
     * @param n
     */
    void generateRandomNumbers(size_t n);

    /**
     * @brief Converts a string to a sketch and stores it in sketch
     *
     * @param s
     * @param sketch
     * @param kMers preallocated for speed during multithreaded initialization
     * (size should be at least s.size()-k+1)
     * @param hashes preallocated for speed during multithreaded initialization
     * (size should be at least n*(s.size()-k+1))
     */
//    void string2Sketch(const std::string &s, kMer_t *sketch,
//                       std::vector<kMer_t> &kMers, std::vector<kMer_t> &hashes);
    /**
     * @Author Lhc
     * @Description //TODO 
     * @Date 10:15 2022/6/2
     * @Param
     * @return 
     * @return null
     **/
    void string2Sketch(const std::string &s, kMer_t *sketch,
                       std::vector<mer_info> &kMers,xxhash &hash);
//    /**
//     * @brief Calculates n hashes of kMer and stores them in hashes
//     *
//     * @param kMer
//     * @param hashes
//     */
//    void hashKMer(const kMer_t kMer, std::vector<kMer_t> &hashes);
    /**
     * @brief 查询与s相同hash值最多的几个，警告：必须保证s满足不同kmer个数大于等于l
     * @param s
     * @param results
     */
    void getFilteredReads(const std::string &s,
                          std::vector<read_t> &results,size_t cnt) override;


    OrderHashReadFilter();

    ~OrderHashReadFilter() override;

    /**
     * Turns a k-mer in string format to an int
     * @param s k-mer
     * @return int representing k-mer
     */
    static kMer_t kMerToInt(const std::string &s);

    /**
     * @brief Converts a base, one of 'A', 'T', 'C', 'G' into a two-bit int
     *
     * @param base
     * @return char
     */
    static char baseToInt(const char base);

    /**
     * @brief Converts a string to kmers and stores it in kMers
     *
     * @param s
     * @param k
     * @param kMers
     */
//    static void string2KMers(const std::string &s, const size_t k,
//                             std::vector<kMer_t> &kMers);
    /// @brief Converts a string to kmers and stores it in kMer_info with kmer's occ<maxocc
    /// \param s
    /// \param k
    /// \param KMers_info
    ///\return 如果挑选能超过l个满足条件的kmer返回true,否则false
    static bool  string2KMers(const std::string &s,const size_t k,const size_t l,std::vector<mer_info> &KMers_info);
private:
    ReadData *rD;
//    std::vector<size_t> *readPos;
//    std::vector<std::pair<size_t, read_t>> readPosSorted;
    read_t numReads;
    std::vector<BBHashMap> hashTables; // vector of size n to store the n hash maps
    kMer_t *randNumbers = nullptr;

    std::hash<kMer_t> hasher;

    /**
     * @brief Initializes the hash tables given the data calculated in
     * sketches
     *
     * @param sketches
     */
    void populateHashTables(const std::vector<kMer_t> &sketches);

    /**
     * @brief Get the Filtered Reads likely to overlap with the read represented
     * by sketch and stores them in results
     *
     * @param sketch
     * @param
     * @param results
     */
    void getFilteredReads(kMer_t sketch[], std::vector<read_t> &results,size_t cnt);

};

#endif //MY_ORDERHASH_COMPRESS_READFILTER_H
