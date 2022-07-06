//
// Created by USER on 2022/6/1.
//
#include <random>
#include "ReadFilter.h"
#include <algorithm>
#include <chrono>
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
        std::vector<mer_info> kMersVec(maxNumkMers);
        std::string readStr;
        xxhash hash;
#pragma omp for
        for (read_t i = 0; i < numReads; ++i){
            rD.getRead(i, readStr);
            string2Sketch(readStr, sketches.data() + i * n, kMersVec, hash);
        }
    } // pragma omp parallel
    populateHashTables(sketches);
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
void OrderHashReadFilter::string2Sketch(const std::string &s, kMer_t *sketch, std::vector<mer_info> &kMers,
                                        xxhash &hash) {
    ssize_t numKMers = s.length() - k + 1;
    if (numKMers < 0)
        return;
//    std::cout<<s.substr(0,5)<<" "<<"yes"<<std::endl;
    if(!string2KMers(s, k,l, kMers))//如果没有至少L个少于max_occ的kmer,认为自重复后面直接压缩。
    {
        return;
    }
    //为每个hash函数挑选前l个最小的mer_info

    for(size_t i=0;i<n;i++) {
        //通过xxhash计算hash值
        for (size_t j = 0; j < numKMers; j++) {
            //还是只hashkmer
            kMers[j].hash=hasher((randNumbers[i]^kMers[j].kmer)+kMers[j].occ);
//            kMers[j].hash=hasher((randNumbers[i]^kMers[j].kmer));
//            //初始化随机种子
//            hash.reset(randNumbers[i]);
//            //拆入kmer和occ
//            hash.update(&kMers[j].kmer, sizeof(kMers[j].kmer));
//            hash.update(&kMers[j].occ, sizeof(kMers[j].occ));
//            //计算hash值
//            kMers[j].hash = hash.digest();
        }

        //按hash值排序
        std::partial_sort(kMers.begin(), kMers.begin() + l, kMers.begin() + numKMers,
                          [&](const mer_info &x, const mer_info &y) { return x.hash < y.hash; });
//        sort_mer_info(kMers,l,numKMers);

        //前l个按pos排序
        std::sort(kMers.begin(), kMers.begin() + l,
                  [&](const mer_info &x, const mer_info &y) { return x.pos < y.pos; });
        //把前l个mer_info哈希成一个值，等效于minhash中的某一个最小hash值
        hash.reset(randNumbers[i]);
        for (size_t j = 0; j < l; j++) {
            hash.update(&kMers[j].kmer, sizeof(kMers[j].kmer));
//            hash.update(&kMers[j].occ, sizeof(kMers[j].occ));
//            std::cout<<kMers[j].kmer<<" "<<kMers[j].occ<<std::endl;
        }
        sketch[i] = hash.digest();
    }
}
bool  OrderHashReadFilter::string2KMers(const std::string &s, const size_t k,const size_t l, std::vector<mer_info> &KMers_info) {
    ssize_t maxI = s.length() - k + 1;
    if (maxI <= 0)
        return false;
    std::unordered_map<kMer_t, unsigned> occurrences;
    kMer_t currentKMer = kMerToInt(s.substr(0, k));
    KMers_info[0].init(0,0,currentKMer);
    occurrences[currentKMer]=1;
    const unsigned long long mask = (1ull << (2 * k)) - 1;
     
    for (size_t i = 1; i < (size_t)maxI; ++i) {
        currentKMer =
                ((currentKMer << 2) | OrderHashReadFilter::baseToInt(s[i + k - 1])) &
                mask;
        auto occ = occurrences[currentKMer]++;

        KMers_info[i].init(i,occ,currentKMer);
    }
    //重新扫一遍，当前位置kmer_occ小于max_occ ,与pos交换
    int pos=0;
    for(size_t i=0;i<(size_t)maxI;++i)
    {
        if(occurrences[KMers_info[i].kmer]<=max_occ)
        {
            if(pos==i)pos++;
            else{
                std::swap(KMers_info[i],KMers_info[pos++]);
            }
        }
    }
    KMers_info.erase(KMers_info.begin()+pos,KMers_info.end());
    if(KMers_info.size()<l)
    {
        return false;
    }
    else return true;
}
void OrderHashReadFilter::populateHashTables(const std::vector<kMer_t> &sketches) {
    std::cout << "Starting to populate hash tables" << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    hashTables.resize(n);
#pragma omp parallel for
    for (size_t i = 0; i < n; ++i)
        hashTables[i].initialize(n,numReads,sketches.data(),i,tempDir);

    auto end = std::chrono::high_resolution_clock::now();
    auto duration =
            std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "finished populating hash tables" << std::endl;
    std::cout << duration.count() << " milliseconds passed" << std::endl;
}

void OrderHashReadFilter::getFilteredReads(const std::string &s, std::vector<read_t> &results, size_t cnt) {
    results.clear();
    std::vector<kMer_t> sketch(n);
    size_t numKmers;
    if (s.size() < k - 1)
        numKmers = 0;
    else
        numKmers = s.size() - k + 1;
    std::vector<mer_info> kMersVec(numKmers);
    xxhash hash;
//    std::vector<kMer_t> hashesVec(n); // preallocation
    string2Sketch(s, sketch.data(), kMersVec,hash);

    getFilteredReads(sketch.data(), results,cnt);
}
void OrderHashReadFilter::getFilteredReads(kMer_t *sketch, std::vector<read_t> &results, size_t cnt) {
    std::vector<read_t> matches;
    std::vector<std::pair<size_t,read_t>> tmp_matches;//first:相同哈希值个数，second序列id
    results.clear();
    for (size_t sketchIndex = 0; sketchIndex < n; ++sketchIndex) {
        kMer_t curHash = sketch[sketchIndex];
        hashTables[sketchIndex].pushMatchesInVector(curHash,matches);
    }
    std::sort(matches.begin(), matches.end());
    auto end = matches.end();
    auto next = matches.begin();
    //获取具有相同hash值的候选id,以及相同hash值的个数
    for (auto it = matches.begin(); it != end; it = next) {
        next = std::upper_bound(it, end, *it);
        tmp_matches.push_back(std::move( std::pair<size_t,read_t>(next - it,*it)));
    }
    //部分排序，挑出min(cnt,size)个相同hash值个数最大的read
    std::partial_sort(tmp_matches.begin(),tmp_matches.begin()+std::min(cnt,tmp_matches.size()),tmp_matches.end(),[&](const std::pair<size_t,read_t> &x, const std::pair<size_t,read_t> &y) { return x.first < y.first; });
    for(size_t i=0;i<std::min(cnt,tmp_matches.size());i++)
    {
        results.push_back(tmp_matches[i].second);
    }
}
OrderHashReadFilter::OrderHashReadFilter() {}
OrderHashReadFilter::~OrderHashReadFilter() {
    if (randNumbers)
        delete[] randNumbers;

}
size_t OrderHashReadFilter::que_cnt;
size_t OrderHashReadFilter::max_occ;
kMer_t OrderHashReadFilter::kMerToInt(const std::string &s) {
    size_t l = s.length();
    kMer_t result = 0;
    for (size_t i = 0; i < l; ++i) {
        result <<= 2;
        result |= baseToInt(s[i]);
    }
    return result;
}
// Using the bit operations version of this function provides a 13X improvement
// in speed
char OrderHashReadFilter::baseToInt(const char base) {
    return (base & 0b10) | ((base & 0b100) >> 2);
}