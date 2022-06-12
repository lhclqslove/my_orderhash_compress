//
// Created by USER on 2022/5/31.
//
#ifndef MY_ORDERHASH_COMPRESS_READDATA_H
#define MY_ORDERHASH_COMPRESS_READDATA_H
#include "Read.h"
#include "Types.h"
class ReadData
{
    public:
        /**
         * @brief The filetype to read
         *
         */
        enum Filetype { FASTQ, READ, GZIP};
        /** The average length of the reads **/
        size_t avgReadLen;
        static std::atomic<read_t> index;//用于给新合并的序列创建id
        /** The length of the longest read **/
        size_t maxReadLen;
        std::string tempDir;
        size_t que_cnt;//查询挽回次数
        static size_t sequence_number_threshold;//序列数阈值，当数据集的序列数小于阈值时不在合并，设置为原有序列10%
        /// 加载序列文件存入DNAbitset
        /// \param filename
        /// \param filetype

        void loadFromFile(const char *filename,enum  Filetype filetype=FASTQ);

        /**
        * @brief Turns a base to its complement base
        *
        * A <-> T
        * C <-> G
        *
        * @param base
        * @return char
        */
        static char toComplement(char base);
        /**
        * @brief Get the number of reads
        *
        * @return read_t
        */
        read_t getNumReads();
        /**
        * @brief Returns the read with readId as a string
        *
        * @param readId, readStr
        * @return none
        */
        void getRead(read_t readId, std::string &readStr);
        /***
         * @brief 以智能指针的形式获取my_read::Read对象。
         * @param readId
         * @param ptr
         */
        void getRead(read_t readId,std::shared_ptr<my_read::Read> &ptr);

        void getindex(read_t readId,read_t &index);
        /**
         * @brief Turns a DNA strand to its reverse complement.
         *
         * @tparam Iterator Bidirectional random access iterator that dereferences
         * to char; in particular, dereferences to one of 'A', 'T', 'C', and 'G'
         * @tparam Inserter Insert iterator used to store the reverse complement
         * DNA. Should dereference to char.
         * @param originalBegin Beginning of original string
         * @param originalEnd End of original string (one past end)
         * @param reverseComplement Insert iterator used to store the reverse
         * complement string
         */
        template <typename Iterator, typename Inserter>
        static void toReverseComplement(Iterator originalBegin,
                                        Iterator originalEnd,
                                        Inserter reverseComplement);
        /* Destructor to delete bitset file when applicable */
    private:


        read_t numReads;
        /*存储序列文件到内存中*/
        std::vector<std::shared_ptr<Read>> readData;
        void loadFromFastqFile(const char *fileName,bool gzip_flag);

};
/******************************************************************************/
/* Implementation of public template functions */
template <typename Iterator, typename Inserter>

void ReadData::toReverseComplement(Iterator originalBegin, Iterator originalEnd,
                                   Inserter reverseComplement) {
    if (originalBegin == originalEnd)
        return;
    Iterator it = originalEnd;
    do {
        --it;
        *(reverseComplement++) = toComplement(*it);
    } while (it != originalBegin);
}

#endif //MY_ORDERHASH_COMPRESS_READDATA_H
