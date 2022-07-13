//
// Created by USER on 2022/5/31.
//

#ifndef MY_ORDERHASH_COMPRESS_READ_H
#define MY_ORDERHASH_COMPRESS_READ_H

#include <memory>
#include "dnaToBits.h"
/**
 * @Description:
 * @Author LHC
 * @Date: 2022/5/31 16:54
 * @Version 1.0
 */
namespace my_read
{
    class Read
    {
    public:
        Read(const std::string &s,size_t id,size_t cnt);
        Read(size_t id);
        Read();
        const size_t id;
        std::unique_ptr<DnaBitset> read;//序列本身每一个字符用两位表示
        std::vector<short> w;//序列中间每条边的权值
        size_t cnt;//余下的查询次数
    };

}

#endif //MY_ORDERHASH_COMPRESS_READ_H
