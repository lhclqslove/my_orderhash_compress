//
// Created by USER on 2022/5/31.
//
#include "Read.h"

namespace my_read{
    Read::Read(const std::string &s,size_t id,size_t cnt):id(id) ,cnt(cnt){
        read=std::unique_ptr<DnaBitset>(new DnaBitset(s.c_str(),s.size()));
        w=std::vector<size_t>(s.size()==0?0:s.size()-1,1);
    }
    Read::Read(size_t id):id(id){}
    Read::Read() : id(0) {
        read= nullptr;
    }

}
