//
// Created by USER on 2022/5/31.
//
#include "Read.h"
Read::Read(const std::string &s,int id,size_t cnt):id(id) ,cnt(cnt){
    read=std::unique_ptr<DnaBitset>(new DnaBitset(s.c_str(),s.size()));
    w=std::vector<size_t>(s.size()-1,1);
}