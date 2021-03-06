//
// Created by USER on 2022/5/30.
//
#include "DirectoryUtils.h"
#include <boost/algorithm/string/predicate.hpp>

namespace DirectoryUtils {

    uint32_t read_var_uint32(std::ifstream &fin) {
        uint32_t val = 0;
        uint8_t byte;
        uint8_t shift = 0;
        do {
            fin.read((char*)&byte, sizeof(uint8_t));
            val |= ((byte & 0x7f) << shift);
            shift += 7;
        } while(byte & 0x80);
        return val;
    }
    //每次都7位写出
    void write_var_uint32(const uint32_t val, std::ofstream &fout){
        uint32_t uval = val;
        uint8_t byte;
        while (uval > 127) {
            byte = (uint8_t)(uval & 0x7f) | 0x80;
            fout.write((char*)&byte, sizeof(uint8_t));
            uval >>= 7;
        }
        byte = (uint8_t)(uval & 0x7f);
        fout.write((char*)&byte, sizeof(uint8_t));
    }

    std::string random_string(size_t length){
        auto randchar = []() -> char {
            const char charset[] =
                    "0123456789"
                    "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                    "abcdefghijklmnopqrstuvwxyz";
            const size_t max_index = (sizeof(charset) - 1);
            return charset[rand() % max_index];
        };
        std::string str(length, 0);
        std::generate_n(str.begin(), length, randchar);
        return str;
    }

    void clearDir(const std::string &path) {
        boost::system::error_code ec;
        const boost::filesystem::path dirPath(path);
        boost::filesystem::remove_all(dirPath, ec);
        boost::filesystem::create_directory(dirPath, ec);
    }

    void combineFilesWithExt(const std::string &filestem,
                             const std::string &fileExt, const size_t numFiles,
                             bool delim) {
        std::ofstream outFile(filestem + fileExt, std::ios_base::binary);
        for (size_t i = 0; i < numFiles; i++) {
            std::string fileName = filestem + std::to_string(i) + fileExt;
            std::ifstream inFile(fileName, std::ios_base::binary);
            outFile << inFile.rdbuf();
            if (delim)
                outFile << ".\n";
            inFile.close();
            boost::system::error_code ec;
            const boost::filesystem::path oldFilePath(fileName);
            boost::filesystem::remove(oldFilePath, ec);
        }
        outFile.close();
    }

    void unpack(const std::string &filepath, const std::string &outputDir) {
        std::ifstream inFile(filepath);
        boost::filesystem::path path(filepath);
        const std::string stem = path.filename().stem().string();
        const std::string ext = path.filename().extension().string();
        std::ofstream outFile;
        std::string line;
        bool need2CreateFile = true;
        size_t i = 0;
        const std::string dirAndStem =
                (boost::algorithm::ends_with(outputDir, "/") ? outputDir
                                                             : outputDir + "/") +
                stem;
        while (std::getline(inFile, line)) {
            if (need2CreateFile) {
                need2CreateFile = false;
                outFile.open(dirAndStem + std::to_string(i) + ext);
                ++i;
                outFile << line << '\n';
            } else {
                if (line == ".") {
                    outFile.close();
                    need2CreateFile = true;
                } else
                    outFile << line << '\n';
            }
        }
    }

    void unpack(const std::string &filepath) {
        boost::filesystem::path path(filepath);
        unpack(filepath, path.parent_path().string());
    }


} // namespace DirectoryUtils

