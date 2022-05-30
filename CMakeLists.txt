cmake_minimum_required(VERSION 3.23)
project(my_orderhash_compress VERSION 0.1
        LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)
include(CheckIPOSupported)
include(CheckCXXSourceCompiles)

# Clear All Flags
set(CMAKE_CXX_FLAGS "")
set(CMAKE_CXX_FLAGS_DEBUG "")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "")
set(CMAKE_CXX_FLAGS_RELEASE "")

#release 版本
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
#判断是否支持ipo编译优化
check_ipo_supported(RESULT result)
if(result)
    message(STATUS "Support for Interprocedural Optimization detected.")
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")#构建类型
message(STATUS "CHECKS=${CHECKS} (Use -DCHECKS=True to enable all checks")
message(STATUS "LOG=${LOG} (Use -LOG=True to output log files for analysis")#log条件编译
message(STATUS "Directory: ${CMAKE_CURRENT_SOURCE_DIR}")
#openmp条件库，用于多线程
find_package(OpenMP)

##############################
# malloc_trim用于运行时检测内存使用情况
check_cxx_source_compiles(
        "#include <malloc.h>
     int main() {
         malloc_trim(0);
         return 0;
     }"
        USE_MALLOC_TRIM)
if(USE_MALLOC_TRIM)
    message(STATUS "Support for malloc_trim detected.")
else()
    message(STATUS "Support for malloc_trim not detected.")
endif()

###############################
###############################
# boost
set( Boost_NO_SYSTEM_PATHS ON ) # do not use system boost
add_subdirectory(boost-cmake)

#######################################################
#######################################################
# libbsc 通用压缩器，用于之后的通用压缩
add_library(libbsc
        STATIC
        libbsc/adler32/adler32.cpp
        libbsc/bwt/libsais/libsais.c
        libbsc/bwt/bwt.cpp
        libbsc/coder/coder.cpp
        libbsc/coder/qlfc/qlfc.cpp
        libbsc/coder/qlfc/qlfc_model.cpp
        libbsc/filters/detectors.cpp
        libbsc/filters/preprocessing.cpp
        libbsc/libbsc/libbsc.cpp
        libbsc/lzp/lzp.cpp
        libbsc/platform/platform.cpp
        libbsc/st/st.cpp
        )


add_executable(my_orderhash_compress main.cpp)
