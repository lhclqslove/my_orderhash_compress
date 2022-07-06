cmake_minimum_required(VERSION 3.10)
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

set(LOG True)

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

#xxhash
find_package(xxHash 0.7 CONFIG REQUIRED)
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
target_include_directories(libbsc
        PUBLIC libbsc)

target_compile_options(libbsc
        PUBLIC
        -g -Wall#gcc编译器参数显示所有警告

        # Comment out CFLAGS line below for compatability mode for 32bit file sizes
        # (less than 2GB) and systems that have compilers that treat int as 64bit
        # natively (ie: modern AIX)
        -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64

        # Comment out CFLAGS line below to disable optimizations
        -O3 -fomit-frame-pointer -fstrict-aliasing -ffast-math

        # Comment out CFLAGS line below to disable AVX2 instruction set (performance will suffer)
        -mavx2

        # Comment out CFLAGS line below to disable OpenMP optimizations
        # We disable this here and parallelize at the top level
        # -fopenmp -DLIBBSC_OPENMP_SUPPORT

        # Comment out CFLAGS line below to enable debug output
        -DNDEBUG

        # Comment out CFLAGS line below to disable unaligned memory access
        -DLIBBSC_ALLOW_UNALIGNED_ACCESS
        )

if(OpenMP_CXX_FOUND)
    target_link_libraries(libbsc PUBLIC OpenMP::OpenMP_CXX -DLIBBSC_OPENMP_SUPPORT)
endif()


#######################################################
#minimap2

add_custom_target(
        MINIMAP2
        COMMAND make clean && make
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/minimap2/
)


#######################################################
#fast-lzma2
add_custom_target(
        FASTLZMA2
        COMMAND make clean && make
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/fast-lzma2/
)

#######################################################
add_library(springlib
        STATIC
        include/DirectoryUtils.h
        src/DirectoryUtils.cpp
        include/dnaToBits.h
        src/dnaToBits.cpp
        include/Read.h
        src/Read.cpp
        include/ReadData.h
        src/ReadData.cpp
        include/Types.h
        include/Compressor.h
        src/Compressor.cpp
        include/ReadFilter.h
        src/ReadFilter.cpp
        include/BBHashMap.h
        src/BBHashMap.cpp
        include/BooPHF.h
        include/xxhash.hpp
        include/SafeQueue.h
        include/ThreadPool.h
        include/Consensus.h
        src/Consensus.cpp
        include/ConsensusGraph.h
        src/ConsensusGraph.cpp
        include/bsc_helper.h
        src/bsc.cpp
        include/lzma2_helper.h
        src/lzma2.cpp
        include/OmpMutex.h
        src/OmpMutex.cpp
        include/Edits.h
        src/Edits.cpp
        )
target_include_directories(springlib
        PUBLIC include src
        PUBLIC libbsc
        PUBLIC minimap2 fast-lzma2
        )
target_compile_options(springlib
        PUBLIC
        -W -Wall -Wextra)

target_compile_options(springlib
        PUBLIC
        $<$<CONFIG:Release>:-O3 -msse4.1>)
target_compile_options(springlib
        PUBLIC
        $<$<CONFIG:RelWithDebInfo>:-g -O3 -msse4.1>)
target_compile_options(springlib
        PUBLIC
        $<$<CONFIG:Debug>:-g -O3 -msse4.1 -fsanitize=address -fno-omit-frame-pointer> -static-libasan)
if ( CHECKS )
    target_compile_options(springlib
            PUBLIC
            -DCHECKS)
endif()
if ( LOG )
    target_compile_options(springlib
            PUBLIC
            -DLOG)
endif()
if ( USE_MALLOC_TRIM )
    target_compile_options(springlib
            PUBLIC
            -DUSE_MALLOC_TRIM)
endif()
target_compile_definitions(springlib PUBLIC "$<$<CONFIG:Debug>:DEBUG>")
target_compile_definitions(springlib PUBLIC "$<$<CONFIG:RelWithDebInfo>:DEBUG>")

add_dependencies(springlib MINIMAP2)
add_dependencies(springlib FASTLZMA2)

target_link_libraries(springlib
        PUBLIC
        libbsc
        ${CMAKE_CURRENT_SOURCE_DIR}/minimap2/libminimap2.a
        ${CMAKE_CURRENT_SOURCE_DIR}/fast-lzma2/libfast-lzma2.a
        Boost::filesystem
        Boost::iostreams
        Boost::program_options
        -lz -ldl -lm -lpthread
        xxHash::xxhash
        )

target_link_libraries(springlib
        PUBLIC
        $<$<CONFIG:Debug>: -fsanitize=address -fno-omit-frame-pointer>)
# OpenMp support
if(OpenMP_CXX_FOUND)
    target_link_libraries(springlib PUBLIC OpenMP::OpenMP_CXX)
endif()
add_executable(my_orderhash_compress src/main.cpp)
target_link_libraries(my_orderhash_compress PUBLIC springlib)