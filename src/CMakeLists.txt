cmake_minimum_required(VERSION 3.7)
project(varViewer VERSION 0.2.7 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if (NOT CMAKE_BUILD_TYPE)
    set(DEFAULT_CMAKE_BUILD_TYPE Release)
    set(CMAKE_BUILD_TYPE ${DEFAULT_CMAKE_BUILD_TYPE} CACHE STRING
            "Choose the type of build (default: ${DEFAULT_CMAKE_BUILD_TYPE})" FORCE)
endif ()

set(Boost_USE_STATIC_LIBS ON)
find_package(Boost 1.73 REQUIRED COMPONENTS program_options filesystem system)

find_package(ZLIB REQUIRED)
find_package(BZip2 REQUIRED)
find_package(LibLZMA REQUIRED)
find_package(CURL REQUIRED)
find_library(htslib libhts.a "/ifstj2/B2C_RD_H1/Personal/huangfei/software/bin/lib")
find_library(htslib hts "/ifstj2/B2C_RD_H1/Personal/huangfei/software/bin/lib")

add_executable(varViewer
        varViewer.cpp
        utils.hh utils.cpp
        Variant.hh Variant.cpp
        Reference.hh Reference.cpp
        Operation.hh Operation.cpp
        LinearAlign.hh LinearAlign.cpp
	Aligns.hh Aligns.cpp
        LanePlot.hh LanePlot.cpp
        GenerateSvg.hh GenerateSvg.cpp)

target_include_directories(varViewer PUBLIC
        ${CMAKE_SOURCE_DIR}
        ${LIBLZMA_INCLUDE_DIRS}
        ${CURL_INCLUDE_DIRS}
        ${Boost_INCLUDE_DIRS}
	"/ifstj2/B2C_RD_H1/Personal/huangfei/software/bin/include"
        )

# Set static linking of gcc standard libraries to simplify binary distribution
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(STATIC_FLAGS -static-libgcc -static-libstdc++)
endif ()

target_link_libraries(varViewer PUBLIC
        ${STATIC_FLAGS}
        ${htslib}
        ${LIBLZMA_LIBRARIES}
        ${CURL_LIBRARIES}
        ZLIB::ZLIB
        BZip2::BZip2)


install(TARGETS varViewer RUNTIME DESTINATION bin)
