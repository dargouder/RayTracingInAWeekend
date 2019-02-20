cmake_minimum_required(VERSION 2.6)
project (RTIAW)

find_package(OpenMP)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
file(GLOB_RECURSE SOLS_HEADERS "source/*.h")

set (SOLS_INCLUDE_DIRS "")
foreach (_headerFile ${SOLS_HEADERS})
    get_filename_component(_dir ${_headerFile} PATH)
    list (APPEND SOLS_INCLUDE_DIRS ${_dir})
endforeach()
list(REMOVE_DUPLICATES SOLS_INCLUDE_DIRS)



include_directories(${SOLS_INCLUDE_DIRS})

file(GLOB_RECURSE SOURCES "source/*.cpp")

add_executable(RTIAW ${SOURCES})

target_link_libraries(RTIAW PRIVATE "${OpenMP_CXX_FLAGS}")
target_compile_options(RTIAW PRIVATE "${OpenMP_CXX_FLAGS}")