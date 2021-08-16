#CMakeLists.txt
cmake_minimum_required(VERSION 2.6)
project(pplsd)

find_package(MPI REQUIRED)

include_directories(${MPI_INCLUDE_PATH})

aux_source_directory(. DIR_SRCS)
add_executable(pplsd ${DIR_SRCS})
target_link_libraries(pplsd ${MPI_LIBRARIES})

if(MPI_COMPILE_FLAGS)
  set_target_properties(pplsd PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(pplsd PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()


