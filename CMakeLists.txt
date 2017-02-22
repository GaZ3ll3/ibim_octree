cmake_minimum_required(VERSION 2.8)
project(ibim_octree)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_RELEASE} -g -Ofast -DRUN_OMP -march=native -std=c++11 -fopenmp ")
set(SOURCE_FILES main.cpp Vector.cpp Vector.h Matrix.cpp Matrix.h linalg.h linalg.cpp molecule.h point.h node.h tree.cpp grid.cpp grid.h levelset.cpp levelset.h)
#set(SOURCE_FILES main.cpp Vector.cpp Vector.h Matrix.cpp Matrix.h linalg.h linalg.cpp molecule.h point.h node.h tree.cpp view.cpp view.h grid.cpp grid.h levelset.cpp levelset.h)

#set(OPENBLAS_PATH "/home/lurker/Workspace/lib/openblas")

#include_directories("${OPENBLAS_PATH}/include")
#link_directories("${OPENBLAS_PATH}/lib")
add_executable(ibim_octree ${SOURCE_FILES})
#target_link_libraries(ibim_octree  m "${OPENBLAS_PATH}/lib/libopenblas.so")
target_link_libraries(ibim_octree  m blas)

#target_link_libraries(ibim_octree  m GL GLU glut /home/lurker/Workspace/lib/openblas/lib/libopenblas.so)