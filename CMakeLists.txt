# CMakeList.txt: 2D_Euler_Solver 的 CMake 项目，在此处包括源代码并定义
# 项目特定的逻辑。
#
cmake_minimum_required (VERSION 3.8)

project ("2D_Euler_Solver")

#set(CMAKE_CXX_COMPILER  g++)

#include(CheckCXXCompilerFlag)
#CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
#CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)

#if(COMPILER_SUPPORTS_CXX11)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
#elseif(COMPILER_SUPPORTS_CXX0X)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
#else()
#message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
#endif()

set(CMAKE_CXX_FLAGS    -openmp)
#set(OMP_NUM_THREADS   4)

#包含目录
include_directories("./include")


set(CMAKE_CXX_STANDARD 11)
# 将源代码添加到此项目的可执行文件。
add_executable (2D_Euler_Solver "src/2D_Euler_Solver.cpp" "include/2D_Euler_Solver.h" "src/QlQr_Solver.cpp" "include/Global.h" "include/Flux_Solver.h" "src/Flux_Solver.cpp" "src/Global.cpp" "src/Geometry.cpp" "include/Geometry.h" "src/Spatial_Derivative.cpp" "include/Spatial_Derivative.h" "include/Time_Integral.h" "src/Time_Integral.cpp" "src/Post_Process.cpp" "include/Post_Process.h" "include/Time_Step.h" "src/Time_Step.cpp" "src/Compute_Boundary.cpp" "src/Init_Flow.cpp")

# TODO: 如有需要，请添加测试并安装目标。
