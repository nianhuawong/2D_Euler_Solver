# CMakeList.txt: 2D_Euler_Solver 的 CMake 项目，在此处包括源代码并定义
# 项目特定的逻辑。
#
cmake_minimum_required (VERSION 3.8)

project ("2D_Euler_Solver")

# 将源代码添加到此项目的可执行文件。
add_executable (2D_Euler_Solver "2D_Euler_Solver.cpp" "2D_Euler_Solver.h" "Compute_Half_Node_Q.cpp" "Global_Variables.h" "Compute_Half_Node_Flux.h" "Compute_Half_Node_Flux.cpp" "Global_Variables.cpp" "Geometry.cpp" "Geometry.h")

# TODO: 如有需要，请添加测试并安装目标。
