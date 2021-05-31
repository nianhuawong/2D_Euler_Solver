# CMakeList.txt: 2D_Euler_Solver 的 CMake 项目，在此处包括源代码并定义
# 项目特定的逻辑。
#
cmake_minimum_required (VERSION 3.8)

project ("2D_Euler_Solver")

#包含目录
include_directories("./include")

set(MY_WORKING_DIR ./)

# 将源代码添加到此项目的可执行文件。
add_executable (2D_Euler_Solver "src/2D_Euler_Solver.cpp" "include/2D_Euler_Solver.h" "src/QlQr_Solver.cpp" "include/Global.h" "include/Flux_Solver.h" "src/Flux_Solver.cpp" "src/Global.cpp" "src/Geometry.cpp" "include/Geometry.h" "src/Spatial_Derivative.cpp" "include/Spatial_Derivative.h" "include/Time_Integral.h" "src/Time_Integral.cpp" "src/Post_Process.cpp" "include/Post_Process.h" "include/Time_Step.h" "src/Time_Step.cpp")

# TODO: 如有需要，请添加测试并安装目标。
