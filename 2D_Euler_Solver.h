// 2D_Euler_Solver.h: 标准系统包含文件的包含文件
// 或项目特定的包含文件。

#pragma once

#include <iostream>

class Simulation
{
public:
	Simulation() {}
	~Simulation(){}

public:
	void Run();

};
int iter, numberOfTimeSteps;

void Init_Global_Param() {};
void Generate_Mesh() {};
void Flow_Initialization() {};

// TODO: 在此处引用程序需要的其他标头。
