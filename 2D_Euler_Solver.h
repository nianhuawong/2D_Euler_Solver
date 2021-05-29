#pragma once

class Simulation
{
public:
	Simulation() {}
	~Simulation(){}

public:
	void Run();
};
namespace GLOBAL {
	void Init_Global_Param();
	void Generate_Mesh();
	void Flow_Initialization();
	void Compute_Boundary();
	void Half_Node_Q();
	void Half_Node_Flux();
}

