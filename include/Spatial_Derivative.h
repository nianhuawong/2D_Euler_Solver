#pragma once
#include <vector>
using namespace std;

extern vector< vector< vector< double > > > rhs;

class Spatial_Derivative
{
public:
	Spatial_Derivative();
	~Spatial_Derivative(){}

public:
	void Compute_Spatial_Derivative();
protected:
	
	void Spatial_Derivative_X();
	void Spatial_Derivative_Y();
};