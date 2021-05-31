#pragma once
#include <vector>
#include "Global.h"
using namespace std;

extern VDouble3D rhs;

class Spatial_Derivative
{
public:
	Spatial_Derivative();
	~Spatial_Derivative(){}

public:
	void Compute_Spatial_Derivative();
protected:
	int ist, ied, jst, jed;

	void Spatial_Derivative_X();
	void Spatial_Derivative_Y();
};