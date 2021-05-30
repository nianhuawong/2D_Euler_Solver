#pragma once
#include <vector>
using namespace std;

class Residual
{
public:
	Residual();
	~Residual() {}
protected:
	vector< double > res_L1;
	vector< double > res_L2;
	vector< double > res_Loo;

public:
	void Compute_Residual();

protected:
	void OutputResidual();
	bool Stop_by_Residual();

};