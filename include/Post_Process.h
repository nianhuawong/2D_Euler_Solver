//////////////////////////////////////////////////////////
//														//
//			2D Euler Solver on Structured Grid			//
//					by Wang Nianhua						//
//				Email: nianhuawong@126.com				//
//					   2021.6.7							//
//////////////////////////////////////////////////////////
#pragma once
#include <vector>
#include "Global.h"
using namespace std;

class Residual
{
public:
	Residual();
	~Residual() {}
protected:
	VDouble res_L1;
	VDouble res_L2;
	VDouble res_Loo;
	VInt2D max_index;
	int ist, ied, jst, jed;
public:
	void Compute_Residual();

protected:
	void OutputResidual();
	bool Stop_by_Residual();

};