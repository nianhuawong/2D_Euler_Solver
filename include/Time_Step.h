//////////////////////////////////////////////////////////
//														//
//			2D Euler Solver on Structured Grid			//
//					by Wang Nianhua						//
//				Email: nianhuawong@126.com				//
//					   2021.6.7							//
//////////////////////////////////////////////////////////
#pragma once 

class Time_Step
{
public:
	Time_Step();
	~Time_Step(){}

protected:
	int ist, ied, jst, jed;
public:
	void Compute_Time_Step();

};