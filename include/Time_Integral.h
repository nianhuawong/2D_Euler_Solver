//////////////////////////////////////////////////////////
//														//
//			2D Euler Solver on Structured Grid			//
//					by Wang Nianhua						//
//				Email: nianhuawong@126.com				//
//					   2021.6.7							//
//////////////////////////////////////////////////////////
#pragma once

class Time_Marching_Solver
{
public:
	Time_Marching_Solver();
	~Time_Marching_Solver() {}

public:
	void Time_Marching();
	
protected:
	
};

void Update_Flowfield(int iStage);
void Update_Flowfield_X(int iStage);
void Update_Flowfield_Y(int iStage);
void SolutionFix(VDouble& primitiveVector, int i, int j);