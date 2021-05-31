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

void Update_Flowfield();
void Update_Flowfield_X();
void Update_Flowfield_Y();