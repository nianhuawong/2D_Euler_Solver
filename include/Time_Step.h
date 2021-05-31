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