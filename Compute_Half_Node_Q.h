#pragma once

class Half_Node_Q_Solver
{
public:
	Half_Node_Q_Solver();
	~Half_Node_Q_Solver() {};
protected:

public:
	void Half_Node_Q();
protected:
	double Limiter_Function(double ita);
	void Half_Node_Q_MUSCL();
	void Half_Node_Q_WCNS ();

	void Half_Node_Q_MUSCL_Y();

protected:
	double minmod_limiter(double a, double b);
	double vanleer_limiter(double a, double);
	double superbee_limiter(double a, double);
};

void Half_Node_Q();