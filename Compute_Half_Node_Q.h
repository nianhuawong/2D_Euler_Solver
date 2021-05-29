#pragma once

class Half_Node_Q_Solver
{
public:
	Half_Node_Q_Solver();
	~Half_Node_Q_Solver() {};
protected:
	int M_Dim, N_Dim;
	int method_of_half_q;
	int method_of_limiter;
	double muscl_k;
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