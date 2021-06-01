#pragma once
#include <vector>
#include "Global.h"
using namespace std;

extern VDouble3D qField;
extern VDouble3D qField1;
extern VDouble3D qField2;
extern VDouble3D qField_N0;
extern VDouble3D qField_N1;

class QlQr_Solver
{
public:
	QlQr_Solver();
	~QlQr_Solver() {};
protected:
	int ist, ied, jst, jed;
public:
	void Solve_QlQr();
protected:
	double Limiter_Function(double ita);

	void QlQr_MUSCL();
	void QlQr_MUSCL_X();
	void QlQr_MUSCL_Y();

	void Boundary_QlQr_MUSCL_X();
	void Boundary_QlQr_MUSCL_Y();

	void QlQr_WCNS();
	void QlQr_WCNS_X();
	void QlQr_WCNS_Y();

protected:

};

double minmod_limiter(double a, double b);
double vanleer_limiter(double a, double);
double superbee_limiter(double a, double);

