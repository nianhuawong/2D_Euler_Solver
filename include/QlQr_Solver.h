#pragma once
#include <vector>
using namespace std;

extern vector< vector< vector< double > > > qField;
extern vector< vector< vector< double > > > qField1;
extern vector< vector< vector< double > > > qField2;
extern vector< vector< vector< double > > > qField_N1;
extern vector< vector< vector< double > > > qField_N2;
extern vector< vector< vector< double > > > qField_N3;
class QlQr_Solver
{
public:
	QlQr_Solver();
	~QlQr_Solver() {};
protected:

public:
	void Solve_QlQr();
protected:
	double Limiter_Function(double ita);
	void QlQr_MUSCL();
	void QlQr_WCNS();

	void QlQr_MUSCL_X();
	void QlQr_MUSCL_Y();

protected:

};

double minmod_limiter(double a, double b);
double vanleer_limiter(double a, double);
double superbee_limiter(double a, double);

