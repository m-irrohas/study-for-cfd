#include<iostream>
using namespace std;

double get_e(double p, double rho, double u){
    static const double GAMMA = 1.4;
    double e = p/(GAMMA-1.0)+1.0/2.0*rho*u*u;
    return e;
}


double get_mid_q(double dt, double dx, double Q,  double E, double E_pre){
    double Q_mid = Q-dt/dx*(E-E_pre);
    return Q_mid
}

double get_next_q(double Q, double Q_mid, double E_mid, double E_mid_next, double dt, double dx){
    double Q_next = 
}

void output(double a){
	cout << "答えは" << a << "です。" << endl;
}