#ifndef AUTO_CONTROL_FUNCTIONS_H
#define AUTO_CONTROL_FUNCTIONS_H

#include <string> 
using namespace std; 


class auto_control_functions{ //not sure if I need to make it a class 

private:

//parameters for adaptive feedback control 
struct AFC_outputs {
	double u[6]; 
	double q_new; 
	double antiWindup_new[6]; 
	double theta0_new[16]; 
	double eta_r_new[6]; 
	double deta_r[6]; 
	double d2eta_r[6]; 
}; 

AFC_outputs afc_outputs; 

//output parameters for reference model 
struct RM_outputs {
	double A[18][18]; 
	double B[18][6]; 
}; 

RM_outputs rm_outputs; 

double temp_array1[6]; 
double temp_array2[6]; 
double temp_array3[6]; 
double temp;
double psi; 
double theta; 
double phi; 
double dpsi; 
double dtheta; 
double dphi; 
double J[6][6]; 
double dJ[6][6]; 
double eta[6]; 
double nu[6]; 
double deta[6]; 
double Ac[18][18]; 
double Bc[18][6]; 
double A[18][18]; 
double B[18][6]; 
double beta[6]; 
double y0[18]; 
double y[18]; 
double Iawp[6]; 
double a_eta[6]; 
double a_nu[6]; 



public: 

	AFC_outputs adaptive_feedback_control(double uOpt[6], double xEst[6], double q, double antiWindup[6], double theta0[16], 
		double eta_r[6], double deta_r[6], double d2eta_r[6], double uU, double uL, double beta[6], double dt);

	double reference_model(double beta[6], int Nu1, double dt); 

};
#endif 



