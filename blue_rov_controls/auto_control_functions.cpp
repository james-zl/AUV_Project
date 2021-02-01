
#include "auto_control_functions.h"
#include <iostream> 
#include <string>
#include <cmath>

using namespace std; 


// auto_control_functions::auto_control_functions(){
	
// }



auto_control_functions::AFC_outputs auto_control_functions::adaptive_feedback_control(double uOpt[6], double xEst[6], double q, double antiWindup[6], 
		double theta0[16], double eta_r[6], double deta_r[6], double d2eta_r[6], double uU, double uL, double beta[6], double dt){
	
	//AFC_outputs outputs; 
	psi = xEst[5]; 
	theta = xEst[4]; 
	phi = xEst[3]; 

//Reference Model*******************************************************************************

 	J[0][0] = cos(psi)*cos(theta); 
 	J[0][1] = cos(psi)*sin(phi)*sin(theta) - cos(phi)*sin(psi); 
 	J[0][2] = sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta); 
 	J[1][0] = cos(theta)*sin(psi); 
 	J[1][1] = cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta); 
 	J[1][2] = cos(phi)*sin(psi)*sin(theta) - cos(psi)*sin(phi); 
 	J[2][0] = -sin(theta); 
 	J[2][1] = cos(theta)*sin(phi); 
 	J[2][2] = cos(phi)*cos(theta); 
 	J[3][3] = 1; 
 	J[3][4] = sin(phi)*tan(theta); 
 	J[3][5] = cos(phi)*tan(theta); 
 	J[4][4] = cos(phi); 
 	J[4][5] = -sin(phi); 
 	J[5][4] = sin(phi)/cos(theta);  
 	J[5][5] = cos(phi)/cos(theta); 
 	J[0][3] = J[0][4] = J[0][5] = J[1][3] = J[1][4] = J[1][5] = J[2][3] = J[2][4] = J[2][5] = J[3][0] 
 	= J[3][1] = J[3][2] = J[4][0] = J[4][1] = J[4][2] = J[4][3] = J[5][0] = J[5][1] = J[5][2] = J[5][3] = 0; 

 	for(int i = 0; i < 6; i++)
 		eta[i] = xEst[i]; 
 	for(int i = 0; i < 6; i++)
 		nu[i] = xEst[i+6];  

 	temp = 0.0f; 
 	for(int i = 0; i < 6; i++){
 		for(int j = 0; j < 6; j++){
 			temp = temp + J[i][j] * nu[j]; 
 		}
 		deta[i] = temp; 
 		temp = 0.0f; 	
 	}


 	//to be implemented 
 	//[A, B, ~] = referenceModel(beta, 6, dt);   

 	//y0 = [eta_r; deta_r; d2eta_r];
 	for(int i = 0; i < 6; i++)
 		y0[i] = eta_r[i]; 
 	for(int i = 0; i < 6; i++)
 		y0[i+6] = deta_r[i]; 
 	for(int i = 0; i < 6; i++)
 		y0[i+12] = d2eta_r[i]; 
 	
 	//y = A*y0 + B*uOpt 
 	//A, B are missing, need to implement referenceModel 
 	temp = 0.0f; 
 	for(int i = 0; i < 18; i++){
 		for(int j = 0; j < 18; j++){
 			temp = temp + A[i][j] * y0[j]; 
 		}
 		y[i] = temp; 
 		temp = 0.0f; 	
 	}
 	temp = 0.0f; 
	for(int i = 0; i < 18; i++){
 		for(int j = 0; j < 6; j++){
 			temp = temp + B[i][j] * uOpt[j]; 
 		}
 		y[i] = y[i] + temp; 
 		temp = 0.0f; 	
 	}

 // eta_r = y(1:6); deta_r = y(7:12); d2eta_r = y(13:18);
	for(int i = 0; i < 6; i++)
		eta_r[i] = y[i]; 
	for(int i = 0; i < 6; i++)
		deta_r[i] = y[i+6]; 
	for(int i = 0; i < 6; i++)
		d2eta_r[i] = y[i+12]; 
//**********************************************************************************

//anti windup **********************************************************************
	for(int i = 0; i < 6; i++)
		Iawp[i] = 1; 
	for(int i = 0; i < 6; i++){
		if (antiWindup[i] >= uU and (eta[i] - eta_r[i]) <= 0)
			Iawp[i] = 0; 
		else if (antiWindup[i] <= uL and (eta[i] - eta_r[i]) >= 0)
			Iawp[i] = 0; 
	}
//**********************************************************************************

//PID ******************************************************************************
	//q = q + Iawp.*(eta-eta_r)*dt;
	for(int i = 0; i < 6; i++)
		temp_array1[i] = 0.0; 
	for(int i = 0; i < 6; i++) 
		temp_array1[i] = temp_array1[i] + Iawp[i]*(eta[i] - eta_r[i])*dt; 
	for(int i = 0; i < 6; i++) 
		q[i] = q[i] + temp_array1[i]; 

	double lambda[6] = {0.4, 0.4, 0.7, 4.0, 4.0, 4.0}; 
//a_eta = d2eta_r - 3*lambda.*(deta-deta_r) - 3*lambda.^2.*(eta-eta_r) - lambda.^3.*q; 
//not entirely sure if the following computes the equation above correctly 
	for(int i = 0; i < 6; i++){
		temp_array1[i] = 0.0; 
		temp_array2[i] = 0.0; 
		temp_array3[i] = 0.0; 
	}

	for(int i = 0; i < 6; i++)
		temp_array1[i] = temp_array1[i] + 3*lambda[i]*(deta[i]-deta_r[i]); 
	for(int i = 0; i < 6; i++)
		temp_array2[i] = temp_array2[i] + 3*(lambda[i]*lambda[i])*(eta[i] - eta_r[i]); 
	for(int i = 0; i < 6; i++)
		temp_array3[i] = temp_array3[i] + (lambda[i]*lambda[i]*lambda[i])*q[i]; 
	for(int i = 0; i < 6; i++)
		a_eta[i] = d2eta_r[i] - temp_array1[i] - temp_array2[i] - temp_array3[i]; 

//*********************************************************************************

//Change of reference frame *******************************************************
//compute a_nu 
	psi = eta[5]; 
	theta = eta[4]; 
	phi = eta[3]; 
	dpsi = nu[5]; 
	dtheta = nu[4]; 
	dphi = nu[3]; 

	dJ[0][0] = - dpsi*cos(theta)*sin(psi) - dtheta*cos(psi)*sin(theta); 
	dJ[0][1] = dphi*sin(phi)*sin(psi) - dpsi*cos(phi)*cos(psi) + dphi*cos(phi)*cos(psi)*sin(theta) 
	+ dtheta*cos(psi)*cos(theta)*sin(phi) - dpsi*sin(phi)*sin(psi)*sin(theta); 
	dJ[0][2] = dphi*cos(phi)*sin(psi) + dtheta*cos(phi)*cos(psi)*cos(theta) - 
	dphi*cos(psi)*sin(phi)*sin(theta) - dpsi*cos(phi)*sin(psi)*sin(theta); 
	dJ[0][3] = dJ[0][4] = dJ[0][5] = 0; 
	dJ[1][0] = dpsi*cos(psi)*cos(theta) - dtheta*sin(psi)*sin(theta); 
	dJ[1][1] = dphi*cos(phi)*sin(psi)*sin(theta) - dpsi*cos(phi)*sin(psi) 
	- dphi*cos(psi)*sin(phi) + dpsi*cos(psi)*sin(phi)*sin(theta) + dtheta*cos(theta)*sin(phi)*sin(psi); 
	dJ[1][2] = dpsi*sin(phi)*sin(psi) - dphi*cos(phi)*cos(psi) + dpsi*cos(phi)*cos(psi)*sin(theta) + 
	dtheta*cos(phi)*cos(theta)*sin(psi) - dphi*sin(phi)*sin(psi)*sin(theta); 
	dJ[1][3] = dJ[1][4] = dJ[1][5] = 0; 
	dJ[2][0] = -dtheta*cos(theta); 
	dJ[2][1] = dphi*cos(phi)*cos(theta) - dtheta*sin(phi)*sin(theta); 
	dJ[2][2] = - dphi*cos(theta)*sin(phi) - dtheta*cos(phi)*sin(theta); 
	dJ[2][3] = dJ[2][4] = dJ[2][5] = 0; 
	dJ[3][0] = dJ[3][1] = dJ[3][2] = dJ[3][3] = 0; 
	dJ[3][4] = dphi*cos(phi)*tan(theta) + (dtheta*sin(phi))/cos(theta)^2; 
	dJ[3][5] = (dtheta*cos(phi))/cos(theta)^2 - dphi*sin(phi)*tan(theta); 
	dJ[4][0] = dJ[4][1] = dJ[4][2] = dJ[4][3] = 0; 
	dJ[4][4] = -dphi*sin(phi); 
	dJ[4][5] = -dphi*cos(phi); 
	dJ[5][0] = dJ[5][1] = dJ[5][2] = dJ[5][3] = 0; 
	dJ[5][4] = (dphi*cos(phi)*cos(theta) + dtheta*sin(phi)*sin(theta))/cos(theta)^2; 
	dJ[5][5] = -(dphi*cos(theta)*sin(phi) - dtheta*cos(phi)*sin(theta))/cos(theta)^2; 


	//to be implemented. Use alglib to invert the matrix J 
	//a_nu = inv(J)*(a_eta-dJ*nu);

//*********************************************************************************


//Pramters update law *************************************************************
	double c1[6] = {3, 3, 1, 1, 1, 1}; 
	






}



auto_control_functions::double auto_control_functions::reference_model(double beta[6], int Nu1, double dt){

}



// auto_control_functions::~auto_control_functions(){
// }