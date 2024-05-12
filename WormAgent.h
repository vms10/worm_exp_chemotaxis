/*
 *  WormAgent.h
 *
 */


#include "VectorMatrix.h"
#include "random.h"
#include <iostream>
#include <math.h>

#pragma once

#define GRAD_GAUS				
//#define GRAD_STEP	

// Global constants
const double	StepSize		=	0.01;			// Fastest time-constant is now 0.1 (10ms)
const double	Pi				=	3.1415926;
const double	MaxDist			=	4.25; // Half the radius of the big petri dish (in cm)
const double	MaxVel			=	0.022;			// Forward velocity (in cm/s) CHECK WHAT THE REAL VELOCITY SHOULD BE
const double	ChemDiffConst	=	2*pow(0.7,2);	// Simulated chemical environment according to Ward, 1973 as described in Ferree and Lockery 1999 equation 14.
const double	pos_x_peak		= 3.0;
const double	pos_y_peak		= 0.0;
const double	pos_x_peak_min		= -3.0;
const double	C_max			=	45.0; 
const double	C_back			=	50.0; 
const double	C_min			=	20.0; 
const double	RunDuration		= 600.1; 
const double	Preexposure		= 10000.0; // 5000 when #define GRAD_STEP	
const int   	len_chem    	= 	(int) (RunDuration/StepSize);

inline double tanh(double x, double b) 
{
  return (exp(b*x)-exp(-b*x))/(exp(b*x)+exp(-b*x)); 
}

inline double sig(double x) 
{
  return 1/(1 + exp(-x)); 
}

inline double hill_repressor(double x, double K_hill) 
{
  return 1/(1 + x/K_hill); 
}


// The WormAgent class declaration
class WormAgent {
public:
	// The constructor
	WormAgent(int entero);	// Construct from evolutionary algorithm
	
	// The destructor
	~WormAgent();
	
	// Accessors
	
	double PositionX(void) {return px;};
	void SetPositionX(double newx) {px = newx;};
	double PositionY(void) {return py;};
	void SetPositionY(double newy) {py = newy;};
	double VelocityX(void) {return vx;};
	void SetVelocityX(double newx) {vx = newx;};
	double VelocityY(void) {return vy;};
	void SetVelocityY(double newy) {vy = newy;};
	double Orientation(void) {return orient;};
	void SetOrientation(double newo) {orient = newo;};
	double ChemCon(void) {return chemCon;};
	void SetChemCon(double newc) {chemCon = newc;};
	double DistanceToCentre(void) {return distanceToCentre;};
	
	
	// Control
	void SetParameters(TVector<double> &v);
	void InitialiseSensors();
	void InitialiseSensorHistory();
	void ResetAgentsBody(RandomState &rs);
	void UpdateChemCon( RandomState &rs );
	void UpdateChemConstante( RandomState &rs );
	void UpdateSensors(double stepsize);
	void PrintDetail(ofstream &file, double timestep);
	void Step(double StepSize, RandomState &rs, double timestep);

	double px, py, vx, vy, orient, distanceToCentre;
	double chemCon;
	int timer;
	//TVector<double> onWeight,offWeight,cpgWeight;
	//double chemConHistory[len_chemi];	// RunDuration / StepSize
	TVector<double> chemConHistory;
	double proba_final, proba_completa;
	double delta, alpha, delta_PKG, gamma, beta, delta_ca, b, alpha_glu, alpha_glu_low, alpha_delta, theta,b_inh,b_exc,omega_low,omega_high,V_low, tau, excWeight, omega_0, inhWeight, theta_inh, theta_exc, glu_basal, theta_out, b_out, K_hill;
	double I_se, distanceToPeak, distanceToMin;
	double beta_DAG, delta_DAG, C_breed;
	double cGMP , PKG , Ca , DAG , AIB , der_cGMP, der_PKG, der_Ca, der_DAG, der_AIB;
};

