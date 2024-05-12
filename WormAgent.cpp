/*
 *  WormAgent.cpp
 */

#include "WormAgent.h"
#include "random.h"


// ****************************
// Constructors and Destructors
// ****************************

// The constructor
WormAgent::WormAgent(int entero)
{
	InitialiseSensors();
}

// The destructor
WormAgent::~WormAgent()
{
	InitialiseSensors();
}

// *********
// Setting paremeters
// *********

void WormAgent::SetParameters(TVector<double> &v)
{	
	alpha = v(1);
	K_hill = v(2);
	delta = v(3);
	delta_PKG = v(4);
	gamma = v(5);
	beta = v(6);
	delta_ca = v(7);
	b = v(8); 
	beta_DAG = v(9);
	delta_DAG = v(10);
	alpha_glu = v(11);
	alpha_delta = v(12);
	theta = v(13);
	tau = v(14); 
	inhWeight = v(15);
	b_inh = v(16);
	theta_inh = v(17);
	excWeight = v(18);
	b_exc = v(19);
	theta_exc = v(20);
	omega_low = v(21);
	omega_high = v(22);
	V_low = v(23);
	alpha_glu_low = v(24);
	C_breed = v(25);
}


// *******
// Initialising
// *******

void WormAgent::InitialiseSensors(){
	cGMP = 0.0;
	PKG = 0.0;
	Ca = 0.0;
	DAG = 0.0;
	AIB = 0.0;
	chemConHistory.SetBounds(1, len_chem);
}





void WormAgent::InitialiseSensorHistory()
{
	timer = 0;
	for (int i = timer; i < (int)(0.5*RunDuration/StepSize); i++)
		chemConHistory[i] = C_breed;
	for (int i =  (int)(0.5*RunDuration/StepSize); i <  (int)(RunDuration/StepSize); i++)
		chemConHistory[i] = C_back;
}

// *******
// Resetting
// *******

void WormAgent::ResetAgentsBody(RandomState &rs)
{

double tempangle = rs.UniformRandom(0,2*Pi);

#ifdef GRAD_CONE				
	distanceToCentre = 0.0;	
#endif
#ifdef GRAD_INVCONE				
	distanceToCentre = 0.0;
#endif
#ifdef GRAD_GAUS				
	distanceToCentre = 0.0;	
#endif		
#ifdef GRAD_STEP
	distanceToCentre = 0.0;
#endif
	px = cos(tempangle)*distanceToCentre;
	py = sin(tempangle)*distanceToCentre;
	vx = 0.0;
	vy = 0.0;
	
	cGMP = 0.0;
	PKG = 0.0;
	Ca = 0.0;
	DAG = 0.0;
	AIB = 0.0;

	orient = rs.UniformRandom(0,2*Pi); 
}



// *******
// Updating
// *******

void WormAgent::UpdateChemCon( RandomState &rs )
{
	distanceToPeak = sqrt( pow((px-pos_x_peak),2) + pow((py-pos_y_peak),2) );
	distanceToMin = sqrt( pow((px-pos_x_peak_min),2) + pow((py-pos_y_peak),2) );
#ifdef GRAD_CONE				
	chemCon = -distanceToPeak*gradSteepness; 
#endif
#ifdef GRAD_INVCONE				
	chemCon = distanceToPeak*gradSteepness; 
#endif
#ifdef GRAD_GAUS				
	chemCon = C_max*exp(-(pow(distanceToPeak,2))/ChemDiffConst) +C_back -C_min*exp(-(pow(distanceToMin,2))/ChemDiffConst);
#endif	
#ifdef GRAD_ISO				
	chemCon = - (fabs((distanceToPeak - MaxDist)*gradSteepness)); 
#endif	
#ifdef GRAD_STEP				
	chemCon = chemConHistory[timer]; 
	timer++;
#endif
	
}


void WormAgent::UpdateChemConstante( RandomState &rs )
{			
	chemCon = C_breed; 
}



void WormAgent::UpdateSensors(double stepsize)
{
	der_cGMP = alpha*hill_repressor(chemCon, K_hill) - delta*cGMP; 
	der_PKG = gamma*cGMP - delta_PKG*PKG;
	I_se = cGMP - PKG;
 	der_Ca = beta*tanh(I_se, b) -  delta_ca*Ca;
	der_DAG = beta_DAG*Ca - delta_DAG*DAG;
	if (DAG>=theta){
		glu_basal = alpha_glu_low + alpha_glu;
	}
	else{
		glu_basal = alpha_glu_low;
	}
	der_AIB = (inhWeight*(-sig(b_inh*(glu_basal + alpha_delta*Ca) + theta_inh) + 1) + excWeight*sig(b_exc*(glu_basal + alpha_delta*Ca) + theta_exc) - AIB)/tau;

	cGMP += stepsize * der_cGMP;
	PKG += stepsize * der_PKG;
	Ca += stepsize * der_Ca;
	DAG += stepsize * der_DAG;
	AIB += stepsize * der_AIB;
}



void WormAgent::PrintDetail( ofstream &file, double timestep)
{
	file << timestep-Preexposure << " ";
	file << chemCon << " ";
	file << cGMP << " ";
	file << PKG << " ";
	file << Ca << " ";
	file << DAG << " ";
	file << (AIB-1.1)/0.02 << " ";
	file << glu_basal << " ";
	file << proba_final << " ";
	file << px << " ";
	file << py << " ";
	file << endl;
}

// Step the bacteria simulation
void WormAgent::Step(double StepSize, RandomState &rs, double timestep)
{	
	
	if (AIB<=V_low){
		proba_completa = StepSize*omega_low;
	}
	
	else {
		proba_completa = StepSize*omega_high;
	}
	
	proba_final = proba_completa;
	
	if (rs.UniformRandom(0.0, 1.0) < proba_final){ 
		SetOrientation(rs.UniformRandom(0, 2*Pi));
	}		
	
	
	// Update the velocity 
	vx = cos(orient) * MaxVel;	//MaxVel
	vy = sin(orient) * MaxVel;	//MaxVel
	
	// Move the agent
	px += StepSize * vx;
	py += StepSize * vy;

	if ( (pow(px,2) + pow(py,2)) > pow(MaxDist,2)){  
		while ((pow(px,2) + pow(py,2)) > pow(MaxDist,2)) {
        SetOrientation(rs.UniformRandom(0, 2*Pi));;
        vx = cos(orient) * MaxVel;	//MaxVel
		vy = sin(orient) * MaxVel;	//MaxVel
		
		// Move the agent
		px += StepSize * vx;
		py += StepSize * vy;
    }
	}

}

