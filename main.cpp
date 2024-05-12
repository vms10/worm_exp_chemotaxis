// ***************************************
// Worm Chemotaxis
// ***************************************

#include "WormAgent.h"
#include <sstream>
//#define chem_index
//#define trajectories
#define chem_index_robust

// Global constants

// Global variables
int	VectSize = 25;

// ------------------------------------
// Behavioral Analysis
// ------------------------------------
void ChemIndex()
{
	ofstream ExampleFile;
	double t;
	long IDUM=time(0);
	int repetitions; 
	RandomState rs;
	rs.SetRandomSeed(IDUM);
	ExampleFile.open("/Users/Maria/Desktop/experience_dependent_c_elegans/CI_breed100.txt");
	TVector<double> phenotype;
	phenotype.SetBounds(1, VectSize);
	phenotype(1) = 825; //alpha
	phenotype(2) = 300.0; // K_hill
    phenotype(3) = 50.0; // delta
    phenotype(4) = 0.12; //delta_PKG
    phenotype(5) = 0.12; // gamma
    phenotype(6) = 1.0; // beta
    phenotype(7) = 1.0; // delta_ca
    phenotype(8) = 2.0; // b
    phenotype(9) = 0.7; // beta_DAG
    phenotype(10) = 0.001; // delta_DAG
    phenotype(11) = 1.39978563772776 - 0.05466237942122176; // alfa_glu
    phenotype(12) = 1.0; // alfa_delta
    phenotype(13) = 0.0; // theta
    phenotype(14) = 0.1; // tau
    phenotype(15) = 0.2; // inhWeight
    phenotype(16) = 92.0; // b_inh
    phenotype(17) = -5.0; // theta_inh
    phenotype(18) = 1.0; // excWeight
    phenotype(19) = 27.0; // b_exc
    phenotype(20) = -40.0; // theta_exc
    phenotype(21) = 0.03; // omega_low
    phenotype(22) = 50.3; // omega_high
    phenotype(23) = 0.0993; // V_low
	phenotype(24) = 0.05466237942122176; // alpha_glu_low
	phenotype(25) = 100; // C_breed
    
	WormAgent Worm(VectSize);
	Worm.SetParameters(phenotype);
	for (repetitions = 1; repetitions <= 600; repetitions++)
	{
		Worm.ResetAgentsBody(rs); //pos inicial (0,0) y condiciones iniciales
		#ifdef GRAD_STEP
			Worm.InitialiseSensorHistory();
		#endif
		//Worm.UpdateChemCon(rs);
		//t=0;
		for (t = StepSize; t < Preexposure; t += StepSize)
		{
			Worm.UpdateChemConstante(rs);
			Worm.UpdateSensors(StepSize);
			Worm.Step(StepSize, rs, t);
		}
		Worm.SetPositionX(0.0);
		Worm.SetPositionY(0.0);
		
		//Worm.PrintDetail(ExampleFile, t);

		for (t = Preexposure; t <= (Preexposure + RunDuration); t += StepSize)
		{
			Worm.UpdateChemCon(rs);
			Worm.UpdateSensors(StepSize);
			Worm.Step(StepSize, rs,t);
		}
		Worm.PrintDetail(ExampleFile, t);
		cout << repetitions << " " << endl;
	}

	ExampleFile.close();
}

void ExampleRun()
{
	ofstream ExampleFile;
	double t;
	long IDUM=time(0);
	int repetitions;
	RandomState rs;
	rs.SetRandomSeed(IDUM);
	ExampleFile.open("/Users/Maria/Desktop/experience_dependent_c_elegans/example.txt");
	TVector<double> phenotype;
	phenotype.SetBounds(1, VectSize);
	phenotype(1) = 825; //alpha
	phenotype(2) = 300.0; // K_hill
	phenotype(3) = 50.0; // delta
	phenotype(4) = 0.12; //delta_PKG
	phenotype(5) = 0.12; // gamma
	phenotype(6) = 1.0; // beta
	phenotype(7) = 1.0; // delta_ca
	phenotype(8) = 2.0; // b 
	phenotype(9) = 0.7; // beta_DAG
	phenotype(10) = 0.001; // delta_DAG
	phenotype(11) = 1.39978563772776 - 0.05466237942122176; // alfa_glu
	phenotype(12) = 1.0; // alfa_delta
	phenotype(13) = 0.0; // theta
	phenotype(14) = 0.1; // tau
	phenotype(15) = 0.2; // inhWeight
	phenotype(16) = 92.0; // b_inh
	phenotype(17) = -5.0; // theta_inh
	phenotype(18) = 1.0; // excWeight
	phenotype(19) = 27.0; // b_exc
	phenotype(20) = -40.0; // theta_exc
	phenotype(21) = 0.03; // omega_low
	phenotype(22) = 50.3; // omega_high
	phenotype(23) = 0.0993; // V_low
	phenotype(24) = 0.05466237942122176; // alpha_glu_low
	phenotype(25) = 50; // C_breed
	
	WormAgent Worm(VectSize);
	Worm.SetParameters(phenotype);
	for (repetitions = 1; repetitions <= 10; repetitions++)
	{
		Worm.ResetAgentsBody(rs); //pos inicial (0,0) y condiciones iniciales
		#ifdef GRAD_STEP
			Worm.InitialiseSensorHistory();
		#endif
		for (t = StepSize; t < Preexposure; t += StepSize)
		{
			Worm.UpdateChemConstante(rs);
			Worm.UpdateSensors(StepSize);
			Worm.Step(StepSize, rs, t);
		}
		Worm.SetPositionX(0.0);
		Worm.SetPositionY(0.0);
		
		Worm.PrintDetail(ExampleFile, t);

		for (t = Preexposure; t <= (Preexposure + RunDuration); t += StepSize)
		{
			Worm.UpdateChemCon(rs);
			Worm.UpdateSensors(StepSize);
			Worm.Step(StepSize, rs,t);
			Worm.PrintDetail(ExampleFile, t);
		}
	}
	ExampleFile.close();
}


// ------------------------------------
// robustness analysis
// ------------------------------------
void ChemIndex_robustness()
{
	//ofstream ExampleFile;
	double t;
	long IDUM=time(0);
	int repetitions; 
	RandomState rs;
	rs.SetRandomSeed(IDUM);
	//ExampleFile.open("/Users/Maria/Desktop/experience_dependent_c_elegans/CI_breed100.txt");
	TVector<double> phenotype;
	phenotype.SetBounds(1, VectSize);
	phenotype(1) = 825; //alpha
	phenotype(2) = 300.0; // K_hill
    phenotype(3) = 50.0; // delta
    phenotype(4) = 0.12; //delta_PKG
    phenotype(5) = 0.12; // gamma
    phenotype(6) = 1.0; // beta
    phenotype(7) = 1.0; // delta_ca
    phenotype(8) = 2.0; // b
    phenotype(9) = 0.7; // beta_DAG
    phenotype(10) = 0.001; // delta_DAG
    phenotype(11) = 1.39978563772776 - 0.05466237942122176; // alfa_glu
    phenotype(12) = 1.0; // alfa_delta
    phenotype(13) = 0.0; // theta
    phenotype(14) = 0.1; // tau
    phenotype(15) = 0.2; // inhWeight
    phenotype(16) = 92.0; // b_inh
    phenotype(17) = -5.0; // theta_inh
    phenotype(18) = 1.0; // excWeight
    phenotype(19) = 27.0; // b_exc
    phenotype(20) = -40.0; // theta_exc
    phenotype(21) = 0.03; // omega_low
    phenotype(22) = 50.3; // omega_high
    phenotype(23) = 0.0993; // V_low
	phenotype(24) = 0.05466237942122176; // alpha_glu_low
	phenotype(25) = 50; // C_breed
    
	int parameter;
	for (parameter = 1; parameter < VectSize; parameter++)
	//for (parameter = 1; parameter < 3; parameter++)
	{
		//int chequeo = 1;
		//if (parameter >1)
		//{
		//	chequeo = parameter-1;
		//}
		//cout << phenotype(chequeo) << endl; 
		double old_value = phenotype(parameter);
		phenotype(parameter) = phenotype(parameter)*0.8;
		//phenotype(parameter) = phenotype(parameter)*1;

		double CI_HIGH, CI_LOW, CI_START;
		CI_HIGH = 0.0;
		CI_LOW = 0.0; 
		CI_START = 0.0;

		WormAgent Worm(VectSize);
		Worm.SetParameters(phenotype);
		TVector<double> cis;
		cis.SetBounds(1, 6);
		int count_cis = 1;
		for (repetitions = 1; repetitions <= 600; repetitions++)
		{
			Worm.ResetAgentsBody(rs); //pos inicial (0,0) y condiciones iniciales
			#ifdef GRAD_STEP
				Worm.InitialiseSensorHistory();
			#endif
			//Worm.UpdateChemCon(rs);
			//t=0;
			for (t = StepSize; t < Preexposure; t += StepSize)
			{
				Worm.UpdateChemConstante(rs);
				Worm.UpdateSensors(StepSize);
				Worm.Step(StepSize, rs, t);
			}
			Worm.SetPositionX(0.0);
			Worm.SetPositionY(0.0);
			
			//Worm.PrintDetail(ExampleFile, t);

			for (t = Preexposure; t <= (Preexposure + RunDuration); t += StepSize)
			{
				Worm.UpdateChemCon(rs);
				Worm.UpdateSensors(StepSize);
				Worm.Step(StepSize, rs,t);
			}
			
			
			//Worm.PrintDetail(ExampleFile, t);
			//cout << repetitions << " " << endl;

			if ( (pow(Worm.PositionX(),2) + pow(Worm.PositionY(),2)) < 1.0){
				CI_START+=1;
			}
			
			else if ( ( pow((Worm.PositionX()-3.0),2) + pow(Worm.PositionY(),2) ) < pow(1.05,2) ){
				CI_HIGH+=1;
			}
			else if ( ( pow((Worm.PositionX()+3.0),2) + pow(Worm.PositionY(),2) ) < pow(1.05,2) ) {
				CI_LOW+=1;
			}
			
			if (repetitions%100 == 0){
				double chem_index_final;
				if (CI_START<100.0)
				{
					chem_index_final = (CI_HIGH-CI_LOW)/(100.0-CI_START);
				}
				else
				{
					chem_index_final = 0;
				}
				
				cis(count_cis) = chem_index_final;
				count_cis = count_cis +1;
				CI_HIGH = 0.0;
				CI_LOW = 0.0; 
				CI_START = 0.0;
			}

		}
		cout << old_value << " " << phenotype(parameter)  << " ";
		phenotype(parameter) = old_value; 
		//cout << phenotype(parameter) << endl; 
		for (int idx=1; idx<=6;idx++)
		{
			cout << cis(idx) << " ";
		}
		cout << endl;
	}
		//ExampleFile.close();
}



// ------------------------------------
// The main program
// ------------------------------------
#ifdef chem_index_robust
int main (int argc, const char* argv[])
{
	ChemIndex_robustness();
	return 0;
}
#endif

#ifdef chem_index
int main (int argc, const char* argv[])
{
	ChemIndex();
	return 0;
}
#endif

#ifdef trajectories
int main (int argc, const char* argv[])
{
	ExampleRun();
	return 0;
}
#endif
