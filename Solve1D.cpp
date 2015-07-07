#include "ConsState.h"
#include "SolidSystem.h"
#include <libconfig.h++>
#include <iomanip>
#include <cstdlib>

using namespace std;
using namespace Eigen;
using namespace libconfig;

/* #include <config5cpp/Configuration.h> */

/* using namespace config4cpp; */
/* Global variables */ 
/* Gloabl functions */
void initialize(ConsState*);

int e_size = 13;

double internalEnergy(Matrix3d, double);

void getState(Vector3d &u, Matrix3d &F, double &rho, double& S, Setting&);

string output = "", outputExact = "", baseDir;
/* Solution vector */
//using the Vector3d class for most operations
//Slic variables
int gCs = 2;

double C, disctyX, aNew, tUnits, dX, dt,w;

int N, cells, slLimChoice;

string test;
/* Slic functions */
void setBCs(ConsState*);

void printArray(ConsState*);

double slope_limiter(double);

double ksi_r(double);

ConsState delta_Bar(ConsState*, int);

void solve(ConsState*);

void calcTimeStep(ConsState*);

void slicSolve(ConsState*, ConsState*, ConsState*);

ConsState calcForceFlux(ConsState*, ConsState*, int);

void outputToFile(ConsState*, string);

enum slope_lim{superbee, minbee, van_leer, albalda};//0,1,2,3

enum ICType { Barton1, Barton2, Barton3, Barton4 };

ICType icType;

struct DomainPs {
	int Ni;
	int GNi;
	int GC;
	int starti;
	int endi;
	double Lx;
	double dx;

	DomainPs(const int a_Ni, const int a_GC,
			const double a_Lx) :
		Ni(a_Ni), GNi(a_Ni+2*a_GC), GC(a_GC),
		starti(a_GC), endi(a_Ni+a_GC),
		Lx(a_Lx), dx(a_Lx/a_Ni) {}
	DomainPs() {}
};

int main(int argc, char ** argv)
{
	/* Config cfg; */
	/* Matrix3d F_t; */
	/* Vector3d u_t; */
	/* double S_t, rho_t, e_t; */
	/* try */
	/* { */
	/* 	cfg.readFile("./Configurations/testCases.cfg"); */
	/* } */
	/* catch(const FileIOException &fioex) */
	/* { */
	/* 	std::cerr << "I/O error while reading file." << std::endl; */
	/* 	return(EXIT_FAILURE); */
	/* } */
	/* catch(const ParseException &pex) */
	/* { */
	/* 	std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine() */
	/* 		<< " - " << pex.getError() << std::endl; */
	/* 	return(EXIT_FAILURE); */
	/* } */
	/* try */
	/* { */
	/* 	cfg.lookupValue("slChoice", slLimChoice); */
	/* 	cfg.lookupValue("base", baseDir); */
	/* 	Setting &settingState = cfg.lookup("testCases.testCase1.left"); */
	/* 	cout << "Length " << settingState.getLength() << endl; */
	/* 	int temp = settingState["rho"].getLength(); */

	/* 	//get leftstate */
	/* 	getState(u_t, F_t, rho_t, S_t, settingState); */
	/* 	e_t = internalEnergy(F_t, S_t); */	
	/* 	ConsState ConsStateL(u_t, F_t, rho_t, S_t, e_t); */
	/* 	cout << ConsStateL << endl; */	
	/* 	settingState.lookupValue("rho", rho_t); */
	/* 	cout << rho_t << endl; */
	/* 	cout << "slChoice" <<slLimChoice << endl; */
	/* } */
	/* catch(const SettingNotFoundException &nfex) */
	/* { */
	/* 	cerr << "No 'name' setting in configuration file." << endl; */
	/* } */
	string output, outputInitial;
	cout << "Parameters " << argc << endl;
	if(argc != 5){
		cout << "Insufficient Parameters: exit" << endl;
		exit(1);
	}
	else {
		const string icStr = string(argv[1]);
		if(icStr == "barton1") 
		{
			icType = Barton1;
		}
		else if(icStr == "barton2")
		{
			icType = Barton2;
		}
		else 
		{
			cerr << "Unknown problem type \'" << icStr << "\' - options are barton1, barton2, barton3, barton4, barton5\n";
			exit(1);
		}
		cout << "Solve Test " << icStr << endl;
		cells = atoi(argv[2]);
		outputInitial = argv[3];
		output = (argv[4]);
	}
	/* * Setting& setting = cfg.lookup("tests.test1"); *1/ */
	double domainLength = 0.01;
	N = 2*gCs + cells;
	dX = domainLength/(double)cells;
	ConsState* domain = new ConsState[N];
	C = 0.55, disctyX = 0.5, tUnits = 6.0e-7, slLimChoice = 0, w = 1.;

	/* string base = "/lsc/zeushome/ma595/Dropbox/2013-2014/Code/Solid/Output/"; */
	/* string sOut = "out1D" + convertToStr(cells) + ".dat"; */
	/* string sOutInit = "out1D" + convertToStr(cells) + "_initial.dat"; */
	/* string outputInitial = base + sOutInit; */
	/* output = base + sOut; */ 
	/* solve(domain, output); */

	initialize(domain);
	/* printArray(domain); */
	outputToFile(domain,outputInitial);
	solve(domain);	
	outputToFile(domain, output);
	delete[] domain;
	/* ConsState::checkEOSConstants(); */
	return 0;
}		

void getState(Vector3d &u, Matrix3d &F, double &rho, double& S, Setting& setting)
{
	setting.lookupValue("rho", rho);
	setting.lookupValue("S", S);
	int max = setting["u"].getLength();
	double temp = 0; 
	for (unsigned int i = 0; i < max; i++)
	{
		u(i) = setting["u"][i];
	}
	int count = 0;
	for (unsigned int i = 0; i < max; i++)
	{
		for (unsigned int j = 0; j < max; j++)
		{
			count = i*max + j;
			F(i,j) = setting["F"][count];
		}
	}
}

void solve(ConsState* U)
{
	double timeElapsed = 0;
	int Nt = 0; // Number of steps
	/* initialize(U); */
	ConsState UL[N];
	ConsState UR[N];
	//cout << "TUnits " << tUnits << endl;	
	while(timeElapsed < tUnits)
	{
		//1. calculate time step: CFL and boundary conditions pg 495
		calcTimeStep(U);
		if(Nt < 10)	{dt /= 10;}
		timeElapsed += dt;		
		if(timeElapsed > tUnits)
		{				
			dt = (tUnits - (timeElapsed - dt));	
			//timeElapsed -= dt;
		}
		solve(U, UL, UR);
		++Nt;
		/* if(Nt > 1) */
		/* { */
		/* 	cout << "Debug break" << endl; */		
		/* 	break; */			
		/* } */
		cout.precision(3);
		cout << "Steps " << setw(6) << Nt << setw(6) << "Time " << setw(6) << timeElapsed << setw(15) << "Remaining " << setw(6) << tUnits-timeElapsed << endl;
	}
}

void calcTimeStep(ConsState* U)
{
	double Smax = 0, a, u_1;
	a = 10000;
	/* double c = ConsState::soundSpeed(); */
	for (int i = 0; i < N; i++)
	{
		u_1 = fabs(U[i].u_(0));				
		if(U[i].rho_() > 0)
		{
			if((a + u_1) > Smax)
			{
				Smax = a + u_1; //check
			}
		}		 
		else
		{
			cout << "density is zero or negative: exit" << endl;
			exit(1);	
		}
	}	
	dt = C*dX/Smax;
	if(dt == 0)
	{
		cout << "dt == 0: exit" << endl;
		exit(1);
	}	
	cout << "Smax " << Smax << " dt " << dt << endl;
}


/* //equation of state parameters */
/* //--sound speed */	
/* double K_0; */
/* //--shear wave speed */
/* double b_0, B_0; */
/* //--Heat capacity at constant volume */
/* double c_v; */
/* //--non-linear dependence of K_0 and T on mass density */
/* double alpha, beta, gamma; */ 
/* //--Temperature */
/* double T_0; */

void initialize(ConsState* domain)
{
	//Parameters for the equation of state
	//	ConsStateL

	double rho_0L = 8930, rho_0R = 8930, eL(0), eR(0);
	disctyX = 0.5;
	Matrix3d FL, FR;
	Vector3d uL, uR;
	double SL, SR;

	if(icType == Barton1)
	{
		FL <<	 0.98,  0,   0,
			 0.02,  1., 0.1, 
			 0,  0,   1.;
		uL <<    0, 500.,  1000.;		
		SL = 1000.;
		eL = internalEnergy(FL, SL);	

		FR <<    1.,  0,   0,
			 0,  1., 0.1, 
			 0,  0,   1.;
		uR <<    0,  0,   0;
		SR = 0.;
		eR = internalEnergy(FR, SR);	
	}
	if(icType == Barton2)
	{
		FL << 1., 0, 0,	
			 -0.01, 0.95, 0.02,
			 -0.015, 0, 0.9;
		uL << 2000., 0, 100.;
		SL = 0;
		eL = internalEnergy(FL, SL);

		FR << 1., 0, 0,
			 0.015, 0.95, 0,
			 -0.01, 0, 0.9;
		uR << 0, -30., -10.;
		SR = 0;
		eR = internalEnergy(FR, SR);
	}

	ConsState ConsStateL(uL, FL, rho_0L, SL, eL);
	ConsState ConsStateR(uR, FR, rho_0R, SR, eR);

	/* ConsState flux = ConsStateL.F(); */
	//Populate array 
	unsigned int i;
	for(i = gCs; i < cells*disctyX + gCs; i++)
	{
		domain[i] = ConsStateL;		
	}
	for(i; i < N - gCs; i++)
	{	
		domain[i] = ConsStateR;
	}	
	setBCs(domain);
	//print out initial domain
	// U vector 
	/* printArray(domain); */
}

inline double internalEnergy(Matrix3d F, double S)
{
	Vector3d I;
	double c_0, b_0, c_v, T_0, alpha, beta, gamma;
	c_0 = 4600, b_0 = 2100, c_v = 390, T_0 = 300, alpha = 1.001, beta = 3.001, gamma = 2.001;
	double B_0, K_0;
	B_0 = b_0*b_0;
	K_0 = c_0*c_0-(4./3.)*B_0;
	Matrix3d G = (F.inverse().transpose())*(F.inverse());	
	// Calculate invariants
	I(0) = G.trace();
	I(1) = 0.5*(pow(G.trace(), 2)-(G*G).trace());
	I(2) = G.determinant();

	double U, W;
	U = (K_0/(2.*alpha*alpha))*pow((pow(I(2),alpha/2.))-1,2.)+
		c_v*T_0*pow(I(2),gamma/2.)*(exp(S/c_v)-1.);
	W = (B_0/2.)*pow(I(2),beta/2.)*(pow(I(0),2)/3. - I(1));
	cout << "Internal SlicCode Start " << U + W << endl;
	return U + W;	
}

void printArray(ConsState* domain)
{
	/* inline std::ostream& operator<<(std::ostream& os, ConsState& param) */
	/* { */
	/* param.initialStates(os); */
	/* os << "Current states \n" << */
	/* 	  "Density\t"    << param.rho_() << "\n" */
	/* 	  "Energy\t"     << param.E_() << "\n" */
	/* 	  "U or F vector\t " << "\n"; */
	/* os.precision(3); */
	/* for (int i = 0; i < ConsState::e_size; i++) */
	/* { */
	/* 	os << setw(5) << " |" << setw(5) << param[i]; */
	/* } */	
	/* os << setw(5) << " | "; */
	/* return os; */
	/* } */
	cout.precision(2);
	cout << setw(7) << left << "0" << ' ';
	for(int i = 0; i < e_size; i++)
	{
		cout << setw(3) << " |" << setw(7) << ConsState::name[i];
	}
	cout << setw(7) <<  " |" << endl;
	for(int i = 0; i < 140; i++)
	{
		cout << '-';
	}
	cout << endl;
	for(int i = 0; i < N; i++)
	{
		cout << setw(7) << left << i << ' ' << domain[i] << endl;
	}
}

void setBCs(ConsState* domain)
{
	// Transmissive
	for (int i = 0; i < gCs; i++)
	{
		domain[i] = domain[2*gCs - i - 1];
		domain[cells + i + gCs] = domain[N - gCs - i - 1]; 
	}					
}	

// slope limiters pg 509 and 510
double slope_limiter(double r)
{
	//calculate ksi_r
	double ksi_sl = 0;
	//slope limiters:
	//cout << "Slope limiter " << slLimChoice << endl;
	switch(slLimChoice)
	{		
		case superbee:
			if(r <= 0)
				ksi_sl = 0;
			else if(r > 0 && r <= 0.5)
				ksi_sl = 2*r;
			else if(r >= 0.5 && r <= 1)
				ksi_sl = 1;
			else if (r >= 1)
				ksi_sl = min(r, min(ksi_r(r), 2.));
			break;
		case van_leer:		
			if(r <= 0)
				ksi_sl = 0;
			else
				ksi_sl = min(2*r/(1+r), ksi_r(r));
			break;
		case albalda:
			if(r <= 0)
				ksi_sl = 0;
			else 
				ksi_sl = min(r*(1+r)/(1+r*r), ksi_r(r));			
			break;
		case minbee:	
			//cout << "minbee"	<< endl;
			if(r <= 0)
				ksi_sl = 0;
			else if(r >= 0 && r <= 0.5)
				ksi_sl = r;
			else 
				ksi_sl = min(1., ksi_r(r));
			break;
	}
	return ksi_sl;			
}

double ksi_r(double r)
{
	/* double cNew = aNew*dt/dX; */
	//cout << "courant Number	" << cNew << endl;
	//return 2./(1. + r);
	double beta_fw = 1;//2/(0.1);
	return 2.*beta_fw/(1. - w + (1. + w)*r);
}

/**
 * Calculate the slope limiter and slope delta_i.
 * Update delta_i.
 * @param Ui Pointer to solution domain Ui
 * @param i Domain index of conserved vector  
 * @return updated slope vector delta_i
 */

// delta ibar slope on pg 509
ConsState delta_Bar(ConsState* Ui, int i)
{
	//1. calculate r
	ConsState num = Ui[i] - Ui[i-1];
	ConsState denom = Ui[i+1] - Ui[i];
	ConsState r;
	ConsState delta = 0.5*(1+w)*num;// + 0.5*(1-w)*denom;
	//aNew = Ui[i].soundSpeed();	
	//cout << i << endl;	
	for(int j = 0; j < 3; j++)
	{	
		if(num[j] == 0 && denom[j] == 0)
			r[j] = 1;
		else if(denom[j] == 0)
		{
			denom[j] = 1.e-10;
			r[j] = num[j]/denom[j];			
		}
		else 
			r[j] = num[j]/denom[j];		
	}
	//ConsState r = num/denom;
	ConsState delta_old = delta;	
	double ksi = 0;	
	for(int j = 0; j < e_size; j++)
	{
		ksi = slope_limiter(r[j]);
		//ksi = 1;
		delta[j] = ksi * delta[j];	 
	}
	return delta; 
} 

void solve(ConsState* U, ConsState* UL, ConsState* UR)
{
	int i;
	//cout << "Iteration " << endl;
	//N
	//printArray(U);
	ConsState slope_bar; 
	for(i = gCs - 1; i < N - 1; ++i) 
	{
		//2. calculate left and right extrapolated values (bar) pg 514 and slope pg 506.
		slope_bar = delta_Bar(U, i);
		//cout << "Slope " << slope_bar << endl;
		ConsState UiL = U[i] - 0.5 * slope_bar;
		ConsState UiR = U[i] + 0.5 * slope_bar;
		ConsState UbarAdd = 0.5 * (dt/dX) * (UiL.F() - UiR.F());
		UL[i] = UiL + UbarAdd;
		UR[i] = UiR + UbarAdd;
	}		
	UL[0] = U[0];
	UR[0] = U[0];//is required
	UL[N - 1] = U[N - 1];
	UR[N - 1] = U[N - 1];	
	//3. calculate force flux using LF and RI, and calculate new cell averaged Ui pg 494
	for(int i = 1; i < cells + gCs + 1; i++)
	{
		U[i] = U[i] + (dt/dX)*(calcForceFlux(UL, UR, i - 1) - calcForceFlux(UL, UR, i));
	}
	//cout << "solution " << endl;
	//printArray(U);
	setBCs(U);
}

//Force flux calculation pg 512
ConsState calcForceFlux(ConsState* UL, ConsState* UR, int i)
{
	ConsState F_LF = 0.5*(UR[i].F() + UL[i+1].F()) + 0.5*(dX/dt)*(UR[i] - UL[i+1]);
	ConsState U_RI = 0.5*(UR[i] + UL[i+1]) + 0.5*(dt/dX)*(UR[i].F() - UL[i+1].F());
	ConsState FORCE = 0.5*(F_LF + U_RI.F());
	return FORCE;
}

void outputToFile( ConsState* U, string fileName)
{
	/* int ret = system(("mkdir -p " + fileName).c_str()); */
	/* if(ret!=0) cerr << "Error creating directory?\n"; */
	cerr << "Writing results to \"" << fileName << "\"\n";
	ofstream output;	
	/* output.precision(3); */
	output.open(fileName.c_str());
	output   << " " << '\t'
		<< "density" << '\t'
		<< "u_1" << '\t'
		<< "u_2" << '\t'
		<< "u_3" << '\t'
		<< "sigma_11" << '\t'
		<< "sigma_12" << '\t'
		<< "sigma_13" << '\t'
		<< "sigma_22" << '\t'
		<< "sigma_23" << '\t'
		<< "sigma_33" << '\t'
		<< "I1" << '\t'
		<< "I2" << '\t' 
		<< "I3" << '\t'
		<< "dI_1" << '\t'
		<< "dI_2" << '\t'
		<< "dI_3" << endl;

	vector <double> out;
	for(int i = gCs; i < cells + gCs; i++)
	{
		U[i].prepareOutputs();
		/* output << (double)(i-2)/cells/100.; */
		/* out.push_back(U[i].rho_()/1000.); */
		/* out.push_back(U[i].u_(0)/1000.); */
		/* out.push_back(U[i].u_(1)/1000.); */
		/* out.push_back(U[i].u_(2)/1000.); */
		/* out.push_back(U[i].sigma_(0,0)/1e9); */
		/* out.push_back(U[i].sigma_(0,1)/1e9); */
		/* out.push_back(U[i].sigma_(0,2)/1e9); */
		/* out.push_back(U[i].sigma_(1,1)/1e9); */
		/* out.push_back(U[i].sigma_(1,2)/1e9); */
		/* out.push_back(U[i].sigma_(2,2)/1e9); */
		/* out.push_back(U[i].S_()); */
		/* out.push_back(U[i].getInvariant(0)); */
		/* out.push_back(U[i].getInvariant(1)); */
		/* out.push_back(U[i].getInvariant(2)); */

		output  << (double)(i-2)/cells/100. << '\t'
			<< U[i].rho_()/1e3 << '\t' 
			<< U[i].u_(0)/1000.<< '\t' 
			<< U[i].u_(1)/1000.<< '\t' 
			<< U[i].u_(2)/1000.<< '\t' 
			<< U[i].sigma_(0,0)/1e9<< '\t' 
			<< U[i].sigma_(0,1)/1e9<< '\t' 
			<< U[i].sigma_(0,2)/1e9 << '\t'
			<< U[i].sigma_(1,1)/1e9 << '\t'
			<< U[i].sigma_(1,2)/1e9 << '\t' 
			<< U[i].sigma_(2,2)/1e9 << '\t'
			<< U[i].S_() << '\t'
			<< U[i].getInvariant(0) << '\t'
			<< U[i].getInvariant(1) << '\t'
			<< U[i].getInvariant(2) << '\t' << endl;
		for(int j = 0; j < out.size(); j++)
		{
			output << '\t' << out[j];
		}
		output << endl;
		out.clear();
	}
	/* for(int i = gCs; i < cells + gCs; i++) */
	/* { */						
	/* 	output << (double)(i-2)/cells << '\t' << U[i].rho_() << '\t' << U[i].u_(0) << '\t' << U[i].u_(1) << '\t' << U[i].u_(2) << '\t' << U[i].sigma_(0,0) << '\t' << U[i].sigma_(0,1); */  
	/* for(int i = 0; i < 6; i++) */
	/* { */
	/* 	output << '\t' << U[i].getInvariant(i); */		
	/* } */
	/* output << endl; */

	/* } */
	output.close();	
}

