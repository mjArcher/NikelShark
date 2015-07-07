
#include <libconfig.h++>
#include <iomanip>
#include <cstdlib>

#include "ElasticState.h"
#include "ElasticPrimState.h"
#include "SolidSystem.h"

using namespace std;
using namespace libconfig;
using namespace Eigen;


#define debug_

/* int e_size = -27; */

/* void initialize(Material, Domain); */

/* double internalEnergy(Matrix3d, double); */

/* void setBCs(Material); */

/* /1* void printArray(Material); *1/ */

enum slope_lim{superbee, minbee, van_leer, albalda};//0,1,2,3

enum ICType { Barton1, Barton2, Barton3, Barton4 };

ICType icType;

slope_lim sl;

struct Domain {
	int Ni; //number of x cells
	int GNi; // number of x cells plus ghost cells
	int GC; // number of ghost cells
	int starti; //start of computational domain
	int endi; //end of computational domain
	double Lx; // length of domain
	double dx; // cell size

	vector<ElasticState> sol;
	
	Domain(const int a_Ni, const int a_GC, const double a_Lx) :
		Ni(a_Ni), GNi(a_Ni+2*a_GC), GC(a_GC),	starti(a_GC), endi(a_Ni+a_GC),
		Lx(a_Lx), dx(a_Lx/a_Ni){
			sol = vector<ElasticState>(GNi);
		}
	Domain(){}
	~Domain(){
		/* delete sol[]; */	
	}
};

//use kevin's domain and Material structs as a base to modify later
//work by solving over each material in turn
//(How does 1D ghost fluid work?)
//two or more separate domains which we solve then repopulate 'solution vector'
//system, dom and EOS associated with each material

struct Material {
  const System sys;
  vector<ElasticState> sol;  
  const Domain dom;
  string name;

  Material(const string& a_name, const Domain& a_dom, const ElasticEOS Eos) :
    sys(Eos), 
    sol(a_dom.GNi), dom(a_dom), name(a_name) {}
};


double slopelim(double);

double ksi_r(double);

ElasticState grad(Material, int);


void BCs(Material& mat) {  

/* i boundaries */
#pragma omp parallel for schedule(dynamic)
  for(int i=0; i<mat.dom.starti; ++i) {
    const int imagei = 2*mat.dom.starti-1-i;
    mat.sol[i] = mat.sol[imagei];
  }

#pragma omp parallel for schedule(dynamic)
  for(int i=mat.dom.endi; i<mat.dom.GNi; ++i) {
    const int imagei = 2*mat.dom.endi-1-i;
    mat.sol[i] = mat.sol[imagei];
  }  
}


void ICInterface(Material& mat,
		 const double iface,
		 const ElasticState& left, const ElasticState& right) 
{
	for(int i=mat.dom.starti; i<mat.dom.endi; ++i) {
		const double x = (i-mat.dom.starti+0.5)*mat.dom.dx;
    if(x<iface) 
    {      
      mat.sol[i] = left;
    }
    else 
    {
      mat.sol[i] = right;
    } 
  }
}

void setBCs(Material& mat)
{
	// Transmissive
	for (int i = 0; i < mat.dom.GC; i++)
	{
		mat.sol[i] = mat.sol[2*mat.dom.GC - i - 1];
		mat.sol[mat.dom.Ni+  i + mat.dom.GC] = mat.sol[mat.dom.GNi - mat.dom.GC - i - 1]; 
	}					
}	
/* void ICInterface(Material mat, const double iface, ElasticState left, ElasticState right){ */
/* 	unsigned int i; */
/* 	for(i = gCs; i < cells*disctyX + gCs; i++) */
/* 	{ */
/* 		dom[i] = ElasticStateL; */		
/* 	} */
/* 	for(i; i < N - gCs; i++) */
/* 	{ */	
/* 		dom[i] = ElasticStateR; */
/* 	} */	
/* 	setBCs(domain); */

/* } */ 

void solveXGodunov()
{
	//0. Convert to primitive state and estimate primitive state at the i+1/2 position
	//1. construct acoustic tensor
	//2. Perform eigen decomposition to get eigen values and ortogonal marix
	//3. construct left and right eigen vectors
	//4. Solve linear system by multiplying by one of these eigenvectors to get coefficients
	//5. Reconstruct primitive state
	//6. He Presto

}

ElasticState forceFlux(System sys, vector<ElasticState> left, vector<ElasticState> right, double dt_dX, int i)
{
	//test
	/* const ElasticState eA; */
	/* const ElasticState eB; */
	/* const ElasticState fA = sys.flux(eA); */	
	/* const ElasticState test0 = eA + eB; */
	/* ElasticState test = sys.flux(left[i]) + sys.flux(right[i]); */

	const ElasticState F_lf = 0.5*(sys.flux(right[i]) + sys.flux(left[i+1])) + 0.5*(1./dt_dX)*(right[i]-left[i+1]);
	const ElasticState C_r = 0.5*(right[i] + left[i+1]) + 0.5*(dt_dX)*(sys.flux(right[i]) - sys.flux(left[i+1]));
	ElasticState Fforce = 0.5*(F_lf + sys.flux(C_r));

	/* const ElasticState F_lf = 0.5*(right[i].flux() + left[i+1].flux()) + 0.5*(1./dt_dX)*(right[i] - left[i+1]); */
	/* ElasticState C_r = 0.5*(right[i] + left[i+1]) + 0.5*(dt_dX)*(right[i].flux() - left[i+1].flux()); */
	/* ElasticState Fforce = 0.5*(F_lf + C_r.flux()); */
	return Fforce;
}


void solveX(Material& mat, const double dt)
{
	int cells = mat.dom.GNi;
	vector<ElasticState> left(cells);
	vector<ElasticState> right(cells);
	//primitive state vector
	vector<ElasticPrimState> prim(cells);

	const double dt_dX = dt/mat.dom.dx;

	for(int i = mat.dom.starti - 1; i < mat.dom.endi + 1; ++i) 
	{
		//2. calculate left and right extrapolated values (bar) pg 514 and slope pg 506.
		const ElasticState slopebar = grad(mat, i);
		ElasticState Cleft = mat.sol[i] - 0.5 * slopebar;
		ElasticState Cright = mat.sol[i] + 0.5 * slopebar;
		const ElasticState Cbar = 0.5 * (dt_dX) * (mat.sys.flux(Cleft) - mat.sys.flux(Cright));
/* Cleft.flux() - Cright.flux()); //change these flux functions */
		left[i] = Cleft + Cbar;
		right[i] = Cright + Cbar;
	}		
	
	left[0] = mat.sol[0];
	right[0] = mat.sol[0];//is required
	left[cells - 1] = mat.sol[cells - 1];
	right[cells - 1] = mat.sol[cells - 1];	

	System system = mat.sys;	
	//3. calculate force flux using LF and RI, and calculate new cell averaged Ui pg 494
	for(int i = mat.dom.starti - 1; i < mat.dom.endi + 1; i++)
	{
		mat.sol[i] += dt_dX*(forceFlux(system, left, right, dt_dX, i - 1) - forceFlux(system, left, right, dt_dX, i));
	}
	//printArray(U);
	setBCs(mat);
}

//Force flux calculation pg 512
ElasticState grad(const Material mat, int i)
{
	//1. calculate r
	const double w = 1; // check this 
	ElasticState num = mat.sol[i] - mat.sol[i-1];
	ElasticState den = mat.sol[i+1] - mat.sol[i];
	ElasticState r;
	ElasticState delta = 0.5*(1+w)*num;// + 0.5*(1-w)*denom;
	//aNew = Ui[i].soundSpeed();	
	//cout << i << endl;	
	for(int j = 0; j < 3; j++)
	{	
		if(num[j] == 0 && den[j] == 0)
		{
			r[j] = 1;
		}
		else if(den[j] == 0)
		{
			den[j] = 1.e-10;
			r[j] = num[j]/den[j];			
		}
		else 
			r[j] = num[j]/den[j];		
	}
	//ElasticState r = num/denom;
	ElasticState deltaold = delta;	
	for(int j = 0; j < ElasticState::e_size; j++)
	{
		const double ksi = slopelim(r[j]);
		//ksi = 1;
		delta[j] *= ksi; 	 
	}
	return delta; 
} 


double getMinDt(const Material& mat)
{
	double Smax = 0, a, u_1, dt;
	a = 10000; //o be calculated
	/* double c = ElasticState::soundSpeed(); */
	//! Convert to primitive 


	for (int i = 0; i < mat.dom.GNi; i++) //N
	{
		ElasticPrimState prim = mat.sys.conservativeToPrimitive(mat.sol[i]);
		double primVel = prim[0];	
		u_1 = fabs(primVel);
		double rho = mat.sys.Density(prim);

		if(rho > 0)
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

	dt = mat.dom.dx/Smax;

	if(dt == 0)
	{
		cout << "dt == 0: exit" << endl;
		exit(1);
	}	
	return dt;
	/* cout << "Smax " << Smax << " dt " << dt << endl; */
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

void printArray(Material mat)
{
	/* inline std::ostream& operator<<(std::ostream& os, ElasticState& param) */
	/* { */
	/* param.initialStates(os); */
	/* os << "Current states \n" << */
	/* 	  "Density\t"    << param.rho_() << "\n" */
	/* 	  "Energy\t"     << param.E_() << "\n" */
	/* 	  "U or F vector\t " << "\n"; */
	/* os.precision(3); */
	/* for (int i = 0; i < ElasticState::e_size; i++) */
	/* { */
	/* 	os << setw(5) << " |" << setw(5) << param[i]; */
	/* } */	
	/* os << setw(5) << " | "; */
	/* return os; */
	/* } */

/* const string name[] = {"rhou_1", "rhou_2", "rhou_3", "rhoF_11", "rhoF_12", "rhoF_13", "rhoF_21", "rhoF_22", "rhoF_23", "rho_31", "rho_32", "rho_33", "rhoE"}; */
	
	cout.precision(4);
	/* for(int i = 0; i < ElasticState::e_size; i++) */
	/* { */
	/* 	cout << setw(3) << " |" << setw(7) << name[i]; */
	/* } */
	/* cout << setw(7) <<  " |" << endl; */
	for(int i = 0; i < 140; i++)
	{
		cout << '-';
	}
	cout << endl;
	for(int i = 0; i < mat.dom.GNi; i++)
	{
		cout << setw(7) << left << i << ' ' << mat.sol[i] << endl;
	}

}


/**
 * Calculate the slope limiter and slope delta_i.
 * Update delta_i.
 * @param Ui Pointer to solution domain Ui
 * @param i Domain index of conserved vector  
 * @return updated slope vector delta_i
 */

// delta ibar slope on pg 509

// slope limiters pg 509 and 510
double slopelim(double r)
{
	//calculate ksi_r
	double ksi_sl = 0;
	//slope limiters:
	//cout << "Slope limiter " << slLimChoice << endl;
	switch(sl)
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
	double w = 1;
	return 2.*beta_fw/(1. - w + (1. + w)*r);
}

void outputGnu(string file, Material mat, int outStep)
{
	/* int ret = system(("mkdir -p " + fileName).c_str()); */
	/* if(ret!=0) cerr << "Error creating directory?\n"; */
	cerr << "Writing results to \"" << file << "\"\n";
	ofstream output;	
	/* output.precision(3); */
	output.open(file.c_str());
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
		<< "S " << '\t'
		<< "I1" << '\t'
		<< "I2" << '\t' 
		<< "I3" << '\t'
		<< "dI_1" << '\t'
		<< "dI_2" << '\t'
		<< "dI_3" << endl;
//some work needs doing here
	vector <double> out;
	for(int i = mat.dom.starti; i < mat.dom.endi; i++)
	{
		//stress (system), velocity (primState), entropy (EOS), invariants
		const ElasticState consState = mat.sol[i];	
		double rho = mat.sys.Density(consState);
		const ElasticPrimState primState = mat.sys.conservativeToPrimitive(consState);
		double entropy = primState.S_();
		Vector3d u = primState.u_();
		Vector3d inv = mat.sys.getInvariants(primState.F_());
		Matrix3d sigma = mat.sys.stress(primState);

		output  << (double)(i-2)/mat.dom.Ni/100. << '\t'
			<< rho/1e3 << '\t' 
			<< u(0)/1000.<< '\t' 
			<< u(1)/1000.<< '\t' 
			<< u(2)/1000.<< '\t' 
			<< sigma(0,0)/1e9 << '\t'
			<< sigma(0,1)/1e9 << '\t'
			<< sigma(0,2)/1e9 << '\t'
			<< sigma(1,1)/1e9 << '\t'
			<< sigma(1,2)/1e9 << '\t'
			<< sigma(2,2)/1e9 << '\t'
			<< entropy << '\t'
			<< inv[0] << '\t'
			<< inv[1] << '\t'
			<< inv[2] << '\t' << endl;


		/* << mat.sol[i].rho_()/1e3 << '\t' */ 
			/* << [i].u_(0)/1000.<< '\t' */ 
			/* << mat.sol[i].u_(1)/1000.<< '\t' */ 
			/* << mat.sol[i].u_(2)/1000.<< '\t' */ 
			/* << mat.sol[i].sigma_(0,0)/1e9<< '\t' */ 
			/* << mat.sol[i].sigma_(0,1)/1e9<< '\t' */ 
			/* << mat.sol[i].sigma_(0,2)/1e9 << '\t' */
			/* << mat.sol[i].sigma_(1,1)/1e9 << '\t' */
			/* << mat.sol[i].sigma_(1,2)/1e9 << '\t' */ 
			/* << mat.sol[i].sigma_(2,2)/1e9 << '\t' */
			/* << mat.sol[i].S_() << '\t' */
			/* << mat.sol[i].getInvariant(0) << '\t' */
			/* << mat.sol[i].getInvariant(1) << '\t' */
			/* << mat.sol[i].getInvariant(2) << '\t' << endl; */

		for(unsigned int j = 0; j < out.size(); j++)
		{
			output << '\t' << out[j];
		}
		output << endl;
		out.clear();
	}
	output.close();	
}

void initialize(Domain &dom, double &tf, double &outFreq)
{
	/* double rho_0L = 8930, rho_0R = 8930, eL(0), eR(0); */

	/* ElasticState flux = ElasticStateL.flux(); */

	// put this code into new function which takes domainInfo and solution vector as parameters
	//Populate array 
	//print out initial domain
	// U vector 
	/* printArray(domain); */
}

void unitTests(Material mat, ElasticPrimState iPrimL, ElasticPrimState iPrimR)
{
	cout << "perform simple solution array tests " << endl;
	/* printArray(mat); */

#ifdef debug_

	mat.sys.Eos.checkEosConstants();

	cout << "Prim stateL" << endl;
	cout << iPrimL << endl;
	cout << "Density" << mat.sys.Density(iPrimL) << endl;
	cout << "Strain tensor L " << mat.sys.strainTensor(iPrimL.F_()) << endl;
	cout << "INVARIANTS L " << endl;
	cout << mat.sys.getInvariants(iPrimL.F_()) << endl;	
	cout << " -----------------------------------------------------------" << endl;

	cout << "Prim stateR" << endl;
	cout << iPrimR << endl;
	cout << "Density" << mat.sys.Density(iPrimR) << endl;
	cout << "Strain tensor L " << mat.sys.strainTensor(iPrimR.F_()) << endl;
	cout << "INVARIANTS R " << endl;
	cout << mat.sys.getInvariants(iPrimR.F_()) << endl;	
	cout << " -----------------------------------------------------------" << endl;

	ElasticState iStateL = mat.sys.primitiveToConservative(iPrimL);
	ElasticState iStateR = mat.sys.primitiveToConservative(iPrimR);

	cout <<"Cons stateL" << endl;
	cout << iStateL << endl;
	cout << "Density" << mat.sys.Density(iStateL) << endl;
	cout << endl;
	cout << " -----------------------------------------------------------" << endl;

	cout <<"Cons stateR" << endl;
	cout << iStateR << endl;
	cout << "Density" << mat.sys.Density(iStateR) << endl;
	cout << endl;
	cout << " -----------------------------------------------------------" << endl;

	ElasticPrimState iPrimL2 = mat.sys.conservativeToPrimitive(iStateL);
	ElasticPrimState iPrimR2 = mat.sys.conservativeToPrimitive(iStateR);

	cout << "Prim stateL" << endl;
	cout << iPrimL2 << endl;
	cout << "Density" << mat.sys.Density(iPrimL2) << endl;
 	cout << endl;
	cout << "INVARIANTS L " << endl;
	cout << mat.sys.getInvariants(iPrimL.F_()) << endl;	

	cout << " -----------------------------------------------------------" << endl;
	cout << "Prim stateR" << endl;
	cout << iPrimR2 << endl;
	cout << "Density" << mat.sys.Density(iPrimR2) << endl;
	cout << endl;
	cout << "INVARIANTS R " << endl;
	cout << mat.sys.getInvariants(iPrimR.F_()) << endl;	
	cout << " -----------------------------------------------------------" << endl;

//---------------------
//

	cout << "Check EOS constants" << endl;
	mat.sys.Eos.checkEosConstants();
#endif
}

//Elastic1D barton1 100 bartonInitialOut bartonFinalOut
int main(int argc, char ** argv)
{
	string outFile, outInitFile;
	int Ni;
	if(argc != 5){
		cout << "Insufficient Parameters: exit" << endl;
		cout << "Number of parameters " << argc << endl;
		exit(1);
	}
	else {
		const string icStr = string(argv[1]);
		cout << "Solve Test " << icStr << endl;
		if(icStr == "barton1") 
		{
			icType = Barton1;
		}
		else if(icStr == "barton2")
		{
			icType = Barton2;
		}
		else if(icStr == "barton3")
		{
			icType = Barton3;	
		}
		else 
		{
			cerr << "Unknown problem type \'" << icStr << "\' - options are barton1, barton2, barton3,\n";
			exit(1);
		}
  }
	/* * Setting& setting = cfg.lookup("tests.test1"); *1/ */
	double CFL = 0.6;

	cout << " Solve barton test with " << Ni << " cells " << " and output to " << outInitFile << " and " << outFile << endl;
	Domain dom;

	string matStr;
	
	/* initialize(dom, tf, Ni, outFreq, matstr); */
	
	Material * mat;
	double iface;
	Matrix3d FL, FR;
	Vector3d uL, uR;
	double SL, SR;

	double t, tf, dt, outFreq;

	if(icType == Barton1)
	{
		tf = 6.0e-7;                 // Final time
		outFreq = 1e-1;            // Output frequency    
		iface = 0.005;	
		dom = Domain(Ni, 2, 0.01); // Domain
		matStr = "copper";

		FL <<	 0.98, 0, 0, 0.02, 1., 0.1, 0, 0, 1.;
		uL << 0, 500., 1000.;		
		SL = 1000.;

		FR << 1., 0, 0, 0, 1., 0.1, 0, 0, 1.;
		uR << 0, 0, 0;
		SR = 0.;
		// new Romenski
	}
	if(icType == Barton2)
	{
		tf = 6.0e-7;                 // Final time
		outFreq = 1e-1;            // Output frequency    
		iface = 0.005;	
		dom = Domain(Ni, 2, 0.01); // Domain
		matStr = "copper";

		FL << 1., 0, 0,	-0.01, 0.95, 0.02, -0.015, 0, 0.9;
		uL << 2000., 0, 100.;
		SL = 0;

		FR << 1., 0, 0, 0.015, 0.95, 0, -0.01, 0, 0.9;
		uR << 0, -30., -10.;
		SR = 0;
		// new romenski
	}

#ifdef debug_
	cout << "Completed initial conditions " << endl;
#endif

	//create primitive states first then convert to conservative (using system function)
	// do not need density in Eos for this (or material) since we can convert to 
	
	if(matStr == "copper"){
		ElasticEOS Eos("copper");
		mat = new Material("copper", dom, Eos);
	}
	else {
		cout << "No valid material" << endl;
		exit(1);
	}

	ElasticPrimState iPrimL(uL, FL, SL); //think about const safety for these
	ElasticPrimState iPrimR(uR, FR, SR); 

	ElasticState iStateL = mat->sys.primitiveToConservative(iPrimL);
	ElasticState iStateR = mat->sys.primitiveToConservative(iPrimR);

	ICInterface(*mat, iface, iStateL, iStateR);
	

#ifdef debug_

	unitTests(*mat, iPrimL, iPrimR); //check all initialized parameters 	/* output(dom,outputInitial); */
	printArray(*mat);
	/* exit(1); */
#endif

	t = 0;
	int step = 0, outStep = 0; 

	/* output initially */
	outputGnu(outInitFile, *mat, outStep);

	while(t < tf)
	{
		BCs(*mat);
		//1. calculate time step: CFL and boundary conditions pg 495
		dt = getMinDt(*mat);
		/* if(step < 10)	{dt /= 10;} */
		dt *= CFL;

		t += dt;		

    if(t < outStep * outFreq && t + dt > outStep * outFreq) {
      dt = outStep * outFreq - t;
    }
    
    if(t >= outStep * outFreq) 
    {
      outputGnu(outFile, *mat, outStep);
      cerr << "Saved results to output file " << outStep << "\n";
      ++outStep;
    }

		/* if(t > tf) */
		/* { */				
		/* 	t -= dt; */
		/* 	dt = tf - t; */
		/* 	t += dt; // should equal tf */
		/* } */

		solveX(*mat, dt);

		++step;

		// output every time step
		cout.precision(3);

		cout << " [" << outStep << "] " << setw(6) << step << setw(6) << "Time " << setw(6) << t << 
						setw(15) << " Remaining " << setw(6) <<tf-t<< endl;
		/* break; */
	}
	outputGnu(outFile, *mat, outStep);
#ifdef debug_
	printArray(*mat);
#endif

  cout << outInitFile << " " << outFile << endl;
  return 0;
  
  }		


