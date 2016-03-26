#include <libconfig.h++>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <omp.h>
#include <fenv.h>
#include <ctime>

#include "ElasticState.h"
#include "ElasticPrimState.h"
#include "SolidSystem.h"
#include "InputSolid.h"

using namespace std;
using namespace libconfig;
using namespace Eigen;

/* #define debug_ */

//use kevin's domain and Material structs as a base to modify later
//work by solving over each material in turn
//(How does 1D ghost fluid work?)
//two or more separate domains which we solve then repopulate 'solution vector'
//system, dom and EOS associated with each material

enum slope_lim{superbee, minbee, vanleer, albalda};//0,1,2,3
slope_lim sl;

struct Domain {
	int Ni; //number of x cells
	int GNi; // number of x cells plus ghost cells
	int GC; // number of ghost cells
	int starti; //start of domain
	int endi; //end of domain
	double Lx; // length of domain
	double dx; // cell size

	Domain(const int a_Ni, const int a_GC, const double a_Lx) :
		Ni(a_Ni), GNi(a_Ni+2*a_GC), GC(a_GC),	starti(a_GC), endi(a_Ni+a_GC),
		Lx(a_Lx), dx(a_Lx/a_Ni){}
  Domain(){}
	~Domain(){
	}
};

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
ElasticState grad(const Material&, int);

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
    const ElasticPrimState& left, const ElasticPrimState& right) 
{
  //convert to conserved states
  ElasticState stateL = mat.sys.primitiveToConservative(left);
  ElasticState stateR = mat.sys.primitiveToConservative(right);

  for(int i=mat.dom.starti; i<mat.dom.endi; ++i) {
    const double x = (i-mat.dom.starti+0.5)*mat.dom.dx;
    if(x<iface) 
    {      
      mat.sol[i] = stateL;
    }
    else 
    {
      mat.sol[i] = stateR;
    } 
  }
}

void solveXGodunov()
{
	//0. Convert to primitive state and estimate primitive state at the i+1/2 position
	//1. construct acoustic tensor
	//2. Perform eigen decomposition to get eigen values and ortogonal marix
	//3. construct left and right eigen vectors
	//4. Solve linear system by multiplying by one of these eigenvectors to get coefficients
	//5. Reconstruct primitive state
	//6. Hey Presto
}

ElasticState forceFlux(System& sys, vector<ElasticState>& left, vector<ElasticState>& right, double dt_dX, int i)
{
	const ElasticState F_lf = 0.5*(sys.flux(right[i]) + sys.flux(left[i+1])) + (0.5/dt_dX)*(right[i]-left[i+1]);
	const ElasticState C_r = 0.5*(right[i] + left[i+1]) + 0.5*(dt_dX)*(sys.flux(right[i]) - sys.flux(left[i+1]));
	ElasticState F_force = 0.5*(F_lf + sys.flux(C_r));
	return F_force;
}


void solveXForce(Material& mat, const double dt)
{

}

void solveX(Material& mat, const double dt)
{
	int cells = mat.dom.GNi;
	vector<ElasticState> left(cells);
	vector<ElasticState> right(cells);
	const double dt_dX = dt/mat.dom.dx;

  #pragma omp parallel for schedule(dynamic)
	for(int i = mat.dom.starti - 1; i < mat.dom.endi + 1; ++i) 
	{
		//extrapolated values (bar) pg 514 and slope pg 506.
		const ElasticState slopebar = grad(mat, i);
		ElasticState Cleft = mat.sol[i] - 0.5 * slopebar;
		ElasticState Cright = mat.sol[i] + 0.5 * slopebar;
		const ElasticState Cbar = 0.5 * (dt_dX) * (mat.sys.flux(Cleft) - mat.sys.flux(Cright));
/* Cleft.flux() - Cright.flux()); //change these flux functions */
		left[i] = Cleft + Cbar;
		right[i] = Cright + Cbar;
	}		
	
	left[0] = mat.sol[0];
	right[0] = mat.sol[0]; //is required
	left[1] = mat.sol[1];
	right[1] = mat.sol[1]; //is required
	left[cells - 1] = mat.sol[cells - 1];
	right[cells - 1] = mat.sol[cells - 1];	
	left[cells - 2] = mat.sol[cells - 2];
	right[cells - 2] = mat.sol[cells - 2];	

	System sys = mat.sys;	
	//3. calculate force flux using LF and RI, and calculate new cell averaged Ui pg 494
  #pragma omp parallel for schedule(dynamic)
	for(int i = mat.dom.starti - 1; i < mat.dom.endi + 1; i++)
	{
		mat.sol[i] += dt_dX*(forceFlux(sys, left, right, dt_dX, i - 1) - forceFlux(sys, left, right, dt_dX, i));
	}
	//printArray(U);
	BCs(mat);
}

//Force flux calculation pg 512
ElasticState grad(const Material& mat, int i)
{
	//1. calculate r
	const double w = 1; // check this 
	ElasticState num = mat.sol[i] - mat.sol[i-1];
	ElasticState den = mat.sol[i+1] - mat.sol[i];
	ElasticState r;
	ElasticState delta = 0.5*(1+w)*num;// + 0.5*(1-w)*denom;
	//aNew = Ui[i].soundSpeed();	
	//cout << i << endl;	
	for(int j = 0; j < ElasticState::e_size; j++)
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

//need to parallelise
double getMinDt(const Material& mat)
{
	double Smax = 0, a, u_1;
	a = 5500; //to be calculated
	/* double c = ElasticState::soundSpeed(); */
	//! Convert to primitive 

  double sharedmindt = std::numeric_limits<double>::max(); 
/* #pragma omp parallel */
  {
    double mindt = std::numeric_limits<double>::max();
/* #pragma omp for nowait schedule(dynamic) */
    for (int i = mat.dom.starti; i < mat.dom.endi; i++) 
    {
      /* ElasticPrimState prim = mat.sys.conservativeToPrimitive(mat.sol[i]); */
      /* double primVel = prim[0]; */	
      /* u_1 = fabs(primVel); */
      /* double rho = mat.sys.Density(prim); */

      ElasticState state = mat.sol[i];
      double momx = state[0];
      double rho = mat.sys.Density(state);
      double primVel = momx/rho;

      u_1 = fabs(primVel);
      /* double rho = mat.sys.Density(prim); */
      

      if(rho > 0) {
        if((a + u_1) > Smax)
        {
          Smax = a + u_1; 
        	mindt = mat.dom.dx/Smax;
        }
      }		 
      else { cout << "density is zero or negative: exit" << endl; exit(1);	}
    }	
/* #pragma omp critical */
    {
      sharedmindt = std::min(sharedmindt, mindt);
    }
  }

  return sharedmindt;
/* 	if(dt == 0) */
/* 	{ */
/* 		cout << "dt == 0: exit" << endl; */
/* 		exit(1); */
/* 	} */	
	/* cout << "Smax " << Smax << " dt " << dt << endl; */
}

void printArray(Material mat)
{
	cout.precision(4);
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
  /* cout << "Slope limiter " << sl << endl; */ 
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
		case vanleer:		
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
	double beta_fw = 1.;//2/(0.1);
	double w = 1.;
	return 2.*beta_fw/(1. - w + (1. + w)*r);
}

void outputGnu(string file, Material mat, int outStep, double t)
{
	/* int ret = system(("mkdir -p " + fileName).c_str()); */
	/* if(ret!=0) cerr << "Error creating directory?\n"; */
	cerr << "Writing results to \"" << file << "\"\n";
	ofstream output;	
	/* output.precision(3); */
  if(outStep == 0)
  {
    output.open(file.c_str(), ios::app);
    output << "# t = 0 " << '\n'
      << "# Column 1: x-coordinate" << '\n'
      << "# Column 2: density" << '\n'
      << "# Column 3: U" << '\n'
      << "# Column 4: V" << '\n'
      << "# Column 5: W" << '\n'
      << "# Column 6: sigma_11" << '\n'
      << "# Column 7: sigma_12" << '\n'
      << "# Column 8: sigma_13" << '\n'
      << "# Column 9: sigma_21" << '\n'
      << "# Column 10: sigma_22" << '\n'
      << "# Column 11: sigma_23" << '\n'
      << "# Column 12: sigma_31" << '\n'
      << "# Column 13: sigma_32" << '\n'
      << "# Column 14: sigma_33" << '\n'
      << "# Column 15: S " << '\n'
      << "# Column 16: I1" << '\n'
      << "# Column 17: I2" << '\n' 
      << "# Column 18: I3" << '\n'
      << "# Column 19: dI_1" << '\n'
      << "# Column 20: dI_2" << '\n'
      << "# Column 21: dI_3" << '\n' 
	    << "# Column 22: G11" << '\n'  
	    << "# Column 22: G12" << '\n'  
	    << "# Column 22: G13" << '\n'  
	    << "# Column 22: G21" << '\n'  
	    << "# Column 22: G22" << '\n'  
	    << "# Column 22: G23" << '\n'  
	    << "# Column 22: G31" << '\n'  
	    << "# Column 22: G32" << '\n'  
	    << "# Column 22: G33" << '\n' << endl;
    output.close();
  }
  //some work needs doing here
	vector <double> out;
  output.open(file.c_str(), ios::app);
  if(outStep != 0)
    output << "# t = " << t << endl;
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
  	const Matrix3d G = mat.sys.strainTensor(primState.F_());
    // convert cell into 
    // phys domain is mat.dom.Lx/cell 
    // (i-2)*dx + dx/2
    // dx*((i-mat.dom.starti) + 0.5)
    //get invariants
    
    
		output << (double)mat.dom.dx*((i-mat.dom.starti)+0.5) << '\t'
			<< rho/1e3 << '\t' 
			<< u(0)/1000.<< '\t' 
			<< u(1)/1000.<< '\t' 
			<< u(2)/1000.<< '\t' 
			<< sigma(0,0)/1e9 << '\t'
			<< sigma(0,1)/1e9 << '\t'
			<< sigma(0,2)/1e9 << '\t'
			<< sigma(1,0)/1e9 << '\t'
			<< sigma(1,1)/1e9 << '\t'
			<< sigma(1,2)/1e9 << '\t'
			<< sigma(2,0)/1e9 << '\t'
			<< sigma(2,1)/1e9 << '\t'
			<< sigma(2,2)/1e9 << '\t'
			<< entropy/1e6 << '\t'
			<< inv[0] << '\t'
			<< inv[1] << '\t'
			<< inv[2] << '\t' 
      << G(0,0) << '\t'
      << G(0,1) << '\t'
      << G(0,2) << '\t'
      << G(1,0) << '\t'
      << G(1,1) << '\t'
      << G(1,2) << '\t'
      << G(2,0) << '\t'
      << G(2,1) << '\t'
      << G(2,2) << '\t'
      << endl;

		/* for(unsigned int j = 0; j < out.size(); j++) */
		/* { */
		/* 	output << '\t' << out[j]; */
		/* } */
		/* out.clear(); */
	}
  output << '\n' << std::endl;
	output.close();	
}

void advance(Material& mat, const double dt)
{
// apply curl constraint (2D) // geometric -
// cylindrical //spherical bcs //plasticity
  //series of advance functions: levelset, geometric bcs, 
  solveX(mat, dt);
}


int solveSystem(InputSolid inputSolid, Material* mat){
  
  double t(0), dt(0), tend(inputSolid.end_time);
  const double CFL(inputSolid.input_CFL);
  //output names 
  double outFreq(tend/inputSolid.frequency);
	int step(0), outStep(0); 
  std::string outDir(inputSolid.filePath), outName(inputSolid.fileName);
  std::string outFile(outDir + outName);
  //delete existing output
  ofstream myfile;
  myfile.open(outFile, ios::trunc);
  myfile.close();

	while(t < tend)
	{
		BCs(*mat);
		//1. calculate time step: CFL and boundary conditions pg 495
		dt = getMinDt(*mat);
		/* if(step < 20)	{dt /= 10;} */
		dt *= CFL;


    if(t < outStep * outFreq && t + dt > outStep * outFreq) {
      dt = outStep * outFreq - t;
    }
    
    if(t >= outStep * outFreq) 
    {
      outputGnu(outFile, *mat, outStep, outStep*outFreq);
      cerr << "Saved results to output file " << outStep << "\n";
      ++outStep;
    }

		/* if(t > tf) */
		/* { */				
		/* 	t -= dt; */
		/* 	dt = tf - t; */
		/* 	t += dt; // should equal tf */
		/* } */


		cout.precision(4);
		cout << " [" << outStep << "] " << setw(6) << step << setw(6) << "Time " << setw(6) << t << " dt " << setw(6) << dt << 
						setw(15) << " Remaining " << setw(6) <<tend-t<< endl;

    advance(*mat, dt);

		t += dt;
    /* if(step == 1) */
    /* { */
    /*   outputGnu(outFile, *mat, outStep, t); */
    /*   break; */
    /* } */

		++step;
	}

	outputGnu(outFile, *mat, outStep, tend); //final output

  //total time
#ifdef debug_
	printArray(*mat);
#endif
  return 0;
}

void unitTests(Material mat, ElasticPrimState iPrimL, ElasticPrimState iPrimR)
{
	cout << "perform simple solution array tests " << endl;
	/* printArray(mat); */

#ifdef debug_
	cout.precision(6);

	mat.sys.Eos.checkEosConstants();

	cout << "Prim stateL" << endl;
	cout << iPrimL << endl;
	cout << "Density" << mat.sys.Density(iPrimL) << endl;
	cout << "Strain tensor L " << mat.sys.strainTensor(iPrimL.F_()) << endl;
	cout << "INVARIANTS L " << endl;
	cout << mat.sys.getInvariants(iPrimL.F_()) << endl;	

	cout << "Prim stateR" << endl;
	cout << iPrimR << endl;
	cout << "Density" << mat.sys.Density(iPrimR) << endl;
	cout << "Strain tensor L " << mat.sys.strainTensor(iPrimR.F_()) << endl;
	cout << "INVARIANTS R " << endl;
	cout << mat.sys.getInvariants(iPrimR.F_()) << endl;	

	ElasticState iStateL = mat.sys.primitiveToConservative(iPrimL);
	ElasticState iStateR = mat.sys.primitiveToConservative(iPrimR);

	cout <<"Cons stateL" << endl;
	cout << iStateL << endl;
	cout << "Density" << mat.sys.Density(iStateL) << endl;
	cout << endl;

	cout <<"Cons stateR" << endl;
	cout << iStateR << endl;
	cout << "Density" << mat.sys.Density(iStateR) << endl;
	cout << endl;

	ElasticPrimState iPrimL2 = mat.sys.conservativeToPrimitive(iStateL);
	ElasticPrimState iPrimR2 = mat.sys.conservativeToPrimitive(iStateR);

	cout << "Prim stateL" << endl;
	cout << iPrimL2 << endl;
	cout << "Density" << mat.sys.Density(iPrimL2) << endl;
 	cout << endl;
	cout << "INVARIANTS L " << endl;
	cout << mat.sys.getInvariants(iPrimL.F_()) << endl;	

	cout << "Prim stateR" << endl;
	cout << iPrimR2 << endl;
	cout << "Density" << mat.sys.Density(iPrimR2) << endl;
	cout << endl;
	cout << "INVARIANTS R " << endl;
	cout << mat.sys.getInvariants(iPrimR.F_()) << endl;	

//---------------------
	cout << "Check EOS constants" << endl;
	mat.sys.Eos.checkEosConstants();
#endif
}

void getLimiter(InputSolid inputSolid)
{
  std::string limStr = inputSolid.limiter;
  if(limStr == "superbee") 
  {
    sl = superbee;
  }
  else if(limStr == "minbee") 
  {
    sl = minbee;
  }
  else if(limStr == "vanleer")
  {
    sl = vanleer;
  }
  else if(limStr == "albalda")
  {
    sl = albalda;
  }
  else 
  {
    cerr << "Unknown limiter \'" << sl << "\' - options are superbee, minbee, vanleer, albalda\n";
    exit(1);
  }
}

int main(int argc, char ** argv)
{
  // Enable floating point error checking
  feenableexcept(FE_INVALID);                                                                                                                                                                                    
  feenableexcept(FE_DIVBYZERO);

  char* icStr = argv[1];
  printf("Settings file: %s\n", icStr);
  InputSolid inputSolid;
  inputSolid.readConfigFile(icStr);

  //create and initialize domain
	Domain dom = Domain(inputSolid.cellCountX, 2, inputSolid.xMax); 
  ElasticEOS Eos(inputSolid.matL);
  Material* mat = new Material(inputSolid.matL, dom, Eos);

  //create left and right states and initialise
  ElasticPrimState primStateL(inputSolid.uL, inputSolid.FL, inputSolid.SL);
  ElasticPrimState primStateR(inputSolid.uR, inputSolid.FR, inputSolid.SR);
  double iface = inputSolid.iface;
  ICInterface(*mat, iface, primStateL, primStateR);

  //get limiter
  getLimiter(inputSolid);

  #ifdef debug_
  unitTests(*mat, primStateL, primStateR);
  exit(1);
  #endif
  //solvesystem
  double begin = omp_get_wtime();
  solveSystem(inputSolid, mat);
  double end = omp_get_wtime();
  double elapsed_secs = double(end - begin);
  std::cout << "TIME " << elapsed_secs << std::endl;
  delete mat;
}



