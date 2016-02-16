#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;
#include "ElasticState.h"
#include "ElasticPrimState.h"
#include "System.h"
#include "LevelSet.h"
#include "PlasticModel/PlasticModel.h"
#include "PlasticModel/PerfectPlastic.h"
#include "PlasticModel/NullPlastic.h"
#include "PlasticModel/ZerilliArmstrong.h"
#include "PlasticModel/PathDependent.h"
#include "PlasticModel/SimpleRate.h"
#include "ElasticEOS/MooneyRivlin.h"
#include "ElasticEOS/Romenski.h"
#include "ElasticEOS/Osborne.h"
#include "ElasticEOS/VinetRose.h"
#include "ElasticEOS/MGShock.h"
#include "ElasticEOS/Wilkins.h"
#include "ElasticEOS/StVenantKirchoff.h"
#include <omp.h>
#include <fenv.h>

#define HDT_SOURCEUPDATE
//#define TRANS_SOURCEUPDATE
//#define LIMIT_CONSERVATIVE

enum InterfaceType { Stick, Slip, Weld };
enum ICType { HBLowSpeed, HBHighSpeed, BerylliumShell, Friction, TaylorAnvil, TaylorAnvilTantalum, Radial, BarPlate, SymmetricTaylor, BendingBeam, ExtrapTest };

System::InterpType interpMethod;

const int XDir = 0;
const int YDir = 1;

typedef TinyVector<int,2> IV;
typedef TinyVector<double,2> DV;
typedef TinyVector<double,3> DV3;
typedef TinyVector<DV,2> Rect;

typedef TinyVector<double,23> StateVector;

//! Lagrangian stations
vector<DV> stations;

struct Domain {
  int Ni;
  int Nj;
  int GNi;
  int GNj;
  int GC;
  int starti;
  int startj;
  int endi;
  int endj;
  double Lx;
  double Ly;
  double dx;
  double dy;

  Domain(const int a_Ni, const int a_Nj, const int a_GC,
	 const double a_Lx) :
    Ni(a_Ni), Nj(a_Nj), GNi(a_Ni+2*a_GC), GNj(a_Nj+2*a_GC), GC(a_GC),
    starti(a_GC), startj(a_GC), endi(a_Ni+a_GC), endj(a_Nj+a_GC), 
    Lx(a_Lx), Ly((a_Lx*a_Nj)/a_Ni), dx(a_Lx/a_Ni), dy(a_Lx/a_Ni) {}

  Domain() {}
};

struct Material {
  System sys;
  LevelSet<2> ls;
  blitz::Array<ElasticState, 2> sol;  
  blitz::Array<double, 2> straindot;

  string name;
  double mass0, energy0;

  Material(const string& a_name, const Domain& a_dom, const ElasticEoS* Eos, const PlasticModel* plastic) :
    sys(Eos, plastic), 
    ls(IV(a_dom.Ni, a_dom.Nj), DV(a_dom.Lx, a_dom.Ly), a_dom.GC), 
    sol(a_dom.GNi, a_dom.GNj), straindot(a_dom.GNi, a_dom.GNj), name(a_name) {}
};

typedef vector<Material*>::iterator MaterialIterator;
typedef vector<Material*>::const_iterator MaterialConstIterator;

//! Compute primitive state right of interface for vacuum problem in the x-axis
DV3 VacuumState(const System& sys, ElasticPrimState& pri, const int LR=1)
{
#if 1
  const double eps =1e-7*sys.pressureScale();
  SymmetricMatrix sigma = sys.Stress(pri);
  int count = 0;
  while((fabs(sigma(0,0))>eps || fabs(sigma(0,1))>eps || fabs(sigma(0,2))>eps) && count < 50) 
  {    
    const TinyVector< ElasticPrimState, 3> Re = sys.InterfaceDecompose(pri, sigma, XDir, LR);
    pri -= dot(Re, sigma.col(XDir));
    sigma = sys.Stress(pri);
    ++count;
  }

  if(count>=50) 
  {
    cerr << "Stresses non-zero after 50 iterations in VacuumState!\n";
    cerr << "Final sigma = " << sigma << "\n";
    exit(1);
  }    
  return pri.velocity;
#else
  SymmetricMatrix sigma = sys.Stress(pri);
  const TinyVector< ElasticPrimState, 3> Re = sys.InterfaceDecompose(pri, sigma, XDir, LR);
  const ElasticPrimState p2 = pri-dot(Re, sigma.col(XDir));
  const SymmetricMatrix s2 = sys.Stress(p2);
  const TinyVector< ElasticPrimState, 3> Re2 = sys.InterfaceDecompose(p2, s2, XDir, LR);
  
  pri -= 2.0*dot(Re2, sigma.col(XDir));
  return p2.velocity;
#endif
}

//! Compute states for interface problem in the x-axis
DV3 InterfaceState(const System& sysL, ElasticPrimState& L,
		   const System& sysR, ElasticPrimState& R,
		   const InterfaceType ifaceType)
{
  const ElasticPrimState L0 = L;
  const ElasticPrimState R0 = R;
  const SymmetricMatrix sigmaL0 = sysL.Stress(L);
  const SymmetricMatrix sigmaR0 = sysR.Stress(R);
  SymmetricMatrix sigmaL = sigmaL0;
  SymmetricMatrix sigmaR = sigmaR0;
  
  const double eps =1e-7*sysL.pressureScale();
  const double veps = 1e-4;
  int count = 0;
  while((fabs(sigmaL(0,0)-sigmaR(0,0))>eps || fabs(sigmaL(1,0)-sigmaR(1,0))>eps || fabs(sigmaL(2,0)-sigmaR(2,0))>eps || fabs(L.velocity[0]-R.velocity[0])>veps || ((ifaceType == Stick || ifaceType == Weld)&&(fabs(L.velocity[1]-R.velocity[1])>veps||(fabs(L.velocity[2]-R.velocity[2])>veps)))) && count<50) 
  {    
    // Compute eigenvectors for left and right problem
     const TinyVector< ElasticPrimState, 3> ReL = sysL.InterfaceDecompose(L, sigmaL, XDir,  1);
    const TinyVector< ElasticPrimState, 3> ReR = sysR.InterfaceDecompose(R, sigmaR, XDir, -1);
    
    SymmetricMatrix dvL, dvR;
    for(int i=0; i<3; ++i) 
    {    
      for(int j=i; j<3; ++j) 
      {    
	dvL(i,j) = ReL[i].velocity[j];
	dvR(i,j) = -ReR[i].velocity[j];
      }
    }
    
    DV3 sigma_tilde(0.0);
    
    if(ifaceType == Stick || ifaceType == Weld) {
      
      // Stick stress at interface (match all stresses)
      const SymmetricMatrix denom= dvL + dvR;
      const DV3 dvsigmaL = dvL * sigmaL.col(XDir);  
      const DV3 dvsigmaR = dvR * sigmaR.col(XDir);
      const DV3 dvel = R.velocity - L.velocity;
      sigma_tilde = denom.Inverse() * (dvsigmaL + dvsigmaR + dvel);
      
    } else {
      
      // Slip stress at interface (match only x stress)
      const double denom= dvL(XDir,XDir) + dvR(XDir,XDir);
      const double dvsigmaL = dot(dvL.row(XDir), sigmaL.col(XDir));
      const double dvsigmaR = dot(dvR.row(XDir), sigmaR.col(XDir));
      const double dvel = R.velocity[XDir] - L.velocity[XDir];
      
      sigma_tilde(XDir) = (dvsigmaL + dvsigmaR + dvel) / denom;
      
    }
    
    // Update states to interface state
    const DV3 dsigmaL = sigma_tilde - sigmaL.col(XDir);
    L += dot(ReL, dsigmaL);
    
    const DV3 dsigmaR = sigma_tilde - sigmaR.col(XDir);
    R += dot(ReR, dsigmaR);

    sigmaL = sysL.Stress(L);    
    sigmaR = sysR.Stress(R);    
    ++count;
  }  
  
  if(count>=50) 
  {
    cerr << "Warning: Stresses or velocities don't match after 50 iterations!\n";
//    cerr << "sigmaL = " << sigmaL << "\n";
//    cerr << "sigmaR = " << sigmaR << "\n";
//    cerr << "velL = " << L.velocity << "\n";
//    cerr << "velR = " << R.velocity << "\n";
    cerr << "L0 = " << L0 << "\n";
    cerr << "R0 = " << R0 << "\n";
  }    
  
  const DV3 vel = 0.5*(L.velocity+R.velocity);
  return vel;
}

//! Compute the ghost cell state for a material-vacuum interface
DV CalcGhostStateVacuum(const System& sys, ElasticState& cns, const DV& normal) 
{
  // Rotate state to align with x-axis
  const SquareMatrix R = rotateToX(normal);
  const SquareMatrix RT= R.Transpose();
    
  ElasticPrimState ghost = sys.ToPrimitive(Rotate(cns, R, RT));

  // Calculate Vacuum Star State 
  DV3 vel = VacuumState(sys, ghost);

  // Rotate the state back
  ghost.Rotate(RT, R);
  vel = RT*vel;

  cns = sys.ToConservative(ghost);  

  return DV(vel[0], vel[1]);
}

//! Compute ghost cell state determined by an interface problem between two materials
double TestImpact(const System& sysL, const System& sysR, 
		  const ElasticState &cnsL, const ElasticState &cnsR,
		  const DV& normal)
{
  // Rotate states to align with x-axis
  const SquareMatrix R = rotateToX(normal);
  const SquareMatrix RT= R.Transpose();
  ElasticPrimState ghostL = sysL.ToPrimitive(Rotate(cnsL, R, RT));
  ElasticPrimState ghostR = sysR.ToPrimitive(Rotate(cnsR, R, RT));

  // Calculate free surface velocities
  const DV3 velL = VacuumState(sysL, ghostL);
  const DV3 velR = VacuumState(sysR, ghostR, -1);
  return velL[XDir]-velR[XDir];
}

//! Compute ghost cell states determined by an interface problem between two materials
DV CalcGhostStateInteraction(const System& sysL, const System& sysR, 
			       ElasticState &cnsL, const ElasticState &cnsR,
			       const DV& normal, const InterfaceType ifaceType)
{
  // Rotate state to align with x-axis
  const SquareMatrix R = rotateToX(normal);
  const SquareMatrix RT= R.Transpose();

  ElasticPrimState ghostL = sysL.ToPrimitive(Rotate(cnsL,R,RT));
  ElasticPrimState ghostR = sysR.ToPrimitive(Rotate(cnsR,R,RT));

  // Calculate Interface States
  DV3 vel = InterfaceState(sysL, ghostL, sysR, ghostR, ifaceType);

  // Rotate the states back
  ghostL.Rotate(RT, R);
  vel = RT*vel;

  cnsL = sysL.ToConservative(ghostL);

  return DV(vel[0], vel[1]);
}

//! Evaluate the spatial divergence of the elastic deformation gradient.
TinyVector<double, 2> divRhoFe(const Material& mat, const Domain& dom, const IV& pos) {

  if(!mat.ls.isInterior(pos)) {
    return TinyVector<double,2>(0.);
  }
  
  static const IV ioff(1,0);
  static const IV joff(0,1);

  const IV posim = pos-ioff;
  const IV posip = pos+ioff;  
  const IV posjm = pos-joff;
  const IV posjp = pos+joff;
  
  int di=2, dj=2;
  SquareMatrix rhoFW, rhoFE, rhoFN, rhoFS;
  if(mat.ls.isInterior(posim)) {
    rhoFW = mat.sol(posim).rhoFe;
  } else {
    --di;
    rhoFW = mat.sol(pos).rhoFe;
  }
  if(mat.ls.isInterior(posip)) {
    rhoFE = mat.sol(posip).rhoFe;
  } else {
    --di;
    rhoFE = mat.sol(pos).rhoFe;
  }
  if(mat.ls.isInterior(posjp)) {
    rhoFN = mat.sol(posjp).rhoFe;
  } else {
    --dj;
    rhoFN = mat.sol(pos).rhoFe;
  }
  if(mat.ls.isInterior(posjm)) {
    rhoFS = mat.sol(posjm).rhoFe;
  } else {
    --dj;
    rhoFS = mat.sol(pos).rhoFe;
  }
  

  TinyVector<double, 2> div;
  for(int k=0; k<2; ++k) {
    div[k] = 0.0;
    if(di>0) div[k] += (rhoFE(0,k)-rhoFW(0,k))/(di*dom.dx);
    if(dj>0) div[k] += (rhoFN(1,k)-rhoFS(1,k))/(dj*dom.dy);
  }
  
  return div;
}

//! Set the domain boundary ghost cells in the x direction
void xBCs(const Material& mat, const Domain& dom, blitz::Array<ElasticState, 2>& sol) {  
  
  /* i boundaries */
#pragma omp parallel for schedule(dynamic)
  for(int i=0; i<dom.starti; ++i) {
    const int imagei = 2*dom.starti-1-i;
    for(int j=dom.startj; j<dom.endj; ++j) {
      sol(i,j) = sol(imagei,j);
      sol(i,j).Reflect(XDir);
    }    
  }
  
#pragma omp parallel for schedule(dynamic)
  for(int i=dom.endi; i<dom.GNi; ++i) {
    const int imagei = 2*dom.endi-1-i;
    for(int j=dom.startj; j<dom.endj; ++j) {
      sol(i,j) = sol(imagei,j);
      sol(i,j).Reflect(XDir);
    }  
  }
}

void yBCs(const Material& mat, const Domain& dom, blitz::Array<ElasticState, 2>& sol) {  
  
  /* j boundaries */
#pragma omp parallel for schedule(dynamic)
  for(int j=0; j<dom.startj; ++j) {
    const int imagej = 2*dom.startj-1-j;
    for(int i=dom.starti; i<dom.endi; ++i) {
      sol(i,j) = sol(i,imagej);
      sol(i,j).Reflect(YDir);
    }      
  }    
  
#pragma omp parallel for schedule(dynamic)
  for(int j=dom.endj; j<dom.GNj; ++j) {
    const int imagej = 2*dom.endj-1-j;
    for(int i=dom.starti; i<dom.endi; ++i) {
      sol(i,j) = sol(i,imagej);
      sol(i,j).Reflect(YDir);
    }    
  }
}

//! Set the domain boundary ghost cells in both x and y direction
void BCs(Material& mat, const Domain& dom) {  
  xBCs(mat,dom,mat.sol);
  yBCs(mat,dom,mat.sol);
}

//! Compute the solution update due to fluxes in the x-direction
void solveX(blitz::Array<ElasticState, 2>& fluxX, blitz::Array<ElasticState, 2>& source, blitz::Array<ElasticState, 2>& sol, const Material& mat, const Domain& dom, const double dt) {
  xBCs(mat, dom, sol);

  const double dt_dx = dt/dom.dx;
  const int si = dom.starti;
  const int sjC = dom.startj;
  const int ei = dom.endi;
  const int ejC = dom.endj;

#pragma omp parallel for schedule(dynamic)
  for(int jC=sjC; jC<ejC; ++jC) {

    vector<bool> interior(dom.GNi);

    // Cell primitive states
    vector<ElasticPrimState> prim(dom.GNi-1);

    // Left and right reconstructed states
    vector<ElasticPrimState> left(dom.GNi-1);
    vector<ElasticPrimState> right(dom.GNi-1);

    // Godunov states
    vector<ElasticPrimState> god(dom.GNi-1);

    for(int i=si; i<ei; ++i) {
      interior[i] = mat.ls.isInterior(IV(i,jC));
    }
    interior[si-1] = false;
    interior[si-2] = false;
    interior[si-3] = false;
    interior[ei] =   false;
    interior[ei+1] = false;
    interior[ei+2] = false;

    // Compute primitive states
    for(int i=si-2; i<=ei+1; ++i) {
      if((i<ei&&interior[i+2]) || interior[i+1] || interior[i] || interior[i-1] || (i>=si&&interior[i-2])) 
      {	
	prim[i] = mat.sys.ToPrimitive(sol(i,jC));
      }      
    }
    
    // Limit slopes and compute reconstructed states
    for(int i=si-1; i<ei+1; ++i) {
      if(interior[i-1]||interior[i]||interior[i+1]) 
      {
#ifdef LIMIT_CONSERVATIVE
	const ElasticState slope = mat.sys.ComputeLimitedSlope(sol(i-1,jC), sol(i,jC), sol(i+1,jC), XDir);
	ElasticState Cleft = sol(i,jC)-0.5*slope;
	ElasticState Cright = sol(i,jC)+0.5*slope;
	const ElasticState fluxDiff = mat.sys.Flux(Cleft,XDir) - mat.sys.Flux(Cright,XDir);
#else
	const ElasticPrimState slope = mat.sys.ComputeLimitedSlope(prim[i-1], prim[i], prim[i+1], XDir, dt_dx);
	const ElasticPrimState Pright = prim[i]+0.5*slope;	
	const ElasticPrimState Pleft  = prim[i]-0.5*slope;
	
	const ElasticState fluxDiff = mat.sys.Flux(Pleft,XDir) - mat.sys.Flux(Pright,XDir);

	// Half-dt update
	ElasticState Cleft = mat.sys.ToConservative(Pleft);
	ElasticState Cright = mat.sys.ToConservative(Pright);
#endif
	
	const ElasticState update = (0.5*dt_dx)*fluxDiff;
	Cleft += update;
	Cright += update;	
	left[i] = mat.sys.ToPrimitive(Cleft);
	right[i] = mat.sys.ToPrimitive(Cright);
      }      
    }      

    // Compute interface Godunov states and fluxes
    for(int i=si; i<=ei; ++i) {
      if(interior[i-1]||interior[i]) {	
	god[i] = mat.sys.GodunovState(right[i-1], left[i], XDir);
	fluxX(i,jC) = mat.sys.Flux(god[i], XDir);
      }      
    }

    // Compute stability sources
    for(int i=si; i<ei; ++i) {
      if(interior[i])
	{		  
#ifdef HDT_SOURCEUPDATE
	  const ElasticState hC = sol(i,jC) - 0.5 * dt_dx * (fluxX(i+1,jC)-fluxX(i,jC));      
	  const ElasticPrimState pC = mat.sys.ToPrimitive(hC);
#else 
	  const ElasticPrimState pC = prim[i];
#endif
	  source(i,jC) += dt_dx*mat.sys.Source(god[i], god[i+1], pC, XDir);
	}      
    }    
  }
}

//! Compute the first order solution update due to fluxes in the x-direction
void solveXfirst(blitz::Array<ElasticState, 2>& sol, blitz::Array<ElasticState, 2>& res, const Material& mat, const Domain& dom, const double dt, const double dt_dx) {
  xBCs(mat, dom, sol);
  const int si = dom.starti;
  const int sjC = dom.startj;
  const int ei = dom.endi;
  const int ejC = dom.endj;

#pragma omp parallel for schedule(dynamic)
  for(int jC=sjC; jC<ejC; ++jC) {

    // Cell states
    vector<bool> interior(dom.GNi);
    vector<ElasticPrimState> prim(dom.GNi-1);
    vector<ElasticState> pfluxX(dom.GNi-1);
    vector<ElasticPrimState> god(dom.GNi-1);

    for(int i=si; i<ei; ++i) {
      interior[i] = ((  jC-2>=sjC && mat.ls.isInterior(IV(i,jC-2)))
		     ||(jC-1>=sjC && mat.ls.isInterior(IV(i,jC-1)))
		     ||(mat.ls.isInterior(IV(i,jC)))
		     ||(jC+1<ejC && mat.ls.isInterior(IV(i,jC+1)))
		     ||(jC+2<ejC && mat.ls.isInterior(IV(i,jC+2))));
    }
    interior[si-1] = false;
    interior[si-2] = false;
    interior[ei] = false;
    interior[ei+1] = false;

    // Compute primitive states
    for(int i=si-1; i<=ei; ++i) {
      if(interior[i-1]||interior[i]||interior[i+1]) {	
	prim[i] = mat.sys.ToPrimitive(sol(i,jC));
      }
    }
    
    // Compute interface Godunov states and fluxes
    for(int i=si; i<=ei; ++i) {
      if(interior[i-1]||interior[i]) {	
	god[i] = mat.sys.GodunovState(prim[i-1], prim[i], XDir);
	pfluxX[i] = mat.sys.Flux(god[i], XDir);
      }      
    }

    // Update, including sources
    for(int i=si; i<ei; ++i) {
      if(interior[i]) {		  
	res(i,jC) = sol(i,jC) - dt_dx * (pfluxX[i+1]-pfluxX[i]) ; 
#ifdef TRANS_SOURCEUPDATE
	res(i,jC) += dt_dx*mat.sys.Source(god[i], god[i+1], prim[i], XDir);
	//	res(i,jC) += dt * mat.sys.GeometricSource(mat.ls.getCellPos(IV(i,jC)), sol(i,jC));
#endif
      }      
    }    
  }
}

//! Compute the first order solution update due to fluxes in the y-direction
void solveYfirst(blitz::Array<ElasticState, 2>& sol, blitz::Array<ElasticState, 2>& res, const Material& mat, const Domain& dom, const double dt, const double dt_dy) {
  yBCs(mat, dom, sol);
  const int siC = dom.starti;
  const int sj = dom.startj;
  const int eiC = dom.endi;
  const int ej = dom.endj;

#pragma omp parallel for schedule(dynamic)
  for(int iC=siC; iC<eiC; ++iC) {

    // Cell states
    vector<bool> interior(dom.GNj);
    vector<ElasticPrimState> prim(dom.GNj-1);
    vector<ElasticState> pfluxY(dom.GNj-1);
    vector<ElasticPrimState> god(dom.GNj-1);

    for(int j=sj; j<ej; ++j) {
      interior[j] = ((  iC-2>=siC && mat.ls.isInterior(IV(iC-2,j)))
		     ||(iC-1>=siC && mat.ls.isInterior(IV(iC-1,j)))
		     ||(mat.ls.isInterior(IV(iC,j)))
		     ||(iC+1<eiC && mat.ls.isInterior(IV(iC+1,j)))
		     ||(iC+2<eiC && mat.ls.isInterior(IV(iC+2,j))));
    }
    interior[sj-1] = false;
    interior[sj-2] = false;
    interior[ej]   = false;
    interior[ej+1] = false;
    
    // Compute primitive states
    for(int j=sj-1; j<=ej; ++j) {
      if(interior[j-1]||interior[j]||interior[j+1]) {	
	prim[j] = mat.sys.ToPrimitive(sol(iC,j));
      }
    }
    
    // Compute interface Godunov states and fluxes
    for(int j=sj; j<=ej; ++j) {
      if(interior[j-1]||interior[j]) {	
	god[j] = mat.sys.GodunovState(prim[j-1], prim[j], YDir);
	pfluxY[j] = mat.sys.Flux(god[j], YDir);
      }      
    }

    // Update, including stability sources
    for(int j=sj; j<ej; ++j) {
      if(interior[j])
	{		  
	  res(iC,j) = sol(iC,j) - dt_dy * (pfluxY[j+1]-pfluxY[j]) ; 
#ifdef TRANS_SOURCEUPDATE
	  res(iC,j) += dt_dy*mat.sys.Source(god[j], god[j+1], prim[j], YDir);
	  //	  res(iC,j) += dt * mat.sys.GeometricSource(mat.ls.getCellPos(IV(iC,j)), sol(iC,j));
#endif
	}      
    }    
  }
}

//! Compute the solution update due to fluxes in the y-direction
void solveY(blitz::Array<ElasticState, 2>& fluxY, blitz::Array<ElasticState, 2>& source, blitz::Array<ElasticState, 2> &sol, const Material& mat, const Domain& dom, const double dt) {
  yBCs(mat, dom, sol);

  const double dt_dx = dt/dom.dy;
  const int si = dom.starti;
  const int sj = dom.startj;
  const int ei = dom.endi;
  const int ej = dom.endj;

#pragma omp parallel for schedule(dynamic)
  for(int iC=si; iC<ei; ++iC) {

    // Cell states
    vector<bool> interior(dom.GNj);

    // Slopes
    vector<ElasticPrimState> prim(dom.GNj);

    // Interface Godunov states
    vector<ElasticPrimState> god(dom.GNj-1);

    // Left and right reconstructed states
    vector<ElasticPrimState> left(dom.GNj-1);
    vector<ElasticPrimState> right(dom.GNj-1);

    for(int j=sj; j<ej; ++j) {
      interior[j] = mat.ls.isInterior(IV(iC,j));
    }    
    interior[sj-1] = false;
    interior[sj-2] = false;
    interior[sj-3] = false;
    interior[ej] = false;
    interior[ej+1] = false;
    interior[ej+2] = false;

    // Compute primitive states
    for(int j=sj-2; j<=ej+1; ++j) {
      if((j<ej&&interior[j+2]) || interior[j+1] || interior[j] || interior[j-1] || (j>=sj&&interior[j-2])) 
      {	
	prim[j] = mat.sys.ToPrimitive(sol(iC,j));
      }      
    }
    
    // Limit slopes and compute reconstructed states
    for(int j=sj-1; j<ej+1; ++j) {
      if(interior[j-1]||interior[j]||interior[j+1]) 
      {
#ifdef LIMIT_CONSERVATIVE
	const ElasticState slope = mat.sys.ComputeLimitedSlope(sol(iC,j-1), sol(iC,j), sol(iC,j+1), YDir);
	ElasticState Cleft = sol(iC,j)-0.5*slope;
	ElasticState Cright = sol(iC,j)+0.5*slope;
	const ElasticState fluxDiff = mat.sys.Flux(Cleft,YDir) - mat.sys.Flux(Cright,YDir);
#else
	const ElasticPrimState slope = mat.sys.ComputeLimitedSlope(prim[j-1], prim[j], prim[j+1], YDir, dt_dx);
	const ElasticPrimState Pright = prim[j]+0.5*slope;	
	const ElasticPrimState Pleft  = prim[j]-0.5*slope;
	const ElasticState fluxDiff = mat.sys.Flux(Pleft,YDir) - mat.sys.Flux(Pright,YDir);
	
	// Half-dt update
	ElasticState Cleft = mat.sys.ToConservative(Pleft);
	ElasticState Cright = mat.sys.ToConservative(Pright);
#endif

	const ElasticState update = (0.5*dt_dx)*fluxDiff;
	Cleft += update;
	Cright += update;	
	left[j] = mat.sys.ToPrimitive(Cleft);
	right[j] = mat.sys.ToPrimitive(Cright);
      }      
    }      
    
    // Compute interface Godunov states and fluxes
    for(int j=sj; j<=ej; ++j) {
      if(interior[j-1]||interior[j]) {	
	god[j] = mat.sys.GodunovState(right[j-1],left[j],YDir);
	fluxY(iC,j) = mat.sys.Flux(god[j], YDir);
      }      
    }
    
    // Compute stability sources
    for(int j=sj; j<ej; ++j) {
      if(interior[j]) {	
#ifdef HDT_SOURCEUPDATE
	const ElasticState hC = sol(iC,j) - 0.5 * dt_dx * (fluxY(iC,j+1)-fluxY(iC,j));      
	const ElasticPrimState pC = mat.sys.ToPrimitive(hC);
#else 
	const ElasticPrimState pC = prim[j];
#endif
	source(iC,j) += dt_dx*mat.sys.Source(god[j], god[j+1], pC, YDir);
      }      
    }    
  }
}

void solveXunsplit(blitz::Array<ElasticState, 2>& fluxX, blitz::Array<ElasticState, 2>& source, blitz::Array<ElasticState, 2>& sol, const Material& mat, const Domain& dom, const double dt) {
  const double dt_dy = dt/dom.dy;

  blitz::Array<ElasticState, 2> psol(dom.GNi, dom.GNj);
  solveYfirst(sol, psol, mat, dom, 0.5*dt, 0.5*dt_dy);
  solveX(fluxX, source, psol, mat, dom, dt);
}

void solveYunsplit(blitz::Array<ElasticState, 2>& fluxY, blitz::Array<ElasticState, 2>& source, blitz::Array<ElasticState, 2>& sol, const Material& mat, const Domain& dom, const double dt) {
  const double dt_dx = dt/dom.dx;

  blitz::Array<ElasticState, 2> psol(dom.GNi, dom.GNj);
  solveXfirst(sol, psol, mat, dom, 0.5*dt, 0.5*dt_dx);
  solveY(fluxY, source, psol, mat, dom, dt);
}

// Compute geometric source increment for all interior cells
void computeGeometricSource(blitz::Array<ElasticState, 2>& source, const blitz::Array<ElasticState, 2>& sol, const Material& mat, const Domain& dom, const double dt) {
#pragma omp parallel for schedule(dynamic)
  for(int j=dom.startj; j<dom.endj; ++j) {
    for(int i=dom.starti; i<dom.endi; ++i) {
      const IV pos(i,j);
      if(mat.ls.isInterior(pos)) {
	source(i,j) += dt * mat.sys.GeometricSource(mat.ls.getCellPos(pos), sol(i,j));
      }
    }
  }
}


void elasticSolve(Material& mat, const Domain& dom, const double dt, const bool axisym) {

  const double dt_dx = dt/dom.dx;
  const double dt_dy = dt/dom.dy;

  // Compute raw fluxes & sources
  blitz::Array<ElasticState, 2> source(dom.GNi, dom.GNj);
  source = ElasticState::Zero;  

  blitz::Array<ElasticState, 2> fluxY(dom.GNi, dom.GNj);  
  solveYunsplit(fluxY, source, mat.sol, mat, dom, dt);

  blitz::Array<ElasticState, 2> fluxX(dom.GNi, dom.GNj);  
  solveXunsplit(fluxX, source, mat.sol, mat, dom, dt);

  if(axisym) {
    blitz::Array<ElasticState, 2> sol(dom.GNi, dom.GNj);
#pragma omp parallel for schedule(dynamic)
    for(int j=dom.startj; j<dom.endj; ++j) {
      for(int i=dom.starti; i<dom.endi; ++i) {
	if(mat.ls.isInterior(IV(i,j)))
	  {	
	    const ElasticState FN = fluxY(i,j+1);
	    const ElasticState FS = fluxY(i,j);
	    const ElasticState FE = fluxX(i+1,j);
	    const ElasticState FW = fluxX(i,j);
	    sol(i,j) = mat.sol(i,j)-0.5*(dt_dx * (FE-FW) + dt_dy * (FN-FS) - source(i,j));
	  }      
      }
    }
    computeGeometricSource(source, sol, mat, dom, dt);
  }

#pragma omp parallel for schedule(dynamic)
  for(int j=dom.startj; j<dom.endj; ++j) {
    for(int i=dom.starti; i<dom.endi; ++i) {
      if(mat.ls.isInterior(IV(i,j)))
	{	
	  const ElasticState FN = fluxY(i,j+1);
	  const ElasticState FS = fluxY(i,j);
	  const ElasticState FE = fluxX(i+1,j);
	  const ElasticState FW = fluxX(i,j);
	  const SymmetricMatrix eps0 = mat.sys.AlmansiStrain(mat.sol(i,j)).Deviator();
	  mat.sol(i,j) -= dt_dx * (FE-FW) + dt_dy * (FN-FS) - source(i,j);
	  const SymmetricMatrix eps1 = mat.sys.AlmansiStrain(mat.sol(i,j)).Deviator();
	  mat.straindot(i,j) = sqrt(2./3.)*(eps1-eps0).Schur()/dt;
	}      
    }
  }
}

//! Compute level set value at specified position for a circle
double lsball(const DV& pos, const DV& centre, const double radius)
{
  const TinyVector<double, 2> dist = pos-centre;
  const double r = sqrt(dot(dist,dist));
  return r-radius;
}

//! Compute level set value at specified position for a grid-aligned rectangle
double lsrect(const DV& pos, const Rect& rect)
{
  int loc[2];
  int tot=0;
  for(int i=0; i<2; ++i) 
    {
      if(pos[i]<rect[0][i]) {loc[i] = 0; ++tot;}
      else if(pos[i]>rect[1][i]) {loc[i] = 1; ++tot;}  
      else loc[i] = 2;
    }
  
  if(tot==0) 
    { // Inside
      return -min(min(pos[0]-rect[0][0], rect[1][0]-pos[0]),
		  min(pos[1]-rect[0][1], rect[1][1]-pos[1]));
    } else if(tot==1) {
    if(loc[0]!=2) 
      {
	return fabs(rect[loc[0]][0]-pos[0]);
      }
    else 
      {
	return fabs(rect[loc[1]][1]-pos[1]);
      }
  } else 
    {
      const double ldx = rect[loc[0]][0]-pos[0];
      const double ldy = rect[loc[1]][1]-pos[1];
      return sqrt(ldx*ldx+ldy*ldy);
    }        
}

double mag(const TinyVector<double, 2>& V) {
  return sqrt(dot(V,V));
}

double lsline(const DV& pos, const DV& p1, const DV& p2) {

  // Get distance along line
  const DV V = p2-p1;
  const DV delta = pos-p1;
  const double d = dot(delta,V);
  if(d<0.0) {
    return mag(delta);
  } else if(d>dot(V,V)) {
    return mag(pos-p2);
  } else {
    DV N(-V[1],V[0]);
    N /= mag(N);
    return dot(N,delta);
  }
}

bool compdist(double A, double B) {
  return fabs(A)<fabs(B);
}

//! Compute level set value at specified position for a grid-aligned rectangle
double lsquad(const DV& pos, const DV& p1, const DV& p2, const DV& p3, const DV& p4)
{
  vector<double > dist;

  dist.push_back(lsline(pos, p1,p2));
  dist.push_back(lsline(pos, p2,p3));
  dist.push_back(lsline(pos, p3,p4));
  dist.push_back(lsline(pos, p4,p1));

  sort(dist.begin(), dist.end(), compdist);
  
  return dist[0];
}

//! Set block initial condition travelling at specified velocity
void ICBlock(Material& mat, 
	     const Domain& dom,
	     const Rect& rect,
	     const DV3& vel) {

  ElasticPrimState pri(SquareMatrix(1,0,0,0,1,0,0,0,1),
			      vel, 0, 0);

  const ElasticState cns = mat.sys.ToConservative(pri);  

  mat.ls.reset();
  for(int j=dom.startj-3; j<dom.endj+3; ++j) {
    for(int i=dom.starti-3; i<dom.endi+3; ++i) {
      const LevelSet<2>::Coord pos = mat.ls.getCellPos(IV(i,j));
      mat.ls.setPhi(IV(i,j), lsrect(pos, rect));
    }    
  }
  mat.ls.reinitialize();

  // Create the initial solution vector
  for(int j=dom.startj; j<dom.endj; ++j) {
    for(int i=dom.starti; i<dom.endi; ++i) {
      if(mat.ls.isInterior(IV(i,j))) {
	mat.sol(i,j) = cns;
      }      
    }    
  }
}

//! Set bending beam initial condition
void ICBeam(Material& mat, 
	     const Domain& dom,
	     const Rect& rect) {

  ElasticPrimState pri(SquareMatrix(1,0,0,0,1,0,0,0,1),
		       DV3(0.), 0, 0);


  mat.ls.reset();
  for(int j=dom.startj-3; j<dom.endj+3; ++j) {
    for(int i=dom.starti-3; i<dom.endi+3; ++i) {
      const LevelSet<2>::Coord pos = mat.ls.getCellPos(IV(i,j));
      mat.ls.setPhi(IV(i,j), lsrect(pos, rect));
    }    
  }
  mat.ls.reinitialize();

  // Create the initial solution vector
  for(int j=dom.startj; j<dom.endj; ++j) {
    for(int i=dom.starti; i<dom.endi; ++i) {
      if(mat.ls.isInterior(IV(i,j))) {
	const LevelSet<2>::Coord pos = mat.ls.getCellPos(IV(i,j));
	const double x=pos[0]-3.5;
	const double A=0.004337*0.235974*56.6369;
	const double B=0.004337*0.235974*57.6355;
	const double Omega1=0.7883902;
	const double uy=A*(sinh(Omega1*(x+3))+sin(Omega1*(x+3)))-
	                B*(cosh(Omega1*(x+3))+cos(Omega1*(x+3)));
	pri.velocity[1]=uy;
	const ElasticState cns = mat.sys.ToConservative(pri);  
	mat.sol(i,j) = cns;
      }      
    }    
  }
}

//! Set block initial condition travelling at specified velocity
void ICQuad(Material& mat, 
	    const Domain& dom,
	    const DV& p1,
	    const DV& p2,
	    const DV& p3,
	    const DV& p4,
	    const DV3& vel) {

  ElasticPrimState pri(SquareMatrix(1,0,0,0,1,0,0,0,1),
			      vel, 1e-4, 0);

  const ElasticState cns = mat.sys.ToConservative(pri);  

  mat.ls.reset();
  for(int j=dom.startj-3; j<dom.endj+3; ++j) {
    for(int i=dom.starti-3; i<dom.endi+3; ++i) {
      const LevelSet<2>::Coord pos = mat.ls.getCellPos(IV(i,j));
      mat.ls.setPhi(IV(i,j), lsquad(pos, p1,p2,p3,p4));
    }    
  }
  mat.ls.reinitialize();

  // Create the initial solution vector
  for(int j=dom.startj; j<dom.endj; ++j) {
    for(int i=dom.starti; i<dom.endi; ++i) {
      if(mat.ls.isInterior(IV(i,j))) {
	mat.sol(i,j) = cns;
      }      
    }    
  }
}

//! Set block initial condition travelling at specified velocity
void ICBall(Material& mat,
	    const Domain& dom,
	    const DV& centre,
	    const double radius,
	    const DV3& vel) {
  
  ElasticPrimState pri(SquareMatrix(1,0,0,0,1,0,0,0,1),			      
			      vel, 1e-4, 0);

  const ElasticState cns = mat.sys.ToConservative(pri);  

  mat.ls.reset();
  for(int j=dom.startj-3; j<dom.endj+3; ++j) {
    for(int i=dom.starti-3; i<dom.endi+3; ++i) {
      const LevelSet<2>::Coord pos = mat.ls.getCellPos(IV(i,j));
      mat.ls.setPhi(IV(i,j), lsball(pos, centre, radius));
    }    
  }
  mat.ls.reinitialize();

  // Create the initial solution vector
  for(int j=dom.startj; j<dom.endj; ++j) {
    for(int i=dom.starti; i<dom.endi; ++i) {
      if(mat.ls.isInterior(IV(i,j))) {
	mat.sol(i,j) = cns;
      }      
    }    
  }
}

//! Set Zaleszak's circle initial condition for extrapolation test
void ICZaleszak(Material& mat,
		const Domain& dom,
		const DV& centre,
		const double radius,
		const double width) {
  
  ElasticPrimState pri(SquareMatrix(1,0,0,0,1,0,0,0,1),			      
		       DV3(0,0,0), 1e-4, 0.0);
  DV bl = centre;
  bl[0]-=width/2;
  DV tr = centre;
  tr[0]+=width/2;
  tr[1]+=radius*2;

  const Rect rect(bl, tr);

  mat.ls.reset();
  for(int j=dom.startj-3; j<dom.endj+3; ++j) {
    for(int i=dom.starti-3; i<dom.endi+3; ++i) {
      const LevelSet<2>::Coord pos = mat.ls.getCellPos(IV(i,j));
      const double lsb = lsball(pos, centre, radius);
      const double lsr = -lsrect(pos, rect);
 
     if(lsb<0. && lsr<0.) 
      {	
	mat.ls.setPhi(IV(i,j), std::max(lsb,lsr));
      }
      else if(lsb>0.0 && lsr>0.0)
      {
	mat.ls.setPhi(IV(i,j), std::min(lsb, lsr));
      }     
      else if(lsb>=0.0)
      {
	mat.ls.setPhi(IV(i,j), lsb);
      }
      else
      {
	mat.ls.setPhi(IV(i,j), lsr);
      }
      
    }    
  }
  mat.ls.reinitialize();

  // Create the initial solution vector
  for(int j=dom.startj; j<dom.endj; ++j) {
    for(int i=dom.starti; i<dom.endi; ++i) {
      const LevelSet<2>::Coord pos = mat.ls.getCellPos(IV(i,j)) - centre;
      const double alpha = atan2(pos[1],pos[0]);
      pri.Fe(0,1) = 0.2*sin(alpha);
      pri.Fe(1,0) = 0.2*cos(alpha);
      
      const ElasticState cns = mat.sys.ToConservative(pri);  
      if(mat.ls.isInterior(IV(i,j))) {
	mat.sol(i,j) = cns;
      }      
    }    
  }
}

double shellF(const double alpha, const double lambda) 
{
  // Integrate F by midpoint rule
  const int res = 1000;
  const double dx = (1.0-lambda)/res;
  double igral = 0.0;
  for(int i=0; i<res; ++i) 
    {
      const double x = lambda + dx*(i+0.5);
      igral += dx*x*log(1.0 + ((2.0+alpha)*alpha)/(x*x));
    }
  return igral;
}

//! Set up a stopping shell initial condition of the specified inner and 
//! outer radius and inner stopping radius
void ICShell(Material& mat,
	     const Domain& dom,
	     const DV& centre,
	     const double R0,
	     const double R1,
	     const double Rstop) 
{
  ElasticPrimState pri(SquareMatrix(1,0,0,0,1,0,0,0,1),			      
			      TinyVector<double, 3>(0, 0, 0),
			      0, 0);

  const ElasticState cns = mat.sys.ToConservative(pri);  
  const double rho0 = mat.sys.m_eos->rho0();
  const double Y0 = mat.sys.m_plastic->yieldStress(0.0, 0.0, 0.0);

  const double alpha = (R1-R0) / R0;
  const double lambda = Rstop / R0;

  const double u0 = sqrt(2.0 * Y0 * shellF(alpha, lambda) 
			 / (sqrt(3.0)*rho0*log(R1/R0)));
  cout << "Y0 = " << Y0 << endl;
  cout << "r'0 / R0 = " << Rstop / R0 << endl;
  cout << "u0 = " << u0 << " m/s = " << u0 / 1E4 << " cm/mus " << endl;
  cout << "alpha = " << alpha << endl;
  cout << "lambda = " << lambda << endl;
  cout << "F = " << shellF(alpha,lambda) << endl;

  mat.ls.reset();
  for(int j=dom.startj; j<dom.endj; ++j) {
    for(int i=dom.starti; i<dom.endi; ++i) {
      const IV idx(i,j);
      const DV pos = mat.ls.getCellPos(idx) - centre;
      const double r = sqrt(dot(pos,pos));
      const double d1 = r-R1;
      const double d0 = R0-r;
      double d=d1;
      if(fabs(d0)<fabs(d)) d=d0;
      mat.ls.setPhi(idx, d);
    }    
  }

  mat.ls.reinitialize();

  // Create the initial solution vector
  for(int j=dom.startj; j<dom.endj; ++j) {
    for(int i=dom.starti; i<dom.endi; ++i) {
      const IV idx(i,j);
      if(mat.ls.isInterior(idx)) {	
	const DV pos = mat.ls.getCellPos(idx) - centre;
	ElasticPrimState priC = pri;
	const double mag = dot(pos,pos);
	if(mag>0.0) 
	  {	  	  
	    const double r = sqrt(mag);
	    const DV norm = pos / r;
	    const double ur = u0 * R0 / r;
	    priC.velocity[0] = -norm[0] * ur;
	    priC.velocity[1] = -norm[1] * ur;
	  }
	mat.sol(i,j) = mat.sys.ToConservative(priC);  
      }      
    }    
  }
}

//! Set up a shell initial condition of the specified inner and 
//! outer radius and inner velocity
void ICShell2(Material& mat,
	      const Domain& dom,
	      const DV& centre,
	      const double R0,
	      const double R1,
	      const double vel) 
{

  ElasticPrimState pri(SquareMatrix(1,0,0,0,1,0,0,0,1),
			      TinyVector<double, 3>(0, 0, 0),
			      1e-4, 0);

  const ElasticState cns = mat.sys.ToConservative(pri);  

  mat.ls.reset();
  for(int j=dom.startj; j<dom.endj; ++j) {
    for(int i=dom.starti; i<dom.endi; ++i) {
      const IV idx(i,j);
      const DV pos = mat.ls.getCellPos(idx) - centre;
      const double r = sqrt(dot(pos,pos));
      const double d1 = r-R1;
      const double d0 = R0-r;
      double d=d1;
      if(fabs(d0)<fabs(d)) d=d0;
      mat.ls.setPhi(idx, d);
    }    
  }
  
  mat.ls.reinitialize();

  // Create the initial solution vector
  for(int j=0; j<dom.GNj; ++j) {
    for(int i=0; i<dom.GNi; ++i) {
      const IV idx(i,j);
      if(mat.ls.isInterior(idx)) {	
	const DV pos = mat.ls.getCellPos(idx) - centre;
	const double r = sqrt(dot(pos,pos));
	ElasticPrimState priC = pri;
	const DV norm = pos / r;
	const double ur = vel * R0 / r;
	priC.velocity[0] = -norm[0] * ur;
	priC.velocity[1] = -norm[1] * ur;
	mat.sol(i,j) = mat.sys.ToConservative(priC);  
      }      
    }    
  }
}

//! Set ghost cell values
void setGhosts(vector<Material*>& mats, const Domain& dom, const double dt, const InterfaceType ifaceType)
{
  // Extrapolate interior state to adjacent cells
  for(MaterialIterator mit = mats.begin(); mit!=mats.end(); ++mit) {
    Material& mat = *(*mit);
    //  mat.ls.extrapolateEntropyToBoundary(mat.sol, mat.sys);    
    BCs(mat, dom);
    mat.ls.extrapolateToAdjacent(mat.sol, mat.sys);
    mat.ls.extrapolateToAdjacent(mat.straindot, mat.sys);
  }

  // TO CONSIDER:
  // when >2 materials are interacting in a cell, the following might need changing
  // e.g. we need to set consistent level set values for all three materials...
  // also, in this case, the ghost states found should be independent of the ordering in which we solve the interactions
  // (which they're not at the moment)

  // Look for boundary interactions between each pair of materials
  for(MaterialIterator mit = mats.begin(); mit!=mats.end(); ++mit) {
    Material& mat = *(*mit);

    const vector<IV> acells = mat.ls.getAdjacentCellVector();
    
    //#pragma omp parallel for schedule(dynamic)  
    for(size_t j=0; j<acells.size(); ++j) {    
      const IV idx = acells[j];
      // Check for interaction with all remaining materials
      bool isInteracting = false;
      for(MaterialIterator m2it = mats.begin(); m2it!=mats.end(); ++m2it) {
	if(*m2it==*mit) continue;
	Material& mat2 = *(*m2it);

	// Current distance between material interfaces
	const double dist = mat.ls.getPhi(idx) + mat2.ls.getPhi(idx);

	// If materials are more than a cell apart, no interaction possible within timestep...
	if(fabs(dist)>dom.dx) continue;

	// Compute an "average" normal using the two materials (remember they're pointing in opp. directions)
	DV normal = mat.ls.getNormal(idx) - mat2.ls.getNormal(idx);
	const double lnorm = dot(normal,normal);

	// If combined normal has length < 1, normals aren't pointing towards each other
	if(lnorm<0.99) continue;

	normal/=sqrt(lnorm);
	
	
	// For welded interfaces, materials interact if they're already "close enough"
	// where "close enough" is determined by weldSeparationTol
	const double weldSeparationTol = 0.1;
	if(ifaceType==Weld && dist < dom.dx * weldSeparationTol) {
	  isInteracting = true;
	}
	else {
	  // Compute rate at which interfaces are closing
	  const double closing = TestImpact(mat.sys, mat2.sys, mat.sol(idx), mat2.sol(idx), normal) + 1e-3;
	  
	  // Determine whether gap will be closed in first half of current time-step
	  if(0.5*closing*dt > dist) {
	    isInteracting = true;	    
	  }
	}
	
	if(isInteracting) {
	  const DV vel = CalcGhostStateInteraction(mat.sys, mat2.sys, mat.sol(idx), mat2.sol(idx), normal, ifaceType);
	  
	  // Compute "correction" velocity such that we move towards matching level set positions
	  // Note that this is an additional source of conservation error...
	  //! TODO - should make this dependent on the vacuum speeds computed 
	  //  previously. i.e. object with greatest unconstrained speed should have largest correction...
	  const double corr = 0.0*dist;
	  const DV corrVel = corr/dt * normal;
	  
	  // Set level-set velocities
	  mat.ls.setVelocity(idx, vel + corrVel);
	  mat2.ls.setVelocity(idx, vel - corrVel);	
	  break;
	}
      }
      if(!isInteracting) {
	// No interaction with other materials, solve material-vacuum problem
	const DV norm = mat.ls.getNormal(idx);
	
	const DV vel = CalcGhostStateVacuum(mat.sys, mat.sol(idx), norm);
	mat.ls.setVelocity(idx, vel);
      }      
    }
  }

  // Extrapolate star state to remaining ghosts
  for(MaterialIterator mit = mats.begin(); mit!=mats.end(); ++mit) {
    Material& mat = *(*mit);
    mat.ls.flagGhosts(1);  
    mat.ls.extrapolateToGhosts(mat.sol, mat.sys);
    mat.ls.extrapolateToGhosts(mat.straindot, mat.sys);
    mat.ls.extrapolateVelocity();
  }

#if 0
  for(MaterialIterator mit = mats.begin(); mit!=mats.end(); ++mit) {
    Material& mat = *(*mit);

    const vector<IV> acells = mat.ls.getGhostCellVector();
    
    //#pragma omp parallel for schedule(dynamic)  
    for(size_t j=0; j<acells.size(); ++j) {    
      const IV idx = acells[j];
      
      // No interaction with other materials, solve material-vacuum problem
      const DV norm = mat.ls.getNormal(idx);
      const DV vel = CalcGhostStateVacuum(mat.sys, mat.sol(idx), norm);
    }
  }
#endif

}

//! Extrapolate interior to adjacent and ghost cells
void extrapolateGhosts(vector<Material*>& mats)
{
  for(MaterialIterator mit = mats.begin(); mit!=mats.end(); ++mit) {
    Material& mat = *(*mit);
    //  mat.ls.extrapolateEntropyToBoundary(mat.sol, mat.sys);
    mat.ls.extrapolateToAdjacent(mat.sol, mat.sys);
    mat.ls.extrapolateToAdjacent(mat.straindot, mat.sys);
    mat.ls.flagGhosts(1);  
    mat.ls.extrapolateToGhosts(mat.sol, mat.sys);
    mat.ls.extrapolateToGhosts(mat.straindot, mat.sys);
  }
  
}

string qStr(const string& str, const int m) {
  ostringstream oss;
  oss << "\"" << str << "_" << m << "\" ";
  return oss.str();
}

string qStr(const string& str) {
  ostringstream oss;
  oss << "\"" << str << "\" ";
  return oss.str();
}

//! Output current solution in a suitable format for tecplot / visit
void outputTec(const string& outdir,
	       const Domain& dom, 
	       const vector<Material*>& mats,
	       const int step) 
{
  // x y rho u v w sig00 sig11 sig22 sig01 sig02 sig12 S
  char fname[20];
  sprintf(fname, "out%03d.tec", step);
  string filename = outdir + string(fname);
  ofstream outFile(filename.c_str());

  outFile << "TITLE=" << qStr("Elastic2D Output") << "\n";
  outFile << "VARIABLES=" << qStr("X") << qStr("Y");

  for(unsigned int m=1; m<=mats.size(); ++m) {
    outFile << qStr("Density",m)
	    << qStr("U",m)
	    << qStr("V",m)
	    << qStr("T",m)
	    << qStr("S11",m)
	    << qStr("S22",m)
	    << qStr("S33",m)
	    << qStr("S12",m)
	    << qStr("P",m)
	    << qStr("Entropy",m)
	    << qStr("e_plast",m)
	    << qStr("Div0",m)
	    << qStr("Div1",m)
	    << qStr("Div2",m)
	    << qStr("phi",m);
  }
  outFile << "\n";

  outFile << "ZONE I="<<dom.Ni<<" J="<<dom.Nj<<" F=POINT\n";

  for(int j=dom.startj; j<dom.endj; ++j) {
    for(int i=dom.starti; i<dom.endi; ++i) {
      const IV C(i,j);
      const LevelSet<2>::Coord pos = mats[0]->ls.getCellPos(C);
      outFile << pos[0] << "\t" << pos[1];

      for(MaterialConstIterator mit=mats.begin(); mit!=mats.end(); ++mit) {
	const Material& mat = *(*mit);

	if(mat.ls.isInteriorGhost(C)) {

	  const ElasticState solij = mat.sol(i,j);
	  const ElasticPrimState psol = mat.sys.ToPrimitive(solij);
	  const SymmetricMatrix sigma = mat.sys.Stress(psol);			
	  const double P = -(sigma(0,0)+sigma(1,1)+sigma(2,2))/3.0;
	  const TinyVector<double, 2> div = divRhoFe(mat,dom,C);
	  outFile << "\t" << mat.sys.Density(psol);
	  outFile << "\t" << psol.velocity[0];
	  outFile << "\t" << psol.velocity[1];
	  outFile << "\t" << mat.sys.Temperature(psol);
	  outFile << "\t" << sigma(0,0);
	  outFile << "\t" << sigma(1,1);
	  outFile << "\t" << sigma(2,2);
	  outFile << "\t" << sigma(0,1);
	  outFile << "\t" << P;
	  outFile << "\t" << psol.S;
	  outFile << "\t" << psol.Epsp;
	  outFile << "\t" << div[0];
	  outFile << "\t" << div[1];
	  outFile << "\t" << div[2];
	}
	else
	  {	  
	    for(int v=0; v<14; ++v) outFile << "\t" << "nan";
	  }

	outFile << "\t" << mat.ls.getPhi(IV(i,j));

      }
      outFile << "\n";
    }
  }
  outFile.close();
}

#include "silo.h"

//! Output current solution in a suitable format for tecplot / visit
void outputSilo(const string& outdir,
		const Domain& dom, 
		const vector<Material*>& mats,
		int step,
		double t) 
{
  char fname[20];
  sprintf(fname, "out%03d.silo", step);
  string filename = outdir + string(fname);

  DBfile* file = DBCreate(filename.c_str(), DB_CLOBBER, DB_LOCAL, "Elastic2D output file", DB_HDF5);
  DBSetCompression("METHOD=GZIP");

  char* coordnames[2] = {"X", "Y"};
  float *nodes[2];
  int dimensions[2] = {dom.Ni+1, dom.Nj+1};
  nodes[0] = new float[dom.Ni+1];
  nodes[1] = new float[dom.Nj+1];
  for(int i=0; i<=dom.Ni; ++i) nodes[0][i] = dom.dx*i;
  for(int j=0; j<=dom.Nj; ++j) nodes[1][j] = dom.dy*j;
  DBoptlist *meshOpts = DBMakeOptlist(10);
  int coord = DB_CARTESIAN;
  char* xunits = "m";
  DBAddOption(meshOpts, DBOPT_COORDSYS, &coord);
  DBAddOption(meshOpts, DBOPT_CYCLE, &step);
  DBAddOption(meshOpts, DBOPT_DTIME, &t);
  DBAddOption(meshOpts, DBOPT_XUNITS, xunits);
  DBAddOption(meshOpts, DBOPT_YUNITS, xunits);
  DBAddOption(meshOpts, DBOPT_ZUNITS, xunits);
  DBPutQuadmesh(file, "mesh", coordnames, nodes, dimensions, 2, DB_FLOAT, DB_COLLINEAR, meshOpts);
  DBFreeOptlist(meshOpts);

  delete[] nodes[0];
  delete[] nodes[1];

  int cells[2] = {dom.Ni, dom.Nj};
  int m=1;

  DBoptlist *baseOpts = DBMakeOptlist(10);
  DBAddOption(baseOpts, DBOPT_COORDSYS, &coord);
  DBAddOption(baseOpts, DBOPT_CYCLE, &step);
  DBAddOption(baseOpts, DBOPT_DTIME, &t);
    
  char **names = new char*[5*mats.size()];
  char **defs = new char*[5*mats.size()];
  int* types = new int[5*mats.size()];

  for(size_t i=0; i<5*mats.size(); ++i) 
  {
    names[i] = new char[100];
    defs[i] = new char[1000];
  }

  for(MaterialConstIterator mit=mats.begin(); mit!=mats.end(); ++mit, ++m) {
    const Material& mat = *(*mit);
    char name[20];

    float *rho = new float[dom.Ni*dom.Nj];
    float *epsdot = new float[dom.Ni*dom.Nj];
    float *S = new float[dom.Ni*dom.Nj];
    float *epsP = new float[dom.Ni*dom.Nj];
    float *T = new float[dom.Ni*dom.Nj];
    float* sigmaOut[4];
    char* sigmaNames[4];
    for(int i=0; i<4; ++i) {
      sigmaNames[i] = new char[20];
      sigmaOut[i] = new float[dom.Ni*dom.Nj];
    }

    char* velNames[2];
    float* vel[2];
    for(int i=0; i<2; ++i) {
      velNames[i] = new char[20];
      vel[i] = new float[dom.Ni*dom.Nj];
    }

    for(int i=0; i<dom.Ni; ++i) {
      for(int j=0; j<dom.Nj; ++j) {
	const IV C(i+dom.starti, j+dom.startj);
	const int cell = i+j*dom.Ni;
	if(mat.ls.isInteriorGhost(C)) 
	{	 
	  const ElasticState solij = mat.sol(i+dom.starti,j+dom.startj);
	  const ElasticPrimState psol = mat.sys.ToPrimitive(solij);
	  const SymmetricMatrix sigma = mat.sys.Stress(psol);			
//	  const TinyVector<double, 2> div = divRhoFe(mat,dom,C);
	  
	  rho[cell] = mat.sys.Density(solij);
	  S[cell] = mat.sys.InternalEnergy(psol);
	  epsP[cell] = psol.Epsp;
	  epsdot[cell] = mat.straindot(i+dom.starti, j+dom.startj);
	  T[cell] = mat.sys.Temperature(psol);
	  vel[0][cell] = psol.velocity[0];
	  vel[1][cell] = psol.velocity[1];
	  sigmaOut[0][cell] = sigma(0,0);
	  sigmaOut[1][cell] = sigma(0,1);
	  sigmaOut[2][cell] = sigma(1,1);
	  sigmaOut[3][cell] = sigma(2,2);
	}
	else 
	{
	  S[cell] = 0;
	  epsP[cell] = 0;
	  epsdot[cell] = 0;
	  T[cell] = 0;
	  rho[cell] = 0;
	  sigmaOut[0][cell] = 0;
	  sigmaOut[1][cell] = 0;
	  sigmaOut[2][cell] = 0;
	  sigmaOut[3][cell] = 0;
	  vel[0][cell] = 0;
	  vel[1][cell] = 0;
	}
	
      }
    }

    sprintf(name, "rho_%d", m);
    DBPutQuadvar1(file, name, "mesh", rho, cells, 2, NULL, 0, DB_FLOAT, DB_ZONECENT, baseOpts); 
    delete[] rho;

    sprintf(name, "T_%d", m);
    DBPutQuadvar1(file, name, "mesh", T, cells, 2, NULL, 0, DB_FLOAT, DB_ZONECENT, baseOpts); 
    delete[] T;

    sprintf(name, "ie_%d", m);
    DBPutQuadvar1(file, name, "mesh", S, cells, 2, NULL, 0, DB_FLOAT, DB_ZONECENT, baseOpts); 
    delete[] S;

    sprintf(name, "epsp_%d", m);
    DBPutQuadvar1(file, name, "mesh", epsP, cells, 2, NULL, 0, DB_FLOAT, DB_ZONECENT, baseOpts); 
    delete[] epsP;

    sprintf(name, "epsdot_%d", m);
    DBPutQuadvar1(file, name, "mesh", epsdot, cells, 2, NULL, 0, DB_FLOAT, DB_ZONECENT, baseOpts); 
    delete[] epsdot;

    sprintf(sigmaNames[0], "stress_xx_%d", m);
    sprintf(sigmaNames[1], "stress_xy_%d", m);
    sprintf(sigmaNames[2], "stress_yy_%d", m);
    sprintf(sigmaNames[3], "stress_zz_%d", m);
    for(int i=0; i<4; ++i)
      DBPutQuadvar1(file, sigmaNames[i], "mesh", sigmaOut[i], cells, 2, NULL, 0, DB_FLOAT, DB_ZONECENT, baseOpts); 

    // Velocity
    sprintf(velNames[0], "v_x_%d", m);
    sprintf(velNames[1], "v_y_%d", m);
    for(int i=0; i<2; ++i)
      DBPutQuadvar1(file, velNames[i], "mesh", vel[i], cells, 2, NULL, 0, DB_FLOAT, DB_ZONECENT, baseOpts); 

    // Add derived variable expressions
    
    sprintf(names[5*(m-1)+0], "stress_%d", m);
    sprintf(defs[5*(m-1)+0], "{{%s, %s, 0}, {%s, %s, 0}, {0, 0, %s}}",
	    sigmaNames[0], sigmaNames[1], sigmaNames[1], 
	    sigmaNames[2], sigmaNames[3]);
    types[5*(m-1)+0] = DB_VARTYPE_TENSOR;

    sprintf(names[5*(m-1)+1], "P_%d", m);
    sprintf(defs[5*(m-1)+1], "-(%s+%s+%s)/3.0", sigmaNames[0], sigmaNames[2], sigmaNames[3]);
    types[5*(m-1)+1] = DB_VARTYPE_SCALAR;

    sprintf(names[5*(m-1)+3], "stress2D_%d", m);
    sprintf(defs[5*(m-1)+3], "{{%s, %s}, {%s, %s}}", sigmaNames[0], sigmaNames[1], sigmaNames[1], sigmaNames[2]);
    types[5*(m-1)+3] = DB_VARTYPE_TENSOR;

    sprintf(names[5*(m-1)+4], "vel_%d", m);
    sprintf(defs[5*(m-1)+4], "{%s, %s}", velNames[0], velNames[1]);
    types[5*(m-1)+4] = DB_VARTYPE_VECTOR;

    for(int i=0; i<2; ++i) {
      delete[] vel[i];
      delete[] velNames[i];
    }
    
    for(int i=0; i<4; ++i) {
      delete[] sigmaOut[i];
      delete[] sigmaNames[i];
    }

   
    // Level set
    float *phi = new float[dom.Ni*dom.Nj];
    float *nx = new float[dom.Ni*dom.Nj];
    float *ny = new float[dom.Ni*dom.Nj];
    for(int i=0; i<dom.Ni; ++i) {
      for(int j=0; j<dom.Nj; ++j) {
	phi[i+j*dom.Ni] = mat.ls.getPhi(IV(i+dom.starti, j+dom.startj));
	nx[i+j*dom.Ni] = mat.ls.getNormal(IV(i+dom.starti, j+dom.startj))[0];
	ny[i+j*dom.Ni] = mat.ls.getNormal(IV(i+dom.starti, j+dom.startj))[1];
      }
    }

    sprintf(name, "phi_%d", m);
    DBPutQuadvar1(file, name, "mesh", phi, cells, 2, NULL, 0, DB_FLOAT, DB_ZONECENT, baseOpts); 
    sprintf(name, "nx_%d", m);
    DBPutQuadvar1(file, name, "mesh", nx, cells, 2, NULL, 0, DB_FLOAT, DB_ZONECENT, baseOpts); 
    sprintf(name, "ny_%d", m);
    DBPutQuadvar1(file, name, "mesh", ny, cells, 2, NULL, 0, DB_FLOAT, DB_ZONECENT, baseOpts); 
    delete[] phi;
    delete[] nx;
    delete[] ny;
  }
  DBFreeOptlist(baseOpts);
  
  DBPutDefvars(file, "defvars", 5*mats.size(), names, types, defs, NULL);

  for(size_t i=0; i<5*mats.size(); ++i) {
    delete[] defs[i];
    delete[] names[i];
  }   
  delete[] defs;
  delete[] names;
  delete[] types;

  DBClose(file);
}

void syncVelocity(Material& mat, const Domain& dom) {
  for(int j=dom.startj; j<dom.endj; ++j) {
    for(int i=dom.starti; i<dom.endi; ++i) {
      const IV pos(i,j);
      if(mat.ls.m_status(i,j)=='b') {
	const double rho = mat.sys.Density(mat.sol(i,j));
	mat.ls.setVelocity(IV(i,j),DV(mat.sol(i,j).momentum[0]/rho,mat.sol(i,j).momentum[1]/rho));
      }
    }
  }
}

//! Apply the plastic correction update to all cells
void plasticSolve(Material& mat, const Domain& dom) {
#pragma omp parallel for schedule(dynamic)
  for(int j=dom.startj; j<dom.endj; ++j) {
    for(int i=dom.starti; i<dom.endi; ++i) {
      const IV pos(i,j);
      if(mat.ls.isInterior(pos)) {
	mat.sys.PlasticUpdate(mat.sol(i,j),mat.straindot(i,j));
      }
    }
  }

}

//! Compute stable timestep based on CFL condition
double getMinDt(Material& mat, const Domain& dom) {

  double sharedmindt = std::numeric_limits<double>::max();

#pragma omp parallel
  {
    double mindt = std::numeric_limits<double>::max();
    

#pragma omp for nowait schedule(dynamic)
    for(int j=dom.startj; j<dom.endj; ++j) {
      for(int i=dom.starti; i<dom.endi; ++i) {			
	if(mat.ls.isInterior(IV(i,j))) {
	  try 
	    {	    
	      const DV3 sp=mat.sys.getMaxWaveSpeeds(mat.sol(i,j));
	      mindt = std::min(dom.dx/sp[0], mindt);
	      mindt = std::min(dom.dy/sp[1], mindt);
	    }
	  catch (char const* exc) 
	    {
	      std::cerr << exc;
	      const LevelSet<2>::Coord pos = mat.ls.getCellPos(IV(i,j));
	      std::cerr << "In material " << mat.name << " at i = " << i << ", j = " << j << ", x = " << pos[0] << ", y = " << pos[1] << "\n\n";
	      std::cerr << "Conservative state = " << mat.sol(i,j) << "\n";
	      std::cerr << "Primitive state = " << mat.sys.ToPrimitive(mat.sol(i,j)) << "\n";	   
	      vector<Material*> mats;
	      mats.push_back(&mat);
	      cerr << "Dumping state for material to current directory as stage 9999\n";
	      extrapolateGhosts(mats);
	      outputSilo("./", dom, mats, 9999,-1.0);
	      exit(1);
	    }
	  
	}
      }      
    }
    
#pragma omp critical
    {
      sharedmindt = std::min(sharedmindt, mindt);
    }
  }

  return sharedmindt;
}


//! Output energy integrals 
void outputIgrals(const string& outdir, const vector<Material*>& mats, const Domain& dom, const double t, const bool axisym) 
{
  double totenergy0 = 0.0;
  for(vector<Material*>::const_iterator mit=mats.begin(); mit!=mats.end(); ++mit) 
    {
      totenergy0 += (**mit).energy0;
    }


  string fname = outdir + "igral.out";
  ofstream outFile;
  if(t==0.0)  
    {    
      outFile.open(fname.c_str(), ios_base::out);
      outFile << "# t\t";
      for(vector<Material*>::const_iterator mit=mats.begin(); mit!=mats.end(); ++mit) 
	{
	  const Material& mat = **mit;
	  outFile << setprecision(10);
	  outFile << mat.name<<"_mass/" << mat.mass0 << "\t";
	  outFile << mat.name<<"_E/" << totenergy0 << "\t";
	  outFile << mat.name<<"_KE/" << totenergy0 << "\t";
	  outFile << mat.name<<"_IE/" << totenergy0 << "\t";
	}
      outFile << "\n";
    } else 
    {    
      outFile.open(fname.c_str(), ios_base::app);
    }
  
  outFile << t << "\t";

  for(vector<Material*>::const_iterator mit=mats.begin(); mit!=mats.end(); ++mit) 
    {
      const Material& mat = **mit;
      double E=0;
      double KE=0;
      double IE=0;
      double M=0;
      for(int i=dom.starti; i<dom.endi; ++i) {
	for(int j=dom.startj; j<dom.endj; ++j) {
	  const IV C(i,j);
	  if(mat.ls.isInteriorGhost(C)) {
	    const ElasticState solij = mat.sol(i,j);
	    const ElasticPrimState psol = mat.sys.ToPrimitive(solij);
	    const double rho = mat.sys.Density(solij);
	    double vf = mat.ls.getCellFraction(C);
	    if(axisym) {
	      const DV pos=mat.ls.getCellPos(C);
	      const double r1=pos[0]+0.5*dom.dx;
	      const double r0=pos[0]-0.5*dom.dx;		
	      if(vf>0.0 && vf<1.0) 
	      {
		double G[4];
		int m=0;
		for(int lj=C[1]-1; lj<=C[1]; ++lj) 
		{
		  for(int li=C[0]-1; li<=C[0]; ++li,++m) 
		  {
		    G[m]=0.0;
		    for(int k=li; k<li+2; ++k) 
		    {
		      for(int l=lj; l<lj+2; ++l) 
		      {      
			G[m]+=0.25*mat.ls.getPhi(IV(k,l));
		      }
		    }
		  }
		}		
		vf = mat.ls.volumeAxi2D(G, r0, r1);
	      }	      
	      else 
	      {
		vf*=M_PI*(r1*r1-r0*r0)*dom.dy;			
	      }
	      
	    } else {
	      vf*=dom.dx*dom.dy;
	    }
	    
	    M += rho*vf;
	    E += solij.rhoE*vf;
	    const double K = 0.5*rho*dot(psol.velocity,psol.velocity);
	    KE += K*vf;
	    IE += (solij.rhoE - K)*vf;
	  }
	
	}    
      }
    
      outFile << setprecision(10);
      outFile << setw(14) << M/mat.mass0 << "\t" << setw(14) << E/totenergy0 << "\t" << setw(14) << KE/totenergy0 << "\t" << setw(14) << IE/totenergy0 << "\t";
    }

  outFile << "\n";
}

//! Compute initial integrals 
void initialIgrals(vector<Material*>& mats, const Domain& dom, const bool axisym) 
{
  cerr << "Computing initial integrals:\n";
  for(vector<Material*>::iterator mit=mats.begin(); mit!=mats.end(); ++mit) 
    {
      Material& mat = **mit;
      double E=0;
      double M=0;
      for(int i=dom.starti; i<dom.endi; ++i) {
	for(int j=dom.startj; j<dom.endj; ++j) {
	  const IV C(i,j);
	  if(mat.ls.isInteriorGhost(C)) {
	    const ElasticState solij = mat.sol(i,j);
	    const ElasticPrimState psol = mat.sys.ToPrimitive(solij);
	    const double rho = mat.sys.Density(solij);
	    double vf = mat.ls.getCellFraction(C);
	    if(axisym) {
	      const DV pos=mat.ls.getCellPos(C);
	      const double r1=pos[0]+0.5*dom.dx;
	      const double r0=pos[0]-0.5*dom.dx;		
	      if(vf>0.0 && vf<1.0) 
	      {
		double G[4];
		int m=0;
		for(int lj=C[1]-1; lj<=C[1]; ++lj) 
		{
		  for(int li=C[0]-1; li<=C[0]; ++li,++m) 
		  {
		    G[m]=0.0;
		    for(int k=li; k<li+2; ++k) 
		    {
		      for(int l=lj; l<lj+2; ++l) 
		      {      
			G[m]+=0.25*mat.ls.getPhi(IV(k,l));
		      }
		    }
		  }
		}		
		vf = mat.ls.volumeAxi2D(G, r0, r1);
	      }	      
	      else 
	      {
		vf*=M_PI*(r1*r1-r0*r0)*dom.dy;			
	      }
	      
	    } else {
	      vf*=dom.dx*dom.dy;
	    }
	    
	    M += rho*vf;
	    E += solij.rhoE*vf;
	  }
	
	}    
      }
    
      mat.mass0 = M;
      mat.energy0 = E;
    }
}

//! Output state at lagrangian stations
void outputStations(const string& outdir, const vector<Material*>& mats, const Domain& dom, const double t, const double dt) 
{
  string fname = outdir+"stations.out";
  ofstream statOut;
  if(t==0.0) 
    {    
      statOut.open(fname.c_str(), ios_base::out);
      statOut << "# t\t";
      for(vector<Material*>::const_iterator mit=mats.begin(); mit!=mats.end(); ++mit) 
	{
	  const Material& mat = **mit;
	  statOut << mat.name<<"_xvel" << "\t";
	  statOut << mat.name<<"_yvel" << "\t";
	  statOut << mat.name<<"_pres" << "\t";
	  statOut << mat.name<<"_dens" << "\t";
	  statOut << mat.name<<"_sigx" << "\t";
	}
      statOut << "\n";
    } else 
    {
      statOut.open(fname.c_str(), ios_base::app);
    }
  
  // Output station states and update station positions
  statOut << t << "\t";
  for(size_t i=0; i<stations.size(); ++i) 
    {           
      DV pos=stations[i];
      IV idx;
      statOut << pos[0] << "\t" << pos[1] << "\t";
    
      idx[0]=int(floor(pos[0]/dom.dx))+dom.GC;
      idx[1]=int(floor(pos[1]/dom.dy))+dom.GC;
      bool found=false;
      for(MaterialConstIterator mit=mats.begin(); mit!=mats.end(); ++mit) {
	const Material& mat = **mit;
	if(mat.ls.isInterior(idx))
	  {
	    found = true;
	
	    // Linearly interpolate to position
	    const double fx = pos[0]/dom.dx-idx[0]-0.5+dom.GC;
	    const double fy = pos[1]/dom.dy-idx[1]-0.5+dom.GC;
	    const double afx = fabs(fx);
	    const double afy = fabs(fy);
	    assert(afx<=0.5);
	    assert(afy<=0.5);
	    IV idxx = idx, idxxy;
	    IV idxy = idx;
	    idxx[0] += (fx<0.?-1:1);
	    idxy[1] += (fy<0.?-1:1);
	    idxxy[0] = idxx[0];
	    idxxy[1] = idxy[1];
	    const ElasticState interp = afx * (afy*mat.sol(idxxy)+(1-afy)*mat.sol(idxx)) +
	      (1-afx) * (afy*mat.sol(idxy)+(1-afy)*mat.sol(idx));
	    const ElasticPrimState psol = mat.sys.ToPrimitive(interp);
	    const SymmetricMatrix sigma = mat.sys.Stress(psol);			
	    statOut << psol.velocity[0] << "\t" << psol.velocity[1] << "\t" << (sigma(0,0)+sigma(1,1)+sigma(2,2))/(-3.0) << "\t" << mat.sys.Density(psol) << "\t" << sigma(0,0) << "\t";
	    pos[0]+=psol.velocity[0]*dt;
	    pos[1]+=psol.velocity[1]*dt;
	    stations[i] = pos;
	    break;
	  }	
      }
      if(!found) 
	{
	  cerr << "Warning: station " << i << " does not lie within any material\n";
	  cerr << "pos = " << pos << endl;
	  cerr << "idx = " << idx << endl;
	  statOut << "0.0\t0.0\t0.0\t0.0\t0.0\t";
	}
    }
  statOut << "\n";
}

//! Compute the inner and outer radius at 0, 45 and 90 degrees
void outputRadius(const string& outdir, const Material& mat, const Domain& dom, const double t) 
{
  string fname = outdir+"radius."+mat.name+".out";
  ofstream outFile;
  if(t==0.0) 
    outFile.open(fname.c_str(), ios_base::out);
  else 
    outFile.open(fname.c_str(), ios_base::app);


  const int Ci = dom.starti;
  const int Cj = dom.startj;
  const DV cpos(0,0);
  outFile << t << "\t";
  { 
    // Radius along x-axis
    int i=Ci;
    for(int j=Cj; j<dom.endj; ++j) {
      const IV CL(i,j);
      const IV CR(i,j+1);
      const double phiL = mat.ls.getPhi(CL);
      const double phiR = mat.ls.getPhi(CR);
      
      if(phiL*phiR<0.0) 
	{	
	  const double frac = phiL / (phiL - phiR);
	  const LevelSet<2>::Coord pos = (1 - frac)*mat.ls.getCellPos(CL) + (frac)*mat.ls.getCellPos(CR) - cpos;
	  const double radius = sqrt(dot(pos,pos));
	  //	const double angle = atan2(pos[1], pos[0]);
	  outFile << radius << "\t";	
	}
    }
  }
  
  {
    // Radius along y-axis
    int j=Cj;
    for(int i=Ci; i<dom.endi; ++i) {
      const IV CL(i,j);
      const IV CR(i+1,j);
      const double phiL = mat.ls.getPhi(CL);
      const double phiR = mat.ls.getPhi(CR);

      if(phiL*phiR<0.0) 
	{	
	  const double frac = phiL / (phiL - phiR);
	  const LevelSet<2>::Coord pos = (1 - frac)*mat.ls.getCellPos(CL) + (frac)*mat.ls.getCellPos(CR) - cpos;
	  const double radius = sqrt(dot(pos,pos));
	  //	const double angle = atan2(pos[1], pos[0]);
	  outFile << radius << "\t";	
	}
    }
  }
  
  {
    // Radius along x=y
    for(int i=Ci; i<dom.endi; ++i) {
      int j=i;
      const IV CL(i,j);
      const IV CR(i+1,j+1);
      const double phiL = mat.ls.getPhi(CL);
      const double phiR = mat.ls.getPhi(CR);

      if(phiL*phiR<0.0) 
	{	
	  const double frac = phiL / (phiL - phiR);
	  const LevelSet<2>::Coord pos = (1 - frac)*mat.ls.getCellPos(CL) + (frac)*mat.ls.getCellPos(CR) - cpos;
	  const double radius = sqrt(dot(pos,pos));
	  //	const double angle = atan2(pos[1], pos[0]);
	  outFile << radius << "\t";	
	}
    }
  }
  
  outFile << "\n";
  outFile.flush();
}

string getDirName(const string& prefix, const Domain& dom, const InterfaceType& ifaceType, const bool axiSym) 
{
  stringstream filenameSS;
  filenameSS << "out/" << prefix << "_" << dom.Ni << "x" << dom.Nj << "_";
  if(ifaceType==Stick) 
    {
      filenameSS << "stick";
    }
  else if(ifaceType==Slip) 
    {
      filenameSS << "slip";
    }
  else if(ifaceType==Weld) 
    {
      filenameSS << "weld";
    }
  if(axiSym) {
    filenameSS << "_axi";
  }

  if(interpMethod == System::ComponentConsInterp)
  {
    filenameSS << "_cons";
  }
  else if(interpMethod == System::ComponentPrimInterp)
  {
    filenameSS << "_prim";
  }
  else if(interpMethod == System::InvariantInterp)
  {
    filenameSS << "_inv";
  }
  else if(interpMethod == System::ColumnLengthInterp)
  {
    filenameSS << "_col";
  }
  else if(interpMethod == System::LengthDensityInterp)
  {
    filenameSS << "_coldens";
  }

  filenameSS << "/";
  return filenameSS.str();
}    

int main(int argc, char **argv){

  // Enable floating point error checking
    feenableexcept(FE_INVALID);
    feenableexcept(FE_DIVBYZERO);

  vector<Material*> mats;
  Domain dom;
  double tf, outFreq;
  string outdir;
  int Ni;
  InterfaceType ifaceType;
  ICType icType;
  double CFL = 0.7;      // CFL condition 
  int reinitFreq = 15;    // Level set reinitialization frequency
  bool axiSym = false;   // Axisymmetric simulation?
  int outInterval = 100000000;
  double fixed_dt = -1.0; // Fixed time-step?

  if(argc!=6) 
    {
      cerr << "Syntax: Elastic2D problem-type x-resolution interface-type axisymmetry interp-method\n";
      exit(1);
    } else 
    {
      const string icStr = string(argv[1]);
      if(icStr == "friction") 
	{
	  icType = Friction;
	}
      else if(icStr == "bending") 
	{
	  icType = BendingBeam;
	}
      else if(icStr == "hblow") 
	{
	  icType = HBLowSpeed;
	}
      else if(icStr == "hbhigh") 
	{
	  icType = HBHighSpeed;
	} 
      else if(icStr == "bar") 
	{
	  icType = BarPlate;
	} 
      else if(icStr == "extrap") 
	{
	  icType = ExtrapTest;
	} 
      else if(icStr == "shell") 
	{
	  icType = BerylliumShell;
	} 
      else if(icStr == "taylor") 
	{
	  icType = TaylorAnvil;
	} 
      else if(icStr == "symmetric") 
	{
	  icType = SymmetricTaylor;
	} 
      else if(icStr == "tantalum") 
	{
	  icType = TaylorAnvilTantalum;
	} 
      else if(icStr == "radial") 
	{
	  icType = Radial;
	} 
      else 
	{
	  cerr << "Unknown problem type \'" << icStr << "\' - options are friction, hblow, hbhigh, taylor, shell, radial\n";
	  exit(1);
	}

      Ni = atoi(argv[2]);

      const string ifaceStr = string(argv[3]);
      if(ifaceStr == "stick") 
	{
	  ifaceType = Stick;
	}
      else if(ifaceStr == "slip") 
	{
	  ifaceType = Slip;
	}
      else if(ifaceStr == "weld") 
	{
	  ifaceType = Weld;
	} else 
	{
	  cerr << "Unknown interface type \'" << ifaceStr << "\' - options are stick, slip, or weld\n";
	  exit(1);
	}

      const string axiStr = string(argv[4]);
      if(axiStr == "true") 
	{
	  axiSym = true;
	}
      else if(axiStr == "false") 
	{
	  axiSym = false;
	}
      else 
	{
	  cerr << "Unknown axisymmetry mode \'" << axiStr << "\' - options are true or false\n";
	  exit(1);
	}
      
      const string itStr = string(argv[5]);
      if(itStr == "cons") 
      {
	interpMethod = System::ComponentConsInterp;
      } 
      else if(itStr == "prim") 
      {      
	interpMethod = System::ComponentPrimInterp;
      }
      else if(itStr == "inv") 
      {      
	interpMethod = System::InvariantInterp;
      }
      else if(itStr == "col") 
      {      
	interpMethod = System::ColumnLengthInterp;
      }
      else if(itStr == "colrho") 
      {      
	interpMethod = System::LengthDensityInterp;
      }
      else 
      {
	cerr << "Unknown interpolation mode \'" << itStr << "\' - options are cons, prim, inv, col, colrho\n";
	exit(1);
      }
    }
  
  
  // Set the initial condition for both materials  
  if(icType == HBLowSpeed) {

    // Howell-Ball low speed bar/plate impact:
    tf = 5e-6;                       // Final time
    outFreq = 1e-7;                  // Output frequency    
    dom = Domain(Ni, (Ni*3)/2, 3, 0.02); // Domain
    fixed_dt = dom.dx / 1.01e4;

    outdir = getDirName("HBLowSpeed", dom, ifaceType, axiSym);

    ElasticEoS* alu_eos = (ElasticEoS*) new Romenski("aluminium");
    PlasticModel* alu_plastic = (PlasticModel*) new PerfectPlastic("aluminium");
    
    // Target:
    mats.push_back(new Material("target", dom, alu_eos, alu_plastic));
    ICBlock(*mats[0], dom, Rect(DV(-0.019, 0.006),DV(0.019, 0.028)), DV3(0,0,0));
    
    // Projectile:
    mats.push_back(new Material("projectile", dom, alu_eos, alu_plastic));
    ICBlock(*mats[1], dom, Rect(DV(-0.006, 0.001),DV(0.006, 0.006)), DV3(0,400,0));

    for(int i=0; i<5; ++i) {      
      stations.push_back(DV(0.0, 7.8125e-3+i*3.625e-3));
    }    

  } else if(icType == Friction) {
    
    // Barton friction-target test configuration (copper plate impacts split aluminium steel target)
    tf = 15e-6;                        // Final time
    outFreq = 1e-7;                   // Output frequency   
    //    outInterval = 5;
    dom = Domain(Ni, (Ni*3)/2, 3, 0.03); // Domain

    outdir = getDirName("Friction", dom, ifaceType, axiSym);

    stations.push_back(DV(0.000, 0.0071));
    stations.push_back(DV(0.015, 0.0071));

    // Aluminium target:
    mats.push_back(new Material("aluminium", dom, new Romenski("aluminium"), new PerfectPlastic("aluminium")));
    ICQuad(*mats[0], dom, DV(-0.01,0.007),DV(-0.01,0.037),DV(0.012605,0.037),DV(0.007395,0.007), DV3(0,0,0.0));

    // Steel target:
    mats.push_back(new Material("steel", dom, new Romenski("steel"), new PerfectPlastic("steel")));
    ICQuad(*mats[1], dom, DV(0.007395,0.007),DV(0.012605,0.037),DV(0.025,0.037),DV(0.025,0.007), DV3(0,0,0.0));
    
    // Flyer plate
    mats.push_back(new Material("copper", dom, new Romenski("copper"), new PerfectPlastic("copper")));
    ICBlock(*mats[2], dom, Rect(DV(-0.027,0.037),DV(0.027,0.041)), DV3(0,-202,0.0));

  } else if(icType == BendingBeam) {
    
    // Barton friction-target test configuration (copper plate impacts split aluminium steel target)
    tf = 120;                        // Final time
    outFreq = 10;                   // Output frequency   
    //    outInterval = 5;
    dom = Domain(Ni, Ni*3/7, 3, 7); // Domain

    outdir = getDirName("BendingBeam", dom, ifaceType, axiSym);

    // Target:
    mats.push_back(new Material("beam", dom, new Osborne("beryllium"), 
				new NullPlastic));
    ICBeam(*mats[0], dom, Rect(DV(0.5, 1),DV(6.5, 2)));

  } else if(icType == HBHighSpeed) {

    // Howell-Ball high speed ball/plate impact:

    tf = 5e-6;                        // Final time
    outFreq = 2e-7;                   // Output frequency   
    dom = Domain(Ni, Ni, 3, 0.026); // Domain
    reinitFreq = 2;

    outdir = getDirName("HBHighSpeed", dom, ifaceType, axiSym);
        
    // Plate:
    mats.push_back(new Material("plate", dom, new MGShock("aluminium"), new PerfectPlastic("aluminium")));
    ICBlock(*mats[0], dom, Rect(DV(-0.025,0.011),DV(0.025,0.013)), DV3(0,0,0.0));
    
    // Projectile:
    mats.push_back(new Material("ball", dom, new MGShock("aluminium"), new PerfectPlastic("aluminium")));

    ICBall(*mats[1], dom, DV(0.0,0.0055), 0.005, DV3(0,3100,0));

  } else if(icType == BarPlate) {

    // Bar-plate low speed impact:

    tf = 8e-5;                        // Final time
    outFreq = 5e-6;                   // Output frequency   
    dom = Domain(Ni, Ni, 3, 0.02601); // Domain

    outdir = getDirName("BarPlate", dom, ifaceType, axiSym);
        
    // Plate:
    mats.push_back(new Material("plate", dom, new Romenski("steel"), new PerfectPlastic("steel")));
    ICBlock(*mats[0], dom, Rect(DV(0.011,-0.025),DV(0.013,0.025)), DV3(0,0,0.0));
    
    // Projectile:
    mats.push_back(new Material("projectile", dom, new Romenski("copper"), new PerfectPlastic("copper")));
    ICBlock(*mats[1], dom, Rect(DV(0.001,-0.0025),DV(0.011,0.0025)), DV3(800,0,0));

  } else if(icType == SymmetricTaylor) {

    // Pagosa Taylor Anvil:
    tf = 100e-6;                      // Final time
    outFreq = 4e-6;                   // Output frequency   
    dom = Domain(Ni, Ni*10, 3, 0.011); // Domain
    reinitFreq = 40;

    outdir = getDirName("SymmetricTaylor", dom, ifaceType, axiSym);
        
    // Rod:
    mats.push_back(new Material("rod", dom, new VinetRose("copper"), new PathDependent("copper")));
    ICBlock(*mats[0], dom, Rect(DV(-0.005,-0.1),DV(0.005,0.1)), DV3(0,-395/2.0,0));

  } else if(icType == TaylorAnvil) {

    // Pagosa Taylor Anvil:
    tf = 100e-6;                        // Final time
    outFreq = 5e-6;                   // Output frequency   
    dom = Domain(Ni, Ni*3, 3, 0.009); // Domain
    reinitFreq = 40;

    outdir = getDirName("TaylorAnvil", dom, ifaceType, axiSym);
        
    // Rod:
    mats.push_back(new Material("rod", dom, new VinetRose("copper"), new ZerilliArmstrong("copper")));
    ICBlock(*mats[0], dom, Rect(DV(-0.00381,-0.0254),DV(0.00381,0.0254)), DV3(0,-190,0));

  } else if(icType == TaylorAnvilTantalum) {

    // Andy's tantalum taylor anvil
    tf = 100e-6;                        // Final time
    outFreq = 5e-6;                   // Output frequency   
    dom = Domain(Ni, Ni*3, 3, 0.009); // Domain
    reinitFreq = 40;

    outdir = getDirName("TantalumAnvil", dom, ifaceType, axiSym);
        
    // Rod:
    mats.push_back(new Material("rod", dom, new VinetRose("tantalum"), new ZerilliArmstrong("tantalum")));
    ICBlock(*mats[0], dom, Rect(DV(-0.00381,-0.0254),DV(0.00381,0.0254)), DV3(0,-190,0));

  } else if(icType == BerylliumShell) {

    if(axiSym) {
      cerr << "Initial conditions not yet done for axisymmetric shell collapse - TODO!\n";
      exit(1);
    }

    // Beryllium shell collapse:    
    // Produces shell which collapses to inner radius 0.05, outer radius 0.0781
    tf = 200e-6;
    outFreq = 25e-6;
    dom = Domain(Ni, Ni, 3, 0.11);
    reinitFreq = 40;

    outdir = getDirName("Beryllium", dom, ifaceType, axiSym);

    // Shell:
    mats.push_back(new Material("shell", dom, new VinetRose("beryllium"), new PerfectPlastic("beryllium"))); 
    ICShell(*mats[0], dom, DV(0,0,0), 0.05/0.625, 0.10, 0.05);
  } else
    if(icType == Radial) {

      // Radial "impact" problem
      tf = 1.0e-6;                 // Final time
      outFreq = 1e-7;            // Output frequency    
      dom = Domain(Ni, Ni, 3, 0.022); // Domain

      outdir = getDirName("Radial2D", dom, ifaceType, axiSym);

      // Centre
      mats.push_back(new Material("inner", dom, new Romenski("copper"), new PerfectPlastic("copper")));
      ICShell2(*mats[0], dom, DV(0.011,0.011), 0.001, 0.005, 0);	     
    
      // Outer
      mats.push_back(new Material("outer", dom, new Romenski("copper"), new PerfectPlastic("copper")));
      ICShell2(*mats[1], dom, DV(0.011,0.011), 0.005, 0.01, 1000);

    } else 
      {
	cerr << "Unknown problem type\n";
	exit(1);   
      }


  int ret = system(("mkdir -p " + outdir).c_str());
  cerr << "Writing results to \"" << outdir << "\" (" << ret << ")\n";

  int step = 0;
  int outStep = 0;

  for(double t=0; t<tf; ++step) {
    // Reinitialize level sets if requested
    if(step%reinitFreq == 0) {
      cerr << "Reinitializing level sets\n";

      for(MaterialIterator mit=mats.begin(); mit!=mats.end(); ++mit) {
	Material& mat=**mit;
	BCs(mat, dom);
	mat.sys.interpMethod = interpMethod;
	mat.ls.extrapolateToAdjacent(mat.sol, mat.sys);
	mat.ls.reinitialize();
	mat.ls.extrapolateToAdjacent(mat.sol, mat.sys);
      }
      
    }    
    
    for(MaterialIterator mit=mats.begin(); mit!=mats.end(); ++mit) {
      Material& mat=**mit;
      BCs(mat, dom);
    }

    // Output current state if requested
    if(t>=outStep*outFreq||step%outInterval==0) 
      {
	extrapolateGhosts(mats);
	cerr << "Saving results to output file " << outStep << "...";
	outputSilo(outdir, dom, mats, outStep, t);
	cerr << " done\n";
	++outStep;
      }

    // Compute timestep
    double dt;
    if(fixed_dt<=0.0) {
      double mindt = 1E100;
      for(MaterialIterator mit=mats.begin(); mit!=mats.end(); ++mit) {
	mindt = min(mindt, getMinDt(**mit, dom));
      }
      dt = CFL * mindt;
    } else {
      dt = fixed_dt;
    }       

    if(t<outStep*outFreq && t+dt>outStep*outFreq) {
      dt = outStep*outFreq-t;
    }
    
    setGhosts(mats, dom, dt, ifaceType);
    //    if(mats.size()==1) {
    //      syncVelocity(*(mats[0]), dom);      
    //    }

    for(MaterialIterator mit=mats.begin(); mit!=mats.end(); ++mit) {
      Material& mat=**mit;
      BCs(mat, dom);
      if(icType==BerylliumShell) outputRadius(outdir, mat, dom, t);    
    }

    if(step==0) initialIgrals(mats, dom, axiSym);
    outputIgrals(outdir, mats, dom, t, axiSym);

    if(!stations.empty()) 
      {      
	outputStations(outdir, mats, dom, t, dt);
      }
    
    // Status output
    cout << "["<<outStep<<"] " << step << ": t = " << t << ", dt = " << dt << endl;

    // Compute updates to each material in turn
    for(MaterialIterator mit=mats.begin(); mit!=mats.end(); ++mit) {
      Material& mat = **mit;
      elasticSolve(mat, dom, dt, axiSym);
      plasticSolve(mat, dom);     
      BCs(mat,dom);
      mat.ls.extrapolateToAdjacent(mat.sol, mat.sys);
      BCs(mat,dom);
      mat.ls.update(dt);
      
    }
    
    t+=dt;
  }

  // Output at final time
  extrapolateGhosts(mats);
  cerr << "Saving results to output file " << outStep << "...";
  outputSilo(outdir, dom, mats, outStep, tf);
  cerr << " done\n";

  return 0;
}
