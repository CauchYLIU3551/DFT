#include "../include/DFT.h"

double _u_(const double * p)
{
  return sin(PI*p[0]) * sin(2*PI*p[1]);
}

double _f_(const double * p)
{ // have to be edited as the random initialization for the initial
  // guess if the electric density.
  return 0.0;
}

void DFT::addBoundaryCondition_Hartree()
{ // This is for the first test example with Dirchlet condition!!!;
  // This is just a temporary function which have to be edited to
  // solve the entire K-S problem; 
  // initialize the finite element functional u_h;
  //u_h->reinit(*fem_space);

  V_hartree = new FEMFunction<double,DIM>(*fem_space);
  
  
  BoundaryFunction<double, DIM> boundary(BoundaryConditionInfo::DIRICHLET, 1, &_f_);
  BoundaryConditionAdmin<double, DIM> boundary_admin(*fem_space);
  boundary_admin.add(boundary);
  boundary_admin.apply(*stiff_Hartree, *V_hartree, *rhs);
}
