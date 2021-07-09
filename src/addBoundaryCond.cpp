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

void DFT::addBoundaryCondition_KS()
{
  const int & n_dof = fem_space->n_dof();
  RegularMesh<DIM>& regular_mesh = irregular_mesh->regularMesh();

  SparsityPattern& sp = stiff_KS->getSparsityPattern();
  const size_t * row_start = sp.get_rowstart_indices();
  const unsigned int * column = sp.get_column_numbers();

  double bnd_value = 0.;
  int k = 0;
  const int& n_boundaryDOFIndex = boundaryDOFIndex.size();
  for(int i = 0;i < n_boundaryDOFIndex;++ i)
    {
      const int& dof_idx = boundaryDOFIndex[i];
      
 
      // the boundary value is set as 0;
      //(*rhs)(dof_idx) = stiff_KS->global_entry(row_start[dof_idx]) * bnd_value;
      for (int j = row_start[dof_idx] + 1;j < row_start[dof_idx + 1];++ j)
	{
	  stiff_KS->global_entry(j) = 0.;
	  mass_KS->global_entry(j) = 0.;
	  k = column[j];
	  const unsigned int *p = std::find(&column[row_start[k] + 1], &column[row_start[k+1]], j);
	  if(p != &column[row_start[k+1]])
	    {
	      //(*rhs)(*p) -= stiff_KS->global_entry(*p) * bnd_value;
	      stiff_KS->global_entry(*p) = 0.;
	      mass_KS->global_entry(*p) = 0.; 
	    }
	}
    } 
}
