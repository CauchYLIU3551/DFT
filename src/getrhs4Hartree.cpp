#include "../include/DFT.h"

void DFT::getRHS_Hartree()
{
  std::cout<<"Begin to compute the rhs vector for Hartree...\n";
  rhs->reinit(fem_space->n_dof(), false);
  FEMSpace<double,DIM>::ElementIterator the_ele = fem_space->beginElement();
  FEMSpace<double,DIM>::ElementIterator end_ele = fem_space->endElement();
  for(;the_ele != end_ele;++ the_ele)
    {
      double vol = the_ele->templateElement().volume();
      const QuadratureInfo<DIM>& qi = the_ele->findQuadratureInfo(3);
      u_int n_q_pnt = qi.n_quadraturePoint();
      std::vector<double> jac = the_ele->local_to_global_jacobian(qi.quadraturePoint());
      std::vector<AFEPack::Point<DIM> > q_pnt = the_ele->local_to_global(qi.quadraturePoint());
      std::vector<std::vector<double> > bas_val = the_ele->basis_function_value(q_pnt);

      std::vector<double> rho_val = rho->value(q_pnt,*the_ele);

      const std::vector<int>& ele_dof = the_ele->dof();
      u_int n_ele_dof = ele_dof.size();

      for (u_int l = 0;l < n_q_pnt;++ l)
	{
	  double Jxw = vol*qi.weight(l)*jac[l];
	  for (u_int i = 0; i < n_ele_dof;++ i)
	    {
	      (*rhs)(ele_dof[i]) += -4.0*PI*bas_val[i][l]*rho_val[l];
	    }
	}
    }
  std::cout<<"The size of rhs is "<<rhs->size()<<std::endl;
}
