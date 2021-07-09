#include "../include/DFT.h"

void DFT::normalize(FEMFunction<double, DIM>& P)
{
  std::cout<<"Begin normalize...\n";
  double norm = 0;
  //rho->reinit(*fem_space);
  //std::cout<<"The size of rho is "<<rho->size()<<std::endl;
  FEMSpace<double, DIM>::ElementIterator
    the_ele = fem_space->beginElement(),
    end_ele = fem_space->endElement();
  for(;the_ele != end_ele;++ the_ele)
    {
      double vol = the_ele->templateElement().volume();
      const QuadratureInfo<DIM>& qi = the_ele->findQuadratureInfo(3);
      u_int n_q_pnt = qi.n_quadraturePoint();
      std::vector<double> jac = the_ele->local_to_global_jacobian(qi.quadraturePoint());
      std::vector<AFEPack::Point<DIM>> q_pnt = the_ele->local_to_global(qi.quadraturePoint());
      std::vector<double> rho_val = P.value(q_pnt, *the_ele);
      // compute the integral of the wavefunction to normalize it.
      for (u_int l = 0;l < n_q_pnt;l ++)
	{
	  double Jxw = qi.weight(l)*jac[l]*vol;
	  norm += Jxw * rho_val[l] * rho_val[l];
	}
    }
 // std::cout<<"This is the initial rho at the vertice of mesh...\n";
  for (u_int i = 0;i < P.size();++ i)
    {
      P[i] /=  sqrt(norm);
   //   std::cout<<P[i]<<" ";
    }
  //std::cout<<std::endl;

}

void DFT::normalize_phi()
{
  normalize(*phi);
}
