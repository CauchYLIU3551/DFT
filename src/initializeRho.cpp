#include "../include/DFT.h"

double inner(std::vector<double>a, std::vector<double>b)
{
  int n = a.size();
  double sum = 0;
  for(int i = 0;i < n;i ++)
    {
      sum += a[i] * b[i];
    }
  return sum;
}

std::vector<double> proj(std::vector<double> u, std::vector<double> v)
{
  double delta = 0;
  delta = inner(u,v) / inner(u,u);
  for(int i=0;i<u.size();i++)
    {
      v[i] = delta*u[i];
    }
  return v;
}
////////////////////////////////////////////////////
// This function GS can compute the Gram-Schmdit process of the matrix
// X to get the orthonormal matrix corresponding to X;
void GS(std::vector<std::vector<double>>& X)
{
  std::vector<std::vector<double>> tempX;
  tempX = X;
  for(u_int i = 1;i < tempX.size();i++)
    {
      for(u_int j = 0;j < i;j ++)
	{
	  std::vector<double> temp_proj;
	  temp_proj = proj(X[j], tempX[i]);
	  for (u_int k = 0;k < tempX[0].size(); k++)
	    {
	      X[i][k] = X[i][k] - temp_proj[k];
	    }
	}
    }

  for (u_int k = 0;k < X.size();k ++)
    {
      double norm = 0;
      norm = sqrt(inner(X[k], X[k]));
      for (u_int i = 0;i < X[0].size();i ++)
	{
	  X[k][i] = X[k][i]/norm;
	}
    }
}

/*
void DFT::initializeRho()
{
  const int& n_dof = fem_space->n_dof();
  //delete rho;
  std::cout<<"give two random wave function phi...\n";
  rho = new Vector<double>(n_dof);
  rho->reinit(n_dof, false);
  std::vector<double> phi_1(n_dof), phi_2(n_dof);

  // Have a try of rho = (phi_1)^2 + (phi_2)^2;

  for (u_int i = 0;i < n_dof;++ i)
    {
      phi_1[i] = ((double) rand() / (RAND_MAX));
      phi_2[i] = ((double) rand() / (RAND_MAX));
    }
  
  // Normalize these two wave function
  std::vector<std::vector<double>> phi;
  // phi.resize(2, false);
  phi.push_back(phi_1);
  phi.push_back(phi_2);
  // std::cout<<"begin GS process...\n";
  GS(phi);
  //  std::cout<<"After GS process...\n";
  phi_1 = phi[0];
  phi_2 = phi[1];

  //////////////////
  // Using phi1 and phi2 to get the electric density rho;
  //std::cout<<"begin compute rho...\n";
  //std::cout<<"The size of rho is : "<<rho->size()<<std::endl;
  //std::cout<<"The n_dof is "<<n_dof<<std::endl;
  for (u_int i = 0;i < n_dof;i ++)
    { //std::cout<<"This is the "<<i <<"-th iterations... \n";
      //std::cout<<"The phi_1["<<i<<"] is "<<phi_1[i]<<"\n";
      //std::cout<<"The phi_2["<<i<<"] is "<<phi_2[i]<<"\n";
      (*rho)[i] = phi_1[i]*phi_1[i] + phi_2[i]*phi_2[i];
      //std::cout<<"The rho["<<i<<"] is "<<rho[i]<<"\n\n";
    }
  std::cout<<"Finish compute rho...\n";
}
*/

void DFT::initializePhi()
{
  phi = new FEMFunction<double,DIM>(*fem_space);
  old_phi = new FEMFunction<double,DIM>(*fem_space);
  rho = new FEMFunction<double,DIM>(*fem_space);

  for (u_int i = 0;i < rho->size();++ i)
    {
      (*phi)[i] = ((double) rand() / (RAND_MAX));
      //  std::cout<<(*phi)[i]<<" ";
    }
  // std::cout<<"The initial rho begin normalize...\n";
  normalize(*phi);
  //rho->reinit(*fem_space);
  //std::cout<<"The size of rho is "<<rho->size()<<std::endl;

  ////////////////////////
  /*
  FEMSpace<double, DIM>::ElementIterator
    the_ele = fem_space->beginElement(),
    end_ele = fem_space->endElement();y
  for(;the_ele != end_ele;++ the_ele)
    {
      double vol = the_ele->templateElement().volume();
      const QuadratureInfo<DIM>& qi = the_ele->findQuadratureInfo(3);
      u_int n_q_pnt = qi.n_quadraturePoint();
      std::vector<double> jac = the_ele->local_to_global_jacobian(qi.quadraturePoint());
      std::vector<AFEPack::Point<DIM>> q_pnt = the_ele->local_to_global(qi.quadraturePoint());
      std::vector<double> rho_val = rho->value(q_pnt, *the_ele);
      // compute the integral of the wavefunction to normalize it.
      for (u_int l = 0;l < n_q_pnt;l ++)
	{
	  double Jxw = qi.weight(l)*jac[l]*vol;
	  norm += Jxw * rho_val[l] * rho_val[l];
	}
    }
  std::cout<<"This is the initial rho at the vertice of mesh...\n";
  for (u_int i = 0;i < rho->size();++ i)
    {
      (*rho)[i] /=  sqrt(norm);
      std::cout<<(*rho)[i]<<" ";
    }
  
  */
  //////////////////////////////

  /*
  for (u_int i = 0;i < rho->size();i ++)
    {
      (*rho)[i] = 2 * (*phi)[i] * (*phi)[i];
      // std::cout<<(*rho)[i]<<" ";
    }*/

  
  // std::cout<<"size of rho is "<<rho->size()<<std::endl;
}
