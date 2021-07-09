#include "../include/DFT.h"

void DFT::solveHartree()
{
  AMGSolver solver(*stiff_Hartree);
  solver.solve(*V_hartree, *rhs, 1.0e-08, 50);
}

void DFT::solveKS()
{
  ///////////////////////////////////////////////
  std::vector<double> tempphi(1,0);
  // In this step-1 case ,it just needs to compute one wavefunction;
  std::vector<std::vector<double>> tmp_sol(phi->size(), tempphi);

  for(int i=0;i<phi->size();i++)
    {
      //tempphi[i] = (*phi)(i);
      //tmp_sol[i][0] = ((double) rand() / (RAND_MAX));
      tmp_sol[i][0] = (*phi)(i);
    }
  ///////////////////////////////////////////////
  //

  ///////////////////////////////////////////////////////////
  //delete old_phi;
  /*
  old_phi = new FEMFunction<double,DIM>(*fem_space);
  old_phi = phi;
  delete phi;
  phi = new FEMFunction<double,DIM>(*fem_space);
*/
  for(int i = 0;i<phi->size();i++)
  {
    (*old_phi)(i) = (*phi)(i);
  }
  ////////////////////////////////////////////////////////////
  //std::cout<<"OldPHI::::: "<<(*old_phi)(0)<<"\n";
  lambda_e.clear();
  lambda_e.resize(phi->size(), 0.0);
  EigenSolver solve(*stiff_KS,*mass_KS);
  solve.LOBPCGSolve(tmp_sol, lambda_e, 1, 100, 1.e-5);
  std::cout<<"lambda is "<< lambda_e[0];
  ///////////////////////////////////////////////
  //std::cout<<"PHI::::: "<<"\n";
  for(int i = 0;i < phi->size();i ++)
    {
      (*phi)(i) = tmp_sol[0][i];
    //  std::cout<<tmp_sol[0][i]<<"\n";
    }


  std::filebuf fb;
  fb.open("stiff_matrix.txt", std::ios::out);
  std::ostream os(&fb);
  stiff_KS->print_formatted(os, 3, true, 0,"0.0",1);
  fb.close();
  std::filebuf fb2;
  fb2.open("mass_matrix.txt", std::ios::out);
  std::ostream os2(&fb2);
  mass_KS->print_formatted(os2, 3, true, 0,"0.0",1);
  fb2.close();

  std::ofstream sparsematrix ("stiff_matrix.1");
  stiff_KS->print(sparsematrix);



  ///////////////////////////////////////////////
}

// compute infi norm of |phi - phi_old|
void DFT::get_res4psi(double& res)
{
  //double norm = 0;
  res = 0.0;
  for (u_int i = 0;i < phi->size();i ++)
    {
      if (fabs((*phi)(i) - (*old_phi)(i)) > res)
	{
	  res = fabs((*phi)(i) - (*old_phi)(i));
	}
    }
}

// computing rho_new = (1-alpha)*old_phi + alpha*phi_1;
void DFT::mixing()
{
  for (u_int i = 0;i < phi->size();i ++)
    {
      (*phi)(i) = (1.0 - alpha)*((*old_phi)(i)) + alpha*((*phi)(i));
    }
}
