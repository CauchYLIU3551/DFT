#include "../include/DFT.h"
//#include "../include/Matrix.h"

///////////////////////////////////////////
// This function computes the matrix A for
// generalized eigen problem of KS
void DFT::buildMatrixA4KS()
{
  std::cout << "build a new stiff_matrix for KS equation..."<<std::endl;
  stiff_KS = new MatrixA(*fem_space, *V_hartree);
  stiff_KS->algebricAccuracy() = 3;
  stiff_KS->build();
}

///////////////////////////////////////////
// This function computes the matrix M for
// generalized eigen problem of KS
void DFT::buildMatrixM4KS()
{
  std::cout << "build a new mass_matrix for KS equation..." << std::endl;
  mass_KS = new MassMatrix<DIM,double>(*fem_space);
  mass_KS->algebricAccuracy() = 3;
  mass_KS->build();
}

void MatrixA::getElementMatrix(const Element<double,DIM>& ele0,
			       const Element<double,DIM>& ele1,
			       const ActiveElementPairIterator<DIM>::State state)
{
  double vol = ele0.templateElement().volume();

  const QuadratureInfo<DIM>& quad_info =
    ele0.findQuadratureInfo(algebricAccuracy());
  std::vector<double> jac =
    ele0.local_to_global_jacobian(quad_info.quadraturePoint());

  // number of quadrature points
  int n_q_pnt = quad_info.n_quadraturePoint();

  // transform the points from template ele to global ele in mesh
  std::vector<AFEPack::Point<DIM> > q_pnt =
    ele0.local_to_global(quad_info.quadraturePoint());

  // compute the gradient of basis functions
  std::vector<std::vector<std::vector<double>>> basis_grad =
    ele0.basis_function_gradient(q_pnt);

  // compute the value of basis functions
  std::vector<std::vector<double>> bas_val = ele0.basis_function_value(q_pnt);
  
  int n_ele_dof = ele0.dof().size();

  // sum the quadrature points
  for (int l = 0;l < n_q_pnt;l ++)
    {
      double Jxw = quad_info.weight(l) * jac[l] * vol;

      // compute the V_Hartree at quadrature points;
      /////
      double V_val = V->value(q_pnt[l],ele0);
      /////

      for (int j = 0;j < n_ele_dof;j ++)
	{
	  for (int k = 0;k < n_ele_dof;k ++)
	    {
	      elementMatrix(j,k) += Jxw*(0.5*innerProduct(basis_grad[j][l],basis_grad[k][l])
					 - kappa * V_val * bas_val[j][l]*bas_val[k][l]);
	      
	    }
	}
      
    }
}
