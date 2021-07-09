#ifndef _DFT_h_
#define _DFT_h_

#include <lac/sparse_matrix.h>
#include <lac/sparsity_pattern.h>
#include <iostream>
#include <math.h>
#include <vector>
#include <string>

#include <EigenSolver/EigenSolver.h>

#include <AFEPack/BilinearOperator.h>
#include <AFEPack/Functional.h>
#include <AFEPack/Operator.h>
#include <AFEPack/EasyMesh.h>
#include <AFEPack/TemplateElement.h>
#include <AFEPack/FEMSpace.h>
#include <AFEPack/Geometry.h>
#include <AFEPack/BoundaryCondition.h>
#include <AFEPack/HGeometry.h>
#include <AFEPack/DGFEMSpace.h>

//#include "./Matrix.h"

#define DIM 3
#define PI 4.0*atan(1)
#define ACC 2
#define kappa 1.0

////////////////////////////////////////////////////////////////////
// This class is used to build the stiffmatrix for KS
class MatrixA : public L2InnerProduct<DIM,double>
{
  //friend class DFT;
 public:
  MatrixA(FEMSpace<double,DIM> &sp, FEMFunction<double,DIM>& v):
  L2InnerProduct<DIM,double>(sp,sp), V(&v){}
  virtual ~MatrixA(){};
 public:
  virtual void
  getElementMatrix(const Element<double,DIM>& ele0,
		   const Element<double,DIM>& ele1,
		   const ActiveElementPairIterator<DIM>::State state);
 private:
  FEMFunction<double,DIM> * V;
};
////////////////////////////////////////////////////////////////////


class DFT
{
  //friend class MatrixA;
 public:
  DFT();
  DFT(HGeometryTree<DIM>* _h_tree,
	   IrregularMesh<DIM>* _ir_mesh);
  
  /**
   * This function reads the input .mesh file to get the mesh info.
   */
  void readMesh(const std::string& fileName,
		const int& refine_times);

  void getBoundaryIndex();
  
  /**
   * This function is used to initialization, including read mesh 
   * data, build the template elent.
   */
  void initialize();

  /**
   * This function initialize the electric density rho to start 
   * self-consistent iteration.
   */
  void initializeRho();

  /**
   * This function is used to build the finite element space for KS
   * equations to get the wavefunctions.
   */
  void buildspace();

  /**
   * This function is used to build the finite element space for 
   * poisson equation to get the Hartree potential to complete 
   * self-consistent iteration.
   */
  //void buildspaceHartree();

  /**
   * This function helps to add boundary mark to the info of the 
   * points at boundary.
   */
  void addBoundMark();

  /**
   * This function builds the stiff matrix of the poisson equation
   * of Hartree potential to build the linear system.
   */
  void buildHartreeMatrix();

  /**
   * This function computes the rhs vector corresponding to the 
   * poisson equation of Hartree potential.
   */
  void getRHS_Hartree();

  /**
   * Apply the boundary condition of the poisson equation into the 
   * stiff matrix and rhs vector.
   * In fact, it is not accurate enough to set the boundary as 
   * Dirchlet Condition with value 0. So it has to be updated in
   * the following tests.
   */
  void addBoundaryCondition_Hartree();
  
  /**
   * This function is used to compute the Hartree potential.
   * In fact, It just solve a linear system of a poisson equation.
   */
  void getHartree();

  //////////////////////////////////////////////////////////////////
  // Begin to build the generalized eigenvalue problems to solve the
  // KS equation.
  // The form of the generalized eigenvalue problem should be
  // Ax = e * Mx;
  /**
   * This function builds the stiff matrix A of the generalized 
   * eigenvalue problem equation of the KS equation.
   */
  void buildMatrixA4KS();

  /**
   * This function builds the mass matrix M of the eigen problem;
   */
  void buildMatrixM4KS();

  /**
   * This function applies the boundary condition into the matrices 
   * A and M in the equation A x = e * M x;
   */
  void addBoundaryCondition_KS();
  
  /**
   * This function can normalize the wave function to satisfy the 
   * normalization condition of the wave functions.
   */
  void normalize(FEMFunction<double,DIM>& P); 

  void printMesh();

  /**
   * This function solve the poisson equation to get V_Hartree 
   */
  void solveHartree();
  
 private:
  //int DIM;
  //Mesh<DIM> mesh;
  HGeometryTree<DIM> * h_tree;
  IrregularMesh<DIM> * irregular_mesh;
  IrregularMesh<DIM> * old_irregular_mesh;
  //IrregularMesh<DIM> * old_irregular_mesh;

  //DGFEMSpace<double, DIM> * fem_space;
  //DGFEMSPace<double, DIM> * old_fem_space;

  StiffMatrix<DIM, double> * stiff_Hartree;
  //StiffMatrix<DIM, double> * stiff_matrix;
  MatrixA * stiff_KS;
  MassMatrix<DIM, double> * mass_KS;
  
  // for regular template 
  TemplateGeometry<DIM> template_geometry;
  CoordTransform<DIM,DIM> coord_transform;
  TemplateDOF<DIM> template_dof;
  BasisFunctionAdmin<double,DIM,DIM> basis_function;
  UnitOutNormal<DIM> unit_out_normal;

  /// for twin_template
  TemplateGeometry<DIM> twin_template_geometry;
  CoordTransform<DIM, DIM> twin_coord_transform;
  TemplateDOF<DIM> twin_template_dof;
  BasisFunctionAdmin<double, DIM, DIM> twin_basis_function;
  UnitOutNormal<DIM> twin_unit_out_normal;

  /// for four_template
  TemplateGeometry<DIM> four_template_geometry;
  CoordTransform<DIM, DIM> four_coord_transform;
  TemplateDOF<DIM> four_template_dof;
  BasisFunctionAdmin<double, DIM, DIM> four_basis_function;
  UnitOutNormal<DIM> four_unit_out_normal;

  /// for boundary 
  TemplateGeometry<DIM-1> triangle_template_geometry;
  CoordTransform<DIM-1,DIM> triangle_to3d_coord_transform;

  TemplateGeometry<DIM-1> twin_triangle_template_geometry;
  CoordTransform<DIM-1,DIM> twin_triangle_to3d_coord_transform;
  
  std::vector<TemplateElement<double,DIM,DIM>> template_element;
  std::vector<TemplateDGElement<DIM-1,DIM> > edge_template_element;
  
  
  DGFEMSpace<double,DIM> * fem_space;
  DGFEMSpace<double,DIM> * old_fem_space;
  
  FEMFunction<double,DIM> * u_h;
  FEMFunction<double,DIM> * V_hartree;
  
  Vector<double> * rhs;
  //Vector<double> * rho;
  FEMFunction<double,DIM> * phi;
  FEMFunction<double,DIM> * rho;

  std::vector<int> boundaryDOFIndex;
  std::vector<int> internalDOFIndex;
  std::vector<int> boundaryDGEleIndex;
  
};


#endif

