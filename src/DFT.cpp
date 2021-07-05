#include "../include/DFT.h"

double _u_(const double * p)
{
  return sin(PI*p[0]) * sin(2*PI*p[1]);
}

double _f_(const double * p)
{ // have to be edited as the random initialization for the initial
  // guess if the electric density.
  return 5*PI*PI*_u_(p);
}


//////////////////////////////////////////
DFT::DFT()
{
}
DFT::DFT(HGeometryTree<DIM>* _h_tree,
	      IrregularMesh<DIM>* _ir_mesh)
{
  h_tree = _h_tree;
  irregular_mesh = _ir_mesh;
  
}

void DFT::initialize()
{
  template_geometry.readData("tetrahedron.tmp_geo");
  coord_transform.readData("tetrahedron.crd_trs");
  template_dof.reinit(template_geometry);
  template_dof.readData("tetrahedron.1.tmp_dof");
  basis_function.reinit(template_dof);
  basis_function.readData("tetrahedron.1.bas_fun");

  template_element.resize(1);
  template_element[0].reinit(template_geometry,
			     template_dof,
			     coord_transform,
			     basis_function);
}

void DFT::buildspace()
{
  RegularMesh<DIM>& regular_mesh = irregular_mesh->regularMesh();
  fem_space->reinit(regular_mesh, template_element);
  int n_element = regular_mesh.n_geometry(DIM);
  fem_space->element().resize(n_element);
  for (u_int i = 0;i < n_element;++ i)
    {
      fem_space->element(i).reinit(*fem_space, i, 0);
    }
  
  fem_space->buildElement();
  fem_space->buildDof();
  fem_space->buildDofBoundaryMark();

  std::cout<<"fem space is built..."<<std::endl;
}

void DFT::buildHartreeMatrix()
{
  if(stiff_Hartree != NULL){
    std::cout << "deleting the current stiff_matrix" << std::endl;
    delete stiff_Hartree;

    stiff_Hartree = new StiffMatrix<DIM, double>(*fem_space);
  }

  std::cout << "build a new stiff_matrix for Hatree potential" << std::endl;
  //  stiff_matrix->reinit(*fem_space);
  stiff_Hartree->algebricAccuracy() = ACC;
  stiff_Hartree->build();
}

void DFT::getRHS_Hartree()
{
  //Vector<double> f_h;
  Operator::L2Discretize(&_f_, *fem_space, *rhs, 3);
}

void DFT::addBoundaryCondition_Hartree()
{
  // This is just a temporary function which have to be edited to
  // solve the entire K-S problem; 
  // initialize the finite element functional u_h;
  u_h->reinit(*fem_space);
  
  BoundaryFunction<double, DIM> boundary(BoundaryConditionInfo::DIRICHLET, 1, &_u_);
  BoundaryConditionAdmin<double, DIM> boundary_admin(*fem_space);
  boundary_admin.add(boundary);
  boundary_admin.apply(*stiff_Hartree, *u_h, *rhs);
}

void DFT::getHartree()
{
  AMGSolver solver(*stiff_Hartree);
  solver.solve(*u_h, *rhs);
  u_h->writeOpenDXData("u.dx");
}