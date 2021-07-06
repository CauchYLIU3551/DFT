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
  rhs = new Vector<double>();
  fem_space = new DGFEMSpace<double, DIM>();
  stiff_matrix = new StiffMatrix<DIM, double>();
}
DFT::DFT(HGeometryTree<DIM>* _h_tree,
	      IrregularMesh<DIM>* _ir_mesh)
{
  h_tree = _h_tree;
  irregular_mesh = _ir_mesh;

  rhs = new Vector<double>();
  fem_space = new DGFEMSpace<double, DIM>();
  stiff_matrix = new StiffMatrix<DIM, double>();
}


/*
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
*/
void DFT::initialize()
{
  template_geometry.readData("tetrahedron.tmp_geo");
  coord_transform.readData("tetrahedron.crd_trs");
  unit_out_normal.readData("tetrahedron.out_nrm");
  template_dof.reinit(template_geometry);
  template_dof.readData("tetrahedron.1.tmp_dof");
  basis_function.reinit(template_dof);
  basis_function.readData("tetrahedron.1.bas_fun");

  twin_template_geometry.readData("twin_tetrahedron.tmp_geo");
  twin_coord_transform.readData("twin_tetrahedron.crd_trs");
  twin_unit_out_normal.readData("twin_tetrahedron.out_nrm");
  twin_template_dof.reinit(twin_template_geometry);
  twin_template_dof.readData("twin_tetrahedron.1.tmp_dof");
  twin_basis_function.reinit(twin_template_dof);
  twin_basis_function.readData("twin_tetrahedron.1.bas_fun");

  four_template_geometry.readData("four_tetrahedron.tmp_geo");
  four_coord_transform.readData("four_tetrahedron.crd_trs");
  four_unit_out_normal.readData("four_tetrahedron.out_nrm");
  four_template_dof.reinit(four_template_geometry);
  four_template_dof.readData("four_tetrahedron.1.tmp_dof");
  four_basis_function.reinit(four_template_dof);
  four_basis_function.readData("four_tetrahedron.1.bas_fun");

  template_element.resize(3);
  template_element[0].reinit(template_geometry,
			     template_dof,
			     coord_transform,
			     basis_function,
			     unit_out_normal);

  template_element[1].reinit(twin_template_geometry,
			     twin_template_dof,
			     twin_coord_transform,
			     twin_basis_function,
			     twin_unit_out_normal);

  template_element[2].reinit(four_template_geometry,
			     four_template_dof,
			     four_coord_transform,
			     four_basis_function,
			     four_unit_out_normal);

  triangle_template_geometry.readData("triangle.tmp_geo");
  triangle_to3d_coord_transform.readData("triangle.to3d.crd_trs");

  twin_triangle_template_geometry.readData("twin_triangle.tmp_geo");
  twin_triangle_to3d_coord_transform.readData("twin_triangle.to3d.crd_trs");

  edge_template_element.resize(2);
  edge_template_element[0].reinit(triangle_template_geometry,
				  triangle_to3d_coord_transform);
  edge_template_element[1].reinit(twin_triangle_template_geometry,
				  twin_triangle_to3d_coord_transform);
}

void DFT::buildspace()
{
  //std::cout<<"Begin to build fem space!\n";
  RegularMesh<DIM>& regular_mesh = irregular_mesh->regularMesh();
  //std::cout<<"build a regular mesh!\n";
  fem_space->reinit(regular_mesh, template_element);
  // std::cout<<"reinitialzie the fem space!\n";
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
  stiff_Hartree = new StiffMatrix<DIM, double>(*fem_space);
  stiff_Hartree->algebricAccuracy() = ACC;
  stiff_Hartree->build();
}

void DFT::getRHS_Hartree()
{
  //Vector<double> f_h;
  //Operator::L2Discretize(&_f_, *fem_space, *rhs, 3);
  const int& n_dof = fem_space->n_dof();

  // The initial hartree rhs can be a random vector of the
  // electric density rho;
  rhs->reinit(n_dof, false);
  for(int i = 0;i < n_dof;++ i)
    {
      (*rhs)(i) = 0.;
    }
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
