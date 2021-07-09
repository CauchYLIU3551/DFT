#include "../include/DFT.h"

//////////////////////////////////////////
DFT::DFT()
{
  rhs = new Vector<double>();
  fem_space = new DGFEMSpace<double, DIM>();
  //stiff_matrix = new StiffMatrix<DIM, double>();
}
DFT::DFT(HGeometryTree<DIM>* _h_tree,
	      IrregularMesh<DIM>* _ir_mesh)
{
  h_tree = _h_tree;
  irregular_mesh = _ir_mesh;

  rhs = new Vector<double>();
  fem_space = new DGFEMSpace<double, DIM>();
  //stiff_matrix = new StiffMatrix<DIM, double>();
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

  u_int n_edge = regular_mesh.n_geometry(DIM-1);
  fem_space->dgElement().resize(n_edge);
  for(u_int i = 0;i < n_edge;i ++)
    {
      fem_space->dgElement(i).reinit(*fem_space,i,0);
    }
  
  std::cout<<"fem space is built..."<<std::endl;
}


void DFT::getHartree()
{
  AMGSolver solver(*stiff_Hartree);
  solver.solve(*u_h, *rhs);
  u_h->writeOpenDXData("u.dx");
}

void DFT::printMesh()
{
  irregular_mesh->regularMesh().writeOpenDXData("D.dx");
}
