#include "../include/DFT.h"

////////////////////////////////////////////////////////////////////
// In fact, this boundary mark is used for the domain as a cube
// in [-0.5,0.5]x[-0.5,0.5]x[-0.5,0.5]
bool onBoundary(const double * p)
{
#if 1 /// for cube
  if (fabs(p[0] - 0.5) < 1.0e-03 || fabs(p[0] + 0.5) < 1.0e-03)
    {
      return true;
    }
  else if (fabs(p[1] - 0.5) < 1.0e-03 || fabs(p[1] + 0.5) < 1.0e-03)
    {
      return true;
    }
  else if (fabs(p[2] - 0.5) < 1.0e-03 || fabs(p[2] + 0.5) < 1.0e-03)
    {
      return true;
    }
  else
    {
      return false;
    }
#else /// for sphere
  
  if(fabs(sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]) - radius_sphere) < 1.0e-03)
    {
      std::cout << "The length is " << sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]) << ", and the boundary mark is " << 1 << std::endl;    
      return true;
    }
  else
    {
      std::cout << "The length is " << sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]) << ", and the boundary mark is " << 0 << std::endl;
      return false;
    }
#endif
}


void DFT::addBoundMark()
{
  RegularMesh<DIM>& mesh = irregular_mesh->regularMesh();

  const int& n_pnt = mesh.n_point();

  for(int i = 0;i < n_pnt;++ i)
    {
      GeometryBM& pnt_geo = mesh.geometry(0, i);
      AFEPack::Point<DIM>& pnt = mesh.point(pnt_geo.vertex()[0]);

      if(onBoundary(pnt))
	{
	  pnt_geo.boundaryMark() = 1;
	}
      else
	{
	  pnt_geo.boundaryMark() = 0;
	}
    }


  std::cout << "boundary points are done..." << std::endl;
  
  const int& n_line = mesh.n_geometry(1);

  for(int i = 0;i < n_line;++ i)
    {
      GeometryBM& line_geo = mesh.geometry(1, i);
      line_geo.boundaryMark() = 1;/// assume it is boundary
      const int& n_vtx = line_geo.n_vertex();

      for(int j = 0;j < n_vtx;++ j)
	{
	  if(mesh.geometry(0, line_geo.vertex()[j]).boundaryMark() == 0)
	    {
	      line_geo.boundaryMark() = 0;
	      break;
	    }
	}
    }

  std::cout << "boundary lines are done..." << std::endl;
  
  const int& n_tri = mesh.n_geometry(2);

  for(int i = 0;i < n_tri;++ i)
    {
      GeometryBM& tri_geo = mesh.geometry(2, i);

      tri_geo.boundaryMark() = 0;

      AFEPack::Point<DIM> bc(DIM);
      bc[0] = 0; bc[1] = 0; bc[2] = 0;
    
      const int& n_vtx = tri_geo.n_vertex();
      for(int j = 0;j < n_vtx;++ j)
	{
	  AFEPack::Point<DIM>& vtx1 = mesh.point(tri_geo.vertex()[j]);
	  for(int k = 0;k < DIM;++ k)
	    {
	      bc[k] += vtx1[k];
	    }
	}

      for(int k = 0;k < DIM;++ k)
	{
	  bc[k] /= n_vtx;
	}
      
      if(fabs(bc[0] + 0.5) < 1.0e-03 || fabs(bc[0] - 0.5) < 1.0e-03 || fabs(bc[1] + 0.5) < 1.0e-03 || fabs(bc[1] - 0.5) < 1.0e-03 || fabs(bc[2] + 0.5) < 1.0e-03 || fabs(bc[2] - 0.5) < 1.0e-03)
	{
	  tri_geo.boundaryMark() = 1;
	}
    }  
  
  std::cout << "boundary triangle are done..." << std::endl;
  
  const int& n_tet = mesh.n_geometry(3);

  for(int i = 0;i < n_tet;++ i)
    {
      GeometryBM& tet_geo = mesh.geometry(3, i);
    
      tet_geo.boundaryMark() = 0;

      const int& n_bnd = tet_geo.n_boundary();
      for(int j = 0;j < n_bnd;++ j)
	{
	  if(mesh.geometry(2, tet_geo.boundary()[j]).boundaryMark() == 1)
	    {
	      tet_geo.boundaryMark() = 1;
	      break;
	    }
	}
    }

  std::cout << "boundary tetrahedron are done..." << std::endl;  
  
}

void DFT::getBoundaryIndex()
{
  boundaryDGEleIndex.clear();
  boundaryDOFIndex.clear();
  internalDOFIndex.clear();

  const int& n_dgEle = fem_space->n_DGElement();

  for(int i = 0;i < n_dgEle;++ i)
    {
      if(fem_space->dgElement(i).boundaryMark() == 0)continue;
      DGElement<double, DIM>& the_dgEle = fem_space->dgElement(i);
      const int& dgEle_idx = the_dgEle.index();

      boundaryDGEleIndex.push_back(dgEle_idx);
    }

  std::set<int> boundaryDOFSet;
  std::set<int> internalDOFSet;

  const int& n_ele = fem_space->n_element();
  for(int i = 0;i < n_ele;++ i)
    {
      Element<double, DIM>& the_ele = fem_space->element(i);
      const int& ele_idx = the_ele.index();

      std::vector<int>& ele_dof = the_ele.dof();
      const int& n_ele_dof = ele_dof.size();

      for(int j = 0;j < n_ele_dof;++ j)
	{
	  if(fem_space->dofBoundaryMark(ele_dof[j]) == 0)
	    {
	      internalDOFSet.insert(ele_dof[j]);
	    }
	  else
	    {
	      boundaryDOFSet.insert(ele_dof[j]);
	    }
	}
    }

  int bndDof_idx = 0;
  int intDof_idx = 0;
  boundaryDOFIndex.resize(boundaryDOFSet.size());
  internalDOFIndex.resize(internalDOFSet.size());
  std::set<int>::iterator it;
  for(it = boundaryDOFSet.begin();it != boundaryDOFSet.end();++ it)
    {
      boundaryDOFIndex[bndDof_idx ++] = *it;
    }

  for(it = internalDOFSet.begin();it != internalDOFSet.end();++ it)
    {
      internalDOFIndex[intDof_idx ++] = *it;
    }

}
