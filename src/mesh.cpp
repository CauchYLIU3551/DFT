#include "../include/DFT.h"

void DFT::readMesh(const std::string& fileName,
			  const int& refine_times)
{

  h_tree = new HGeometryTree<DIM>();
  h_tree->readMesh(fileName);

  irregular_mesh = new IrregularMesh<DIM>();
  irregular_mesh->reinit(*h_tree);
  irregular_mesh->globalRefine(refine_times);
  irregular_mesh->semiregularize();
  irregular_mesh->regularize(false);

}
