#include "../include/DFT.h"

void DFT::buildHartreeMatrix()
{
  /*
  if(stiff_Hartree != NULL){
    std::cout << "deleting the current stiff_matrix" << std::endl;
    delete stiff_Hartree;
    std::cout << "deleting successfully!\n";
    stiff_Hartree = new StiffMatrix<DIM, double>(*fem_space);
  }*/

  std::cout << "build a new stiff_matrix for Hatree potential..." << std::endl;
  stiff_Hartree = new StiffMatrix<DIM, double>(*fem_space);
  stiff_Hartree->algebricAccuracy() = ACC;
  stiff_Hartree->build();
}
