#include "DFT.h"
#include <string>

int main(int argc, char* argv[])
{
  int refined_time = atoi(argv[1]);
  std::string fileName = argv[2];
  std::cout<<"reading parameter from terminal is done..."<<std::endl;
  //////////////////////////////////////////////////////////////
  // declare the class DFT
  DFT myDFT;
  std::cout<<"Declare class DFT is done..."<<std::endl;
  //////////////////////////////////////////////////////////////
  // setup mesh
  myDFT.readMesh(fileName, refined_time);
  std::cout<<"reading mesh is done..."<<std::endl;
  //////////////////////////////////////////////////////////////
  // initialize template element information
  myDFT.initialize();
  std::cout<< "initializing template elements is done..."<<std::endl;
  //////////////////////////////////////////////////////////////
  // build the finite element space.
  myDFT.buildspace();
  
  
  //////////////////////////////////////////////////////////////
  // prepare has been done. It's time to assemble matrix and rhs 
  // to get solution of the problem.
  // assemble the matrix for poisson equation of Hartree potential.
  myDFT.buildHartreeMatrix();

  //////////////////////////////////////////////////////////////
  // initialize the electric density to start self-consistent
  // iteration.
  std::cout<<"initialize the electric density...\n";
  myDFT.initializeRho();
  std::cout<<"Density finish initializing!\n";
  //////////////////////////////////////////////////////////////
  // get the right-hand side vector of the poisson equations via
  // electric density.
  std::cout<<"get the rhs vector of the poisson equation of Hartree potential...\n";
  myDFT.getRHS_Hartree();
  std::cout<<"Apply boundary condition...\n";
  myDFT.addBoundaryCondition_Hartree();
  std::cout<<"Begin to solve the poisson equation for Hartree...\n";
  myDFT.solveHartree();
  //myDFT.printMesh(); 
}
