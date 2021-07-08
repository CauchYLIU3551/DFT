#include "../include/DFT.h"

void DFT::solveHartree()
{
  AMGSolver solver(*stiff_Hartree);
  solver.solve(*V_hartree, *rhs, 1.0e-08, 50);
}
