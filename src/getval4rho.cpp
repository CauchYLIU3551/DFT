#include "../include/DFT.h"

void DFT::getval4rho()
{
    for (u_int i = 0;i < rho->size();i ++)
    {
      (*rho)[i] = 2 * (*phi)[i] * (*phi)[i];
      // std::cout<<(*rho)[i]<<" ";
    }
}
