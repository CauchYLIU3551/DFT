#ifndef _Matrix_h_
#define _Matrix_h_

#include <AFEPack/BilinearOperator.h>

#define DIM 3
///////////////////////////////////////////////////////
// This class is used to build the stiffmatrix for KS
class MatrixA: public StiffMatrix<DIM,double>
{
public:
  MatrixA(FEMSpace<double,DIM> &sp):
    StiffMatrix<DIM,double>(sp){};
  virtual ~MatrixA();
public:
  virtual void
  getElementMatrix(const Element<double,2>& ele0,
		   const Element<double,2>& ele1,
		   const ActiveElementPairIterator<DIM>::State state);
};


#endif
