/* Description
 Testing function degsv for soiving Ax = b
 with 3X3

*/


#include <iostream>
#include <lapacke.h>

using namespace std;

int main(){

  const int N =3, NRHS = 1, LDA = N, LDB = N;
  int ipiv[N], info;
  // double A[LDA*N] = {
  // 		     1, 0, 0,
  //     		     0, 1, 0,
  // 		     0, 0, 1
  // };
   double A[LDA*N] = {
		     2,-4, 0,
      		     3, 1, 9,
		     0, 5, 7
  };

  double B[LDB*NRHS] = {
			1.7, 0.5, -3.2
  };

  info = LAPACKE_dgesv( LAPACK_COL_MAJOR, N, NRHS, A, LDA, ipiv, B, LDB);
  for(int i = 0; i <= N-1; i++){
    cout << B[i] << "\t";
  }
  cout <<"\n";
  
}
