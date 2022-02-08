// Solve 1D Heat equation

// dphi/dt = k(d^2 phi/dx^2)
//   k is a coefficient of diffusivity

//   dt = time step size    // Delta t
//   dx = space btw 2 nodes //Delta x

//   -----------------------------------
//
     // Central Finite Difference
     // d^2 phi / dx^2 \approx (phi[i+1] - 2phi[i] + phi[i-1])/dx^2

     //   Explicit Euler

     //   dphi/dt = RHS => (phi^(n+1) - phi^(n))/dt = RHS^n
     //                                             phi^(n+1) = phi^(n) + dt * RHS^n

     //Boundary Con. => bc[0] bc[10] = 0
     // phi[0] = 0, phi[nx-1] = 0


     //rod model                     - - - - - - - - - - -  
     //                              0 C   Heat Here    0C

#include <iostream>  //cout
#include <cmath>      // math lib.
#include <fstream>   // malloc
#include <stdlib.h>  // sleep
#include <unistd.h>

using namespace std;

void initialize(double *var, double *var_new, int nx) {
  // initialize variables
  for(int i =0; i <= nx-1; i++){
    var[i] = 0.;
    var_new[i] = 0.;
  }
  var[(nx-1)/2] = 1.;
}

void visualize(double *var, int nx){
  for(int i =0; i <= nx-1; i++){
    cout << var[i] << " ";
  }
  cout << "\n";
}

void simulation(double *var, double *var_new, int nx, double c, double dx, double dt){
  double RHS;
  //simualtion
  
  for(int i = 1; i <= nx-2; i++){   //do not put value to B.C., 0 and nx-1
    RHS = c*((var[i+1] - 2.*var[i] + var[i-1])/(dx*dx));
    var_new[i] = var[i]+dt*RHS; // explicit Euler
  }
}

void update(double *var, double *var_new, int nx){
  for(int i = 0; i <= nx-1; i++){
    var[i] = var_new[i];
  }
  
}

int main(){

  int nx = 11;
  double k  = 1.;
  double dx = 1.;
  double dt = 0.01;

  //-------------------------

  double *phi;
  phi = (double *) malloc (nx * sizeof(double));
  double *phi_new;
  phi_new = (double *) malloc (nx * sizeof(double));

  //-------------------------

  initialize(phi, phi_new, nx);
  visualize(phi,nx);

  //begin time loop

  for(int n = 1; n <= 10^6; n++){
    simulation(phi, phi_new, nx, k, dx, dt);
    update(phi, phi_new, nx);
    visualize(phi_new,nx);
    // sleep(1);
  }
  
}
