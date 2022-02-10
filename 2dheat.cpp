// Solve 2D Heat equation

// dphi/dt = k(d^2 phi/dx^2)
//   k is a coefficient of diffusivity

//   dt = time step size    // Delta t
//   dx = space btw 2 nodes //Delta x
//   dy = space btw 2 nodes // Delta y
//   -----------------------------------
//
     // Central Finite Difference
     // d^2 phi/dx^2 \approx (phi[i+1][j] - 2phi[i][j] + phi[i-1][j])/dx^2
     // d^2 phi/dy^2 \approx (phi[i][j+1] - 2phi[i][j] + phi[i][j-1])/dx^2

     //   Explicit Euler

     //   dphi/dt = RHS => (phi^(n+1) - phi^(n))/dt = RHS^n
     //                                             phi^(n+1) = phi^(n) + dt * RHS^n

     //Boundary Con. => bc[0] bc[10] = 0
     // phi[0] = 0, phi[nx-1] = 0


     //rod model                     - - - - - - - - - - -  
     //                              0 C   Heat Here    0C

#include <iostream>  //cout
#include <cmath>     // math lib.
#include <fstream>   // malloc
#include <stdlib.h>  // sleep
#include <unistd.h>
#include <string>

using namespace std;

void initialize(double **var, double **var_new, int nx, int ny) {
  // initialize variables
  for(int i =0; i <= nx-1; i++){
  for(int j =0; j <= ny-1; j++){
  var[i][j] = 0.;
  var_new[i][j] = 0.;
  }
  }
  var[(nx-1)/2][(ny-1)/2] = 1.;
  var[(nx-1)/2-1][(ny-1)/2] = 1.;
  var[(nx-1)/2][(ny-1)/2-1] = 1.;
  var[(nx-1)/2+1][(ny-1)/2] = 1.;
  var[(nx-1)/2][(ny-1)/2+1] = 1.;
  var[(nx-1)/2-1][(ny-1)/2+1] = 1.;
  var[(nx-1)/2+1][(ny-1)/2-1] = 1.;
  var[(nx-1)/2+1][(ny-1)/2+1] = 1.;
  var[(nx-1)/2-1][(ny-1)/2-1] = 1.;
}

void visualize(double **var, int nx, int ny){
  for(int i =0; i <= nx-1; i++){
    for(int j =0; j <= ny-1; j++){
      cout << var[i][j] << "\t ";
    }
    cout << "\n";
  }
  cout << "\n";
}

void simulation(double **var, double **var_new, int nx, int ny, double c, double dx, double dy, double dt){
  double RHS;
  //simualtion
  for(int i = 1; i <= nx-2; i++){   //do not put value to B.C., 0 and nx-1
    for(int j = 1; j <= ny-2; j++){
      RHS = c*((var[i+1][j] - 2.*var[i][j] + var[i-1][j])/(dx*dx)) +
	c*((var[i][j+1] - 2.*var[i][j] + var[i][j-1])/(dy*dy));

      var_new[i][j] = var[i][j] + dt*RHS; // explicit Euler
    }
  }
}

void update(double **var, double **var_new, int nx, int ny){
    for(int i = 0; i <= nx-1; i++){
      for(int j = 0; j <= ny-1; j++){
	var[i][j] = var_new[i][j];
      }
    }
}

void paraview(string fileName, double **var, int nx, int ny, double dx, double dy){
  ofstream myfile;
  myfile.open(fileName);
  //----------------------------------------------//

  //Paraview header
  myfile << "# vtk DataFile Version 2.0\n";
  myfile << "FlowField\n";
  myfile << "ASCII\n";

  //Grid
  myfile << "DATASET STRUCTURED_GRID\n";
  myfile << "DIMENSIONS " << nx << " " << 1 << " " << ny << "\n";
  myfile << "POINTS " << nx*1*ny << " float\n";
  for(int j= 0; j <= ny-1; j++){
    for(int i = 0; i <= nx-1;i++){
      myfile << dx*i << " " << dy*j << " 0\n";
    }
  }

  //Data
  myfile << "\n";
  myfile << "POINT_DATA";
  myfile << " " << nx*ny << "\n";

  myfile << "\n";
  myfile << "SCALARS PHI float 1\n";
  myfile << "LOOKUP_TABLE default\n";

  for(int j=0; j <= ny-1; j++){
    for(int i =0; i <= nx-1; i++){
      myfile << var[i][j] << "\n";
    }
  }
  myfile.close();
}
   
      
  

int main(){

  int nx = 11;
  int ny = 11;
  double k  = 1.;
  double dx = 1.;
  double dy = 1.;
  double dt = 0.001;
  string fileName; 

  //-------------------------

  double **phi;
  phi = (double **) malloc (nx * sizeof(double));
  for (int i=0; i <= nx; i++){
    phi[i] = (double *) malloc (ny *sizeof(double));
  }

  
  double **phi_new;
  phi_new = (double **) malloc (nx * sizeof(double));
  for (int i=0; i <= nx; i++){
    phi_new[i] = (double *) malloc (ny *sizeof(double));
  }

  //-------------------------

  initialize(phi, phi_new, nx, ny);
  // visualize(phi, nx, ny);
  // paraview(phi, nx, ny, dx, dy);

   //begin time loop
   for(int n = 1; n <= 100000; n++){
       simulation(phi, phi_new, nx, ny, k, dx, dy, dt);
       update(phi, phi_new, nx, ny);
       // visualize(phi, nx, ny);
       //sleep(1);
       cout << "n = " << n << "\n";
       if(n%100 == 0){
	 fileName = "var_" + to_string(n) + ".vtk";
         paraview(fileName, phi, nx, ny, dx, dy);
       }
   }
   //  fileName = "var.vtk";
  
}

