


ofstream myfile;
myfile.open("phi_paraview.vtk");

 //Paraview Header
myfile << "# vtk DataFile Version 2.0\n";
myfile << "FlowField\n";
myfile << "ASCII\n";


 //Grid
myfile << "DATASET STRUCTURED_GRID\n";
myfile << "DIMENSIONS " << nx << " " << 1 << " " << nx << "\n";
myfile <<  "POINTS " << nx*1*ny << " float\n";
for(j = 0; j <=ny-1; j++){
  for(i = 0; i <= nx-1; i++){
      myfile << i << " " << j << " 0\n";
  }
 }
 //Data
myfile << "\n";
myfile << "POINT_DATA";
myfile << " " << nx*ny << "\n";

myfile << "\n";
myfile << "SCALARS PHI float 1\n";
myfile << "LOOKUP_TABLE default\n";
for(j = 0; j <=ny-1; j++){
  for(i = 0; i <= nx-1; i++){
    myfile << phi[i][j] << "\n";
  }
 }

myfile.close();
