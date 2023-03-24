
# include <cmath>
# include <ctime>
# include <cstdlib>
# include <cstring>
# include <fstream>
# include <iostream>
# include <iomanip>
# include "constants.hpp"
# include "Libraries/GrowthLib/growthlib.h"
# include "Libraries/FEMLib/fem2D.h"
# include "Libraries/ToolsLib/toolslib.h"

using namespace constants;

bool node_label;
double eh1, el2, xr = 1.0, yb = 0.0, yt = 1.0, xl = 0.0;
double element_area[ELEMENT_NUM], node_xy[2*NODE_NUM], wq[QUAD_NUM], xq[QUAD_NUM*ELEMENT_NUM], yq[QUAD_NUM*ELEMENT_NUM];
double *a_N, *a_M, *f_N, *f_M, *N, *M, *N_old, *M_old;
float COEF_DIFF = 10;
float *DIV_mat;
int ib, ierr_N, ierr_M, job, node, node_show, triangle_show;
int element_node[NNODES*ELEMENT_NUM];
int *pivot_N, *pivot_M, *T, *D, *H;
int *node_boundary=(int *) malloc(sizeof(int)*NODE_NUM);
std::string folder = "Results/Growth/";
std::string time_filename = folder + "rectangle_time.txt";
std::string node_filename = folder + "rectangle_nodes.txt";
std::string triangulation_filename = folder + "rectangle_elements.txt";
std::string N_filename = folder + "N/0000.txt";
std::string M_filename = folder + "M/0000.txt";
std::string T_filename = folder + "T/0000.txt";
std::string D_filename = folder + "D/0000.txt";
std::string H_filename = folder + "H/0000.txt";
std::string DIV_filename = folder + "DIV/0000.txt";

int main ( void )
{
  double time;
  double time_final;
  double time_init;
  int time_step;
  int time_step_num;
  double time_step_size;
  std::ofstream time_unit;
  
  timestamp ( );

  // Set and save coordinates and useful matrices
  set_coordinates(NX, NY, NODE_NUM, xl, xr, yb, yt, node_xy, NNODES, 
                  ELEMENT_NUM, element_node, wq, xq, yq, element_area);
  node_boundary = node_boundary_set (NX, NY, NODE_NUM);
  ib = bandwidth (NNODES, ELEMENT_NUM, element_node, NODE_NUM);
  nodes_write (NODE_NUM, node_xy, node_filename );
  element_write (NNODES, ELEMENT_NUM, element_node, triangulation_filename);

  // Memory allocation
  N = new double[NODE_NUM];
  M = new double[NODE_NUM];
  T = new int[NODE_NUM];
  D = new int[NODE_NUM];
  H = new int[NODE_NUM];
  a_N = new double[(3*ib+1)*NODE_NUM];
  a_M = new double[(3*ib+1)*NODE_NUM];
  f_N = new double[NODE_NUM];
  f_M = new double[NODE_NUM];
  N_old = new double[NODE_NUM];
  M_old = new double[NODE_NUM];
  pivot_N = new int[NODE_NUM];
  pivot_M = new int[NODE_NUM];
  DIV_mat = new float[NODE_NUM];

  // Times
  time_init = 0.0;
  time_final = 0.5;
  time_step_num = GENERATION_IT;
  time_step_size = ( time_final - time_init ) / ( double ) ( time_step_num );
  time = time_init;

  
  // Create nutrient vectors
  N = initial_nutrients ( NODE_NUM, node_xy, NX, NY);
  M = initial_nutrients ( NODE_NUM, node_xy, NX, NY);

  time_unit.open ( time_filename.c_str ( ) );
  time_unit << "  " << std::setw(14) << time << "\n";

  // Create T H y D vectors
  create_vec(NODE_NUM, DIV_mat, DIV);
  create_vec(NODE_NUM, T, 0);
  create_vec(NODE_NUM, H, 1);
  create_vec(NODE_NUM, D, 0);

  // Place a T cell and save matrices
  T[int((2*NX -1)*(2*NY - 1)/2)] = 1;
  H[int((2*NX -1)*(2*NY - 1)/2)] = 0;
  solution_write ( NODE_NUM, N, N_filename );
  solution_write ( NODE_NUM, M, M_filename );
  save_mat(NODE_NUM, T, T_filename);
  
  // Onset of the tumor mass
  for ( time_step = 1; time_step <= time_step_num; time_step++ )
  {
    for ( node = 0; node < NODE_NUM; node++ )
    {
      N_old[node] = N[node];
      M_old[node] = M[node];
    }

    time = ( ( double ) ( time_step_num - time_step ) * time_init
           + ( double ) (                 time_step ) * time_final )
           / ( double ) ( time_step_num             );
    
    // FEM applied to nutrients
    N = get_nutrient_vec(a_N, ALPHA, COEF_DIFF, element_area, element_node, ELEMENT_NUM, 
                        f_N, H, ib, L_N, NNODES, node_boundary, NODE_NUM, node_xy, NX, NY,
                        pivot_N, QUAD_NUM, T, time, time_step_size, N_old, wq, xq, yq);
    M = get_nutrient_vec(a_M, ALPHA, COEF_DIFF, element_area, element_node, ELEMENT_NUM,
                        f_M, H, ib, L_M, NNODES, node_boundary, NODE_NUM, node_xy, NX, NY,
                        pivot_M, QUAD_NUM, T, time, time_step_size, M_old, wq, xq, yq);
    
    // Tumor grows

    grow(M, N, T, D, H, DIV_mat, 2*NX-1, 2*NY-1);

    time_unit << std::setw(14) << time << "\n";

    // Increase filename
    filename_inc ( &T_filename );
    filename_inc ( &N_filename );

    // Save each 5 time units
    if(time_step%5==0){
      save_mat(NODE_NUM, T, T_filename);
      save_mat(NODE_NUM, N, N_filename);
    }
  
    if(metastasis(T, 2*NX - 1, 2*NY - 1)){
      break;
    }

    std::cout<<"ITERATION: "<<time_step<<std::endl;
  }

  // Free memory
  delete [] N;
  delete [] M;
  delete [] T;
  delete [] H;
  delete [] D;
  delete [] a_N;
  delete [] a_M;
  delete [] f_N;
  delete [] f_M;
  delete [] N_old;
  delete [] M_old;
  delete [] pivot_N;
  delete [] pivot_M;
  delete [] node_boundary;

  time_unit.close( );
  std::cout <<"\nFinished\n";
  timestamp( );

  return 0;
}
