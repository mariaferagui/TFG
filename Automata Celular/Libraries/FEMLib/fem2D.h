#include <string>
#include "../../constants.hpp"
#include "../ToolsLib/toolslib.h"
using namespace constants;


double *get_nutrient_vec (double a[], double ALPHA, float COEF_DIFF,  
double element_area[], int element_node[], int element_num, double f[],
int H[], int ib, float lambda, int nnodes, int node_boundary[], int node_num,
double node_xy[], int nx, int ny, int pivot[], int quad_num, int T[], double time,
double time_step_size, double u_old[], double wq[], double xq[], double yq[]);

void set_coordinates( int nx, int ny, int node_num, double xl, double xr, double yb,
double yt, double node_xy[], int nnodes, int element_num, int element_node[], double wq[], double xq[],
  double yq[], double element_area[]);

void adjust_backward_euler ( int node_num, double node_xy[], int nnodes,
  int element_num, int element_node[], int quad_num, double wq[],
  double xq[], double yq[], double element_area[], int ib, double time,
  double time_step_size, double u_old[], double a[], double f[] );
void adjust_boundary ( int node_num, double node_xy[], int node_boundary[],
  int ib, double time, double a[], double f[], int nx, int ny );
void area_set ( int node_num, double node_xy[], int nnodes,
  int element_num, int element_node[], double element_area[] );
void assemble ( int node_num, double node_xy[], int nnodes,
  int element_num, int element_node[], int quad_num,
  double wq[], double xq[], double yq[], double element_area[],
  int ib, double time, double a[], double f[], double u_old[],
  int T[], int H[], float lambda, double ALPHA, float COEF_DIFF);
int bandwidth ( int nnodes, int element_num, int element_node[],
  int node_num );
int dgb_fa ( int n, int ml, int mu, double a[], int pivot[] );
void dgb_print_some ( int m, int n, int ml, int mu, double a[], int ilo,
  int jlo, int ihi, int jhi, std::string title );
double *dgb_sl ( int n, int ml, int mu, double a[], int pivot[],
  double b[], int job );
void element_write ( int nnodes, int element_num, int element_node[],
  std::string triangulation_txt_file_name );
double *initial_nutrients ( int node_num, double node_xy[], int nx, int ny);
// void filename_inc ( std::string *file_name );
void grid_t6 ( int nx, int ny, int nnodes, int element_num, int element_node[] );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
void i4vec_print_some ( int n, int a[], int max_print, std::string title );
int *node_boundary_set ( int nx, int ny, int node_num );
void nodes_plot ( std::string file_name, int node_num, double node_xy[],
  bool node_label );
void nodes_write ( int node_num, double node_xy[], std::string output_filename );
void qbf ( double x, double y, int element, int inode, double node_xy[],
  int element_node[], int element_num, int nnodes,
  int node_num, double *bb, double *bx, double *by );
void quad_a ( double node_xy[], int element_node[],
  int element_num, int node_num, int nnodes, double wq[], double xq[],
  double yq[] );
void quad_e ( double node_xy[], int element_node[],
  int element, int element_num, int nnodes, int node_num, int nqe,
  double wqe[], double xqe[], double yqe[] );
double r8_huge ( );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
int r8_nint ( double x );
void r8vec_print_some ( int n, double a[], int i_lo, int i_hi, std::string title );
double rhs ( double x, double y, double time, double N_old, 
int T, int H, float lambda, double ALPHA );
int s_len_trim ( std::string s );

void timestamp ( );
void xy_set ( int nx, int ny, int node_num, double xl, double xr, double yb,
  double yt, double node_xy[] );

void solution_write ( int node_num, double u[], std::string u_file_name ) ;
void solution_write ( int node_num, float u[], std::string u_file_name ) ;
