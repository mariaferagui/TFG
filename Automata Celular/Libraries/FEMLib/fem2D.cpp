# include <iostream>
# include <iomanip>
# include <fstream>
# include <cmath>
# include <cstring>
#include <string>
# include "fem2D.h"


double *get_nutrient_vec (double a[], double ALPHA, float COEF_DIFF,  
double element_area[], int element_node[], int element_num, double f[],
int H[], int ib, float lambda, int nnodes, int node_boundary[], int node_num,
double node_xy[], int nx, int ny, int pivot[], int quad_num, int T[], double time,
double time_step_size, double u_old[], double wq[], double xq[], double yq[])
  {
    double *x;
    int ierr;
    x = new double[node_num];
    // Assemble A and F matrices
    assemble (node_num, node_xy, nnodes, element_num, element_node, quad_num, 
    wq, xq, yq, element_area, ib,  time, a, f, u_old, T, H, lambda, ALPHA, COEF_DIFF);

    // Modify A and F to get dM/dt and dN/dt
    adjust_backward_euler(node_num, node_xy, nnodes, element_num, element_node,  quad_num, wq,
    xq, yq,  element_area, ib,  time, time_step_size, u_old, a,  f );

    // Modify A and F for given initial conditions
    adjust_boundary ( node_num, node_xy, node_boundary, ib, time, a, f, nx, ny );

    // Solve the system
    ierr = dgb_fa ( node_num, ib, ib, a, pivot);
     if ( ierr != 0 )
    {
      std::cout << "\n";
      std::cout << "FEM2D_HEAT - Fatal error!\n";
      std::cout << "  DGB_FA returned the error condition IERR = " << ierr << ".\n";
      std::cout << "\n";
      std::cout << "  The linear system was not factored, and the\n";
      std::cout << "  algorithm cannot proceed.\n";
      exit ( 1 );
    }
    x = dgb_sl ( node_num, ib, ib, a, pivot, f, 0);
   return x;

  }

void set_coordinates( int nx, int ny, int node_num, double xl, double xr, double yb,
  double yt, double node_xy[], int nnodes, int element_num, int element_node[], double wq[], double xq[],
  double yq[], double element_area[])
{
  int *node_boundary=(int *) malloc(sizeof(int)*node_num);
  xy_set(nx, ny, node_num, xl, xr, yb, yt, node_xy);
  grid_t6 (nx, ny, nnodes, element_num, element_node);
  quad_a (node_xy, element_node, element_num, node_num, nnodes, wq, xq, yq);
  area_set (node_num, node_xy, nnodes, element_num, element_node, element_area);

};

void adjust_backward_euler ( int node_num, double node_xy[], int nnodes,
  int element_num, int element_node[], int quad_num, double wq[],
  double xq[], double yq[], double element_area[], int ib, double time,
  double time_step_size, double u_old[], double a[], double f[] )

//****************************************************************************80
//
//  Purpose:
//
//    ADJUST_BACKWARD_EULER adjusts the system for the backward Euler term.
//
//  Discussion:
//
//    The input linear system
//
//      A * U = F
//
//    is appropriate for the equation
//
//      -Uxx - Uyy = RHS
//
//    We need to modify the matrix A and the right hand side F to
//    account for the approximation of the time derivative in
//
//      Ut - Uxx - Uyy = RHS
//
//    by the backward Euler approximation:
//
//      Ut approximately equal to ( U - Uold ) / dT
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double NODE_XY[2*NODE_NUM], the coordinates of nodes.
//
//    Input, int NNODES, the number of nodes used to form one element.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input, int ELEMENT_NODE[NNODES*ELEMENT_NUM];
//    ELEMENT_NODE(I,J) is the global index of local node I in element J.
//
//    Input, int QUAD_NUM, the number of quadrature points used in assembly.
//
//    Input, double WQ[QUAD_NUM], quadrature weights.
//
//    Input, double XQ[QUAD_NUM*ELEMENT_NUM],
//    YQ[QUAD_NUM*ELEMENT_NUM], the
//    coordinates of the quadrature points in each element.
//
//    Input, double ELEMENT_AREA[ELEMENT_NUM], the area of elements.
//
//    Input, int IB, the half-bandwidth of the matrix.
//
//    Input, double TIME, the current time.
//
//    Input, double TIME_STEP_SIZE, the size of the time step.
//
//    Input, double U_OLD[NODE_NUM], the finite element
//    coefficients for the solution at the previous time.
//
//    Input/output, double A[(3*IB+1)*NODE_NUM], the NODE_NUM
//    by NODE_NUM coefficient matrix, stored in a compressed format.
//
//    Input/output, double F[NODE_NUM], the right hand side.
//
{
  int basis;
  double bi;
  double bj;
  double dbidx;
  double dbidy;
  double dbjdx;
  double dbjdy;
  int element;
  int j;
  int node;
  int quad;
  int test;
  double w;
  double x;
  double y;

  for ( element = 0; element < element_num; element++ )
  {
    for ( quad = 0; quad < quad_num; quad++ )
    {
      x = xq[quad+element*quad_num];
      y = yq[quad+element*quad_num];
      w = element_area[element] * wq[quad];

      for ( test = 0; test < nnodes; test++ )
      {
        node = element_node[test+element*nnodes];

        qbf ( x, y, element, test, node_xy, element_node,
          element_num, nnodes, node_num, &bi, &dbidx, &dbidy );
//
//  Carry the U_OLD term to the right hand side.
//
        f[node] = f[node] + w * bi * u_old[node] / time_step_size;
//
//  Modify the diagonal entries of A.
//
        for ( basis = 0; basis < nnodes; basis++ )
        {
          j = element_node[basis+element*nnodes];

          qbf ( x, y, element, basis, node_xy, element_node,
            element_num, nnodes, node_num, &bj, &dbjdx, &dbjdy );

          a[node-j+2*ib+j*(3*ib+1)] = a[node-j+2*ib+j*(3*ib+1)]
            + w * bi * bj / time_step_size;
        }
      }
    }
  }

  return;
}
//****************************************************************************80

void adjust_boundary ( int node_num, double node_xy[], int node_boundary[],
  int ib, double time, double a[], double f[], int nx, int ny )

//****************************************************************************80
//
//  Purpose:
//
//    ADJUST_BOUNDARY modifies the linear system for boundary conditions.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double NODE_XY[2*NODE_NUM], the coordinates of nodes.
//
//    Input, int NODE_BOUNDARY[NODE_NUM], is
//    0, if a node is an interior node;
//    1, if a node is a Dirichlet boundary node.
//
//    Input, int IB, the half-bandwidth of the matrix.
//
//    Input, double TIME, the current time.
//
//    Input/output, double A[(3*IB+1)*NODE_NUM], the NODE_NUM by NODE_NUM
//    coefficient matrix, stored in a compressed format.
//    On output, A has been adjusted for boundary conditions.
//
//    Input/output, double F[NODE_NUM], the right hand side.
//    On output, F has been adjusted for boundary conditions.
//
{
  int j;
  int jhi;
  int jlo;
  int node;
  double *u_exact;
//
//  Get the exact solution at every node.
//
  u_exact = new double[node_num];

  u_exact = initial_nutrients( node_num, node_xy, nx, ny);

  for ( node = 0; node < node_num; node++ )
  {
    if ( node_boundary[node] != 0 )
    {
      jlo = i4_max ( node - ib, 0 );
      jhi = i4_min ( node + ib, node_num - 1 );

      for ( j = jlo; j <= jhi; j++ )
      {
        a[node-j+2*ib+j*(3*ib+1)] = 0.0;
      }
      a[node-node+2*ib+node*(3*ib+1)] = 1.0;

      f[node] = u_exact[node];
    }
  }

  delete [] u_exact;

  return;
}
//****************************************************************************80

void area_set ( int node_num, double node_xy[], int nnodes,
  int element_num, int element_node[], double element_area[] )

//****************************************************************************80
//
//  Purpose:
//
//    AREA_SET sets the area of each element.
//
//  Discussion:
//
//    The areas of the elements are needed in order to adjust
//    the integral estimates produced by the quadrature formulas.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 April 2004
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double NODE_XY[2*NODE_NUM], the
//    coordinates of the nodes.
//
//    Input, int NNODES, the number of local nodes per element.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input, int ELEMENT_NODE[NNODES*ELEMENT_NUM];
//    ELEMENT_NODE(I,J) is the global index of local node I in element J.
//
//    Output, double ELEMENT_AREA[ELEMENT_NUM], the area of elements.
//
{
  int element;
  int i1;
  int i2;
  int i3;
  double x1;
  double x2;
  double x3;
  double y1;
  double y2;
  double y3;

  for ( element = 0; element < element_num; element++ )
  {
    i1 = element_node[0+element*nnodes];
    x1 = node_xy[0+i1*2];
    y1 = node_xy[1+i1*2];

    i2 = element_node[1+element*nnodes];
    x2 = node_xy[0+i2*2];
    y2 = node_xy[1+i2*2];

    i3 = element_node[2+element*nnodes];
    x3 = node_xy[0+i3*2];
    y3 = node_xy[1+i3*2];

    element_area[element] = 0.5 * fabs
      ( y1 * ( x2 - x3 )
      + y2 * ( x3 - x1 )
      + y3 * ( x1 - x2 ) );
  }

  return;
}
//****************************************************************************80*

void assemble ( int node_num, double node_xy[], int nnodes, int element_num,
  int element_node[], int quad_num, double wq[], double xq[], double yq[],
  double element_area[], int ib, double time, double a[], double f[],
  double u_old[], int T[], int H[], float lambda, double ALPHA, float COEF_DIFF)

//****************************************************************************80*
//
//  Purpose:
//
//    ASSEMBLE assembles the matrix and right-hand side using piecewise quadratics.
//
//  Discussion:
//
//    The matrix is known to be banded.  A special matrix storage format
//    is used to reduce the space required.  Details of this format are
//    discussed in the routine DGB_FA.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 April 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double NODE_XY[2*NODE_NUM], the X and Y coordinates of nodes.
//
//    Input, int NNODES, the number of nodes used to form one element.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input, int ELEMENT_NODE[NNODES*ELEMENT_NUM]; ELEMENT_NODE(I,J) is the global
//    index of local node I in element J.
//
//    Input, int QUAD_NUM, the number of quadrature points used in assembly.
//
//    Input, double WQ[QUAD_NUM], quadrature weights.
//
//    Input, double XQ[QUAD_NUM*ELEMENT_NUM], YQ[QUAD_NUM*ELEMENT_NUM], the X and Y
//    coordinates of the quadrature points in each element.
//
//    Input, double ELEMENT_AREA[ELEMENT_NUM], the area of each element.
//
//    Input, int IB, the half-bandwidth of the matrix.
//
//    Input, double TIME, the current time.
//
//    Output, double A[(3*IB+1)*NODE_NUM], the NODE_NUM by NODE_NUM
//    coefficient matrix, stored in a compressed format.
//
//    Output, double F[NODE_NUM], the right hand side.
//
//  Local parameters:
//
//    Local, double BI, DBIDX, DBIDY, the value of some basis function
//    and its first derivatives at a quadrature point.
//
//    Local, double BJ, DBJDX, DBJDY, the value of another basis
//    function and its first derivatives at a quadrature point.
//
{
  double aij;
  int basis;
  double bi;
  double bj;
  double dbidx;
  double dbidy;
  double dbjdx;
  double dbjdy;
  int element;
  int i;
  int j;
  int node;
  int quad;
  int test;
  double w;
  double x;
  double y;
//
//  Initialize the arrays to zero.
//
  for ( i = 0; i < node_num; i++ )
  {
    f[i] = 0.0;
  }

  for ( j = 0; j < node_num; j++ )
  {
    for ( i = 0; i < 3*ib + 1; i++ )
    {
      a[i+j*(3*ib+1)] = 0.0;
    }
  }
//
//  The actual values of A and F are determined by summing up
//  contributions from all the elements.
//
  for ( element = 0; element < element_num; element++ )
  {
    for ( quad = 0; quad < quad_num; quad++ )
    {
      x = xq[quad+element*quad_num];
      y = yq[quad+element*quad_num];
      w = element_area[element] * wq[quad];

      for ( test = 0; test < nnodes; test++ )
      {
        node = element_node[test+element*nnodes];

        qbf ( x, y, element, test, node_xy, element_node,
          element_num, nnodes, node_num, &bi, &dbidx, &dbidy );

        f[node] = f[node] + w * rhs ( x, y, time, u_old[node], T[node], H[node], lambda, ALPHA) * bi;
        
//  We are about to compute a contribution associated with the
//  I-th test function and the J-th basis function, and add this
//  to the entry A(I,J).
//
//  Because of the compressed storage of the matrix, the element
//  will actually be stored in A(I-J+2*IB+1,J).
//
//  An extra complication: we are storing the array as a vector.
//
//  Therefore, we ACTUALLY store the entry in A[I-J+2*IB+1-1 + J * (3*IB+1)];
//
        for ( basis = 0; basis < nnodes; basis++ )
        {
          j = element_node[basis+element*nnodes];

          qbf ( x, y, element, basis, node_xy, element_node,
            element_num, nnodes, node_num, &bj, &dbjdx, &dbjdy );

          aij = COEF_DIFF*(dbidx * dbjdx + dbidy * dbjdy);

          a[node-j+2*ib+j*(3*ib+1)] = a[node-j+2*ib+j*(3*ib+1)] + w * aij;
        }
      }
    }
  }

  return;
}
//****************************************************************************80

int bandwidth ( int nnodes, int element_num, int element_node[], int node_num )

//****************************************************************************80
//
//  Purpose:
//
//    BANDWIDTH determines the bandwidth of the coefficient matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 April 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Parameters:
//
//    Input, int NNODES, the number of local nodes per element.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input, int ELEMENT_NODE[NNODES*ELEMENT_NUM]; ELEMENT_NODE(I,J) is the global
//    index of local node I in element J.
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Output, int BANDWIDTH, the half bandwidth of the matrix.
//
{
  int element;
  int i;
  int iln;
  int j;
  int jln;
  int nhba;

  nhba = 0;

  for ( element = 0; element < element_num; element++ )
  {
    for ( iln = 0; iln < nnodes; iln++ )
    {
      i = element_node[iln+element*nnodes];
      for ( jln = 0; jln < nnodes; jln++ )
      {
        j = element_node[jln+element*nnodes];
        nhba = i4_max ( nhba, j - i );
      }
    }
  }

  return nhba;
}
//****************************************************************************80
int dgb_fa ( int n, int ml, int mu, double a[], int pivot[] )

//****************************************************************************80
//
//  Purpose:
//
//    DGB_FA performs a LINPACK-style PLU factorization of an DGB matrix.
//
//  Discussion:
//
//    The DGB storage format is used for an M by N banded matrix, with lower bandwidth ML
//    and upper bandwidth MU.  Storage includes room for ML extra superdiagonals,
//    which may be required to store nonzero entries generated during Gaussian
//    elimination.
//
//    The original M by N matrix is "collapsed" downward, so that diagonals
//    become rows of the storage array, while columns are preserved.  The
//    collapsed array is logically 2*ML+MU+1 by N.
//
//    The two dimensional array can be further reduced to a one dimensional
//    array, stored by columns.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 February 2004
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative, and no greater than N-1.
//
//    Input/output, double A[(2*ML+MU+1)*N], the matrix in band storage.
//    On output, A has been overwritten by the LU factors.
//
//    Output, int PIVOT[N], the pivot vector.
//
//    Output, int SGB_FA, singularity flag.
//    0, no singularity detected.
//    nonzero, the factorization failed on the INFO-th step.
//
{
  int col = 2 * ml + mu + 1;
  int i;
  int i0;
  int j;
  int j0;
  int j1;
  int ju;
  int jz;
  int k;
  int l;
  int lm;
  int m;
  int mm;
  double t;

  m = ml + mu + 1;
//
//  Zero out the initial fill-in columns.
//
  j0 = mu + 2;
  j1 = i4_min ( n, m ) - 1;

  for ( jz = j0; jz <= j1; jz++ )
  {
    i0 = m + 1 - jz;
    for ( i = i0; i <= ml; i++ )
    {
      a[i-1+(jz-1)*col] = 0.0;
    }
  }

  jz = j1;
  ju = 0;

  for ( k = 1; k <= n-1; k++ )
  {
//
//  Zero out the next fill-in column.
//
    jz = jz + 1;
    if ( jz <= n )
    {
      for ( i = 1; i <= ml; i++ )
      {
        a[i-1+(jz-1)*col] = 0.0;
      }
    }
//
//  Find L = pivot index.
//
    lm = i4_min ( ml, n-k );
    l = m;

    for ( j = m+1; j <= m + lm; j++ )
    {
      if ( fabs ( a[l-1+(k-1)*col] ) < fabs ( a[j-1+(k-1)*col] ) )
      {
        l = j;
      }
    }

    pivot[k-1] = l + k - m;
//
//  Zero pivot implies this column already triangularized.
//
    if ( a[l-1+(k-1)*col] == 0.0 )
    {
      std::cout << "\n";
      std::cout << "DGB_FA - Fatal error!\n";
      std::cout << "  Zero pivot on step " << k << "\n";
      return k;
    }
//
//  Interchange if necessary.
//
    t                = a[l-1+(k-1)*col];
    a[l-1+(k-1)*col] = a[m-1+(k-1)*col];
    a[m-1+(k-1)*col] = t;
//
//  Compute multipliers.
//
    for ( i = m+1; i <= m+lm; i++ )
    {
      a[i-1+(k-1)*col] = - a[i-1+(k-1)*col] / a[m-1+(k-1)*col];
    }
//
//  Row elimination with column indexing.
//
    ju = i4_max ( ju, mu + pivot[k-1] );
    ju = i4_min ( ju, n );
    mm = m;

    for ( j = k+1; j <= ju; j++ )
    {
      l = l - 1;
      mm = mm - 1;

      if ( l != mm )
      {
        t                 = a[l-1+(j-1)*col];
        a[l-1+(j-1)*col]  = a[mm-1+(j-1)*col];
        a[mm-1+(j-1)*col] = t;
      }
      for ( i = 1; i <= lm; i++ )
      {
        a[mm+i-1+(j-1)*col] = a[mm+i-1+(j-1)*col]
          + a[mm-1+(j-1)*col] * a[m+i-1+(k-1)*col];
      }
    }
  }

  pivot[n-1] = n;

  if ( a[m-1+(n-1)*col] == 0.0 )
  {
    std::cout << "\n";
    std::cout << "DGB_FA - Fatal error!\n";
    std::cout << "  Zero pivot on step " << n << "\n";
    return n;
  }

  return 0;
}
//****************************************************************************80

void dgb_print_some ( int m, int n, int ml, int mu, double a[], int ilo,
  int jlo, int ihi, int jhi, std::string title )

//****************************************************************************80
//
//  Purpose:
//
//    DGB_PRINT_SOME prints some of a DGB matrix.
//
//  Discussion:
//
//    The DGB storage format is used for an M by N banded matrix, with lower bandwidth ML
//    and upper bandwidth MU.  Storage includes room for ML extra superdiagonals,
//    which may be required to store nonzero entries generated during Gaussian
//    elimination.
//
//    The original M by N matrix is "collapsed" downward, so that diagonals
//    become rows of the storage array, while columns are preserved.  The
//    collapsed array is logically 2*ML+MU+1 by N.
//
//    The two dimensional array can be further reduced to a one dimensional
//    array, stored by columns.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows of the matrix.
//    M must be positive.
//
//    Input, int N, the number of columns of the matrix.
//    N must be positive.
//
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative, and no greater than min(M,N)-1..
//
//    Input, double A[(2*ML+MU+1)*N], the SGB matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, std::string TITLE, a title to print.
{
# define INCX 5

  int col = 2 * ml + mu + 1;
  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  if ( 0 < s_len_trim ( title ) )
  {
    std::cout << "\n";
    std::cout << title << "\n";
  }
//
//  Print the columns of the matrix, in strips of 5.
//
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    std::cout << "\n";
    std::cout << "  Col: ";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      std::cout << std::setw(7) << j << "       ";
    }
    std::cout << "\n";
    std::cout << "  Row\n";
    std::cout << "  ---\n";
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );
    i2lo = i4_max ( i2lo, j2lo - mu );

    i2hi = i4_min ( ihi, m );
    i2hi = i4_min ( i2hi, j2hi + ml );

    for ( i = i2lo; i <= i2hi; i++ )
    {
//
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
      std::cout << std::setw(6) << i << "  ";
      for ( j = j2lo; j <= j2hi; j++ )
      {
        if ( ml < i-j || mu < j-i )
        {
          std::cout << "            ";
        }
        else
        {
          std::cout << std::setw(10) << a[i-j+ml+mu+(j-1)*col] << "  ";
        }
      }
      std::cout << "\n";
    }
  }

  std::cout << "\n";

  return;
# undef INCX
}
//****************************************************************************80

double *dgb_sl ( int n, int ml, int mu, double a[], int pivot[],
  double b[], int job )

//****************************************************************************80
//
//  Purpose:
//
//    DGB_SL solves a system factored by DGB_FA.
//
//  Discussion:
//
//    The DGB storage format is used for an M by N banded matrix, with lower bandwidth ML
//    and upper bandwidth MU.  Storage includes room for ML extra superdiagonals,
//    which may be required to store nonzero entries generated during Gaussian
//    elimination.
//
//    The original M by N matrix is "collapsed" downward, so that diagonals
//    become rows of the storage array, while columns are preserved.  The
//    collapsed array is logically 2*ML+MU+1 by N.
//
//    The two dimensional array can be further reduced to a one dimensional
//    array, stored by columns.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 February 2004
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative, and no greater than N-1.
//
//    Input, double A[(2*ML+MU+1)*N], the LU factors from DGB_FA.
//
//    Input, int PIVOT[N], the pivot vector from DGB_FA.
//
//    Input, double B[N], the right hand side vector.
//
//    Input, int JOB.
//    0, solve A * x = b.
//    nonzero, solve A' * x = b.
//
//    Output, double DGB_SL[N], the solution.
//
{
  int col = 2 * ml + mu + 1;
  int i;
  int k;
  int l;
  int la;
  int lb;
  int lm;
  int m;
  double t;
  double *x;

  x = new double[n];

  for ( i = 0; i < n; i++ )
  {
    x[i] = b[i];
  }
//
  m = mu + ml + 1;
//
//  Solve A * x = b.
//
  if ( job == 0 )
  {
//
//  Solve L * Y = B.
//
    if ( 1 <= ml )
    {
      for ( k = 1; k <= n-1; k++ )
      {
        lm = i4_min ( ml, n-k );
        l = pivot[k-1];

        if ( l != k )
        {
          t      = x[l-1];
          x[l-1] = x[k-1];
          x[k-1] = t;
        }
        for ( i = 1; i <= lm; i++ )
        {
          x[k+i-1] = x[k+i-1] + x[k-1] * a[m+i-1+(k-1)*col];
        }
      }
    }
//
//  Solve U * X = Y.
//
    for ( k = n; 1 <= k; k-- )
    {
      x[k-1] = x[k-1] / a[m-1+(k-1)*col];
      lm = i4_min ( k, m ) - 1;
      la = m - lm;
      lb = k - lm;
      for ( i = 0; i <= lm-1; i++ )
      {
        x[lb+i-1] = x[lb+i-1] - x[k-1] * a[la+i-1+(k-1)*col];
      }
    }
  }
//
//  Solve A' * X = B.
//
  else
  {
//
//  Solve U' * Y = B.
//
    for ( k = 1; k <= n; k++ )
    {
      lm = i4_min ( k, m ) - 1;
      la = m - lm;
      lb = k - lm;
      for ( i = 0; i <= lm-1; i++ )
      {
        x[k-1] = x[k-1] - x[lb+i-1] * a[la+i-1+(k-1)*col];
      }
      x[k-1] = x[k-1] / a[m-1+(k-1)*col];
    }
//
//  Solve L' * X = Y.
//
    if ( 1 <= ml )
    {
      for ( k = n-1; 1 <= k; k-- )
      {
        lm = i4_min ( ml, n-k );
        for ( i = 1; i <= lm; i++ )
        {
          x[k-1] = x[k-1] + x[k+i-1] * a[m+i-1+(k-1)*col];
        }
        l = pivot[k-1];

        if ( l != k )
        {
          t      = x[l-1];
          x[l-1] = x[k-1];
          x[k-1] = t;
        }
      }
    }
  }
  for(int node = 0;node<NODE_NUM;node++)changeNegativeValue(x[node]);
  
  return x;
}
//****************************************************************************80*

void element_write ( int nnodes, int element_num, int element_node[],
 std::string output_filename )

//****************************************************************************80*
//
//  Purpose:
//
//    ELEMENT_WRITE writes the elements to a file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NNODES, the number of nodes used to form one element.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input, int ELEMENT_NODE[NNODES*ELEMENT_NUM]; ELEMENT_NODE(I,J) is the global
//    index of local node I in element J.
//
//    Input, std::string OUTPUT_FILENAME, the name of the file
//    in which the data should be stored.
//
{
  int element;
  int i;
  std::ofstream output;

  output.open ( output_filename.c_str ( ) );

  if ( !output )
  {
    std::cout << "\n";
    std::cout << "ELEMENT_WRITE - Warning!\n";
    std::cout << "  Could not write the node file.\n";
    return;
  }

  for ( element = 0; element < element_num; element++ )
  {
    for ( i = 0; i < nnodes; i++ )
    {
      output << std::setw(8)  << element_node[i+element*nnodes] << "  ";
    }
    output << "\n";
  }

  output.close ( );

  return;
}
//****************************************************************************80*

double *initial_nutrients ( int node_num, double node_xy[], int nx, int ny )

//****************************************************************************80
//
//  Purpose:
//
//    initial_nutrients creates the initial nutrient matrix.
//
//  Discussion:
//
//    It is assumed that the user knows the exact solution and its
//    derivatives.  This, of course, is NOT true for a real computation.
//    But for this code, we are interested in studying the convergence
//    behavior of the approximations, and so we really need to assume
//    we know the correct solution.
//
//    As a convenience, this single routine is used for several purposes:
//
//    * it supplies the initial value function H(X,Y,T);
//    * it supplies the boundary value function G(X,Y,T);
//    * it is used by the COMPARE routine to make a node-wise comparison
//      of the exact and approximate solutions.
//    * it is used by the ERRORS routine to estimate the integrals of
//      the L2 and H1 errors of approximation.
//
//    DUDX and DUDY are only needed for the ERRORS calculation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 April 2004
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes at which
//    a value is desired.
//
//    Input, double NODE_XY(2,NODE_NUM), the coordinates of
//    the points where a value is desired.
//
//    Input, double TIME, the current time.
//
//    Output, double U[NODE_NUM], the exact solution.
//
//    Output, double DUDX[NODE_NUM], DUDY[NODE_NUM],
//    the X and Y derivatives of the exact solution.
//
{
  int node; node = 0;
  int var;
  double *u;
  u = new double[node_num];

  node = 0;
  for ( int j = 1; j <= 2 * ny - 1; j++ )
  {
    for ( int i = 1; i <= 2 * nx - 1; i++ )
    {
      if ( (j == 1 ) || (j == 2 * ny - 1))
      {
        u[node] = STARTING_NUTRIENTS;
      }
      else
      {
        u[node] = 0;
      }
      node = node + 1;
    }
  }

  return u;
}

//****************************************************************************80

// void filename_inc ( std::string *filename )

// //****************************************************************************80
// //
// //  Purpose:
// //
// //    FILENAME_INC increments a partially numeric file name.
// //
// //  Discussion:
// //
// //    It is assumed that the digits in the name, whether scattered or
// //    connected, represent a number that is to be increased by 1 on
// //    each call.  If this number is all 9's on input, the output number
// //    is all 0's.  Non-numeric letters of the name are unaffected.
// //
// //    If the name is empty, then the routine stops.
// //
// //    If the name contains no digits, the empty std::string is returned.
// //
// //  Example:
// //
// //      Input            Output
// //      -----            ------
// //      "a7to11.txt"     "a7to12.txt"  (typical case.  Last digit incremented)
// //      "a7to99.txt"     "a8to00.txt"  (last digit incremented, with carry.)
// //      "a9to99.txt"     "a0to00.txt"  (wrap around)
// //      "cat.txt"        " "           (no digits to increment)
// //      " "              STOP!         (error)
// //
// //  Licensing:
// //
// //    This code is distributed under the GNU LGPL license. 
// //
// //  Modified:
// //
// //    22 November 2011
// //
// //  Author:
// //
// //    John Burkardt
// //
// //  Parameters:
// //
// //    Input/output, std::string *FILENAME, the filename to be incremented.
// //
// {
//   char c;
//   int change;
//   int i;
//   int lens;
//   std::string numstr;
//   int num;

//   lens = (*filename).length ( );

//   if ( lens <= 0 )
//   {
//     std::cout << "\n";
//     std::cout << "FILENAME_INC - Fatal error!\n";
//     std::cout << "  The input std::string is empty.\n";
//     exit ( 1 );
//   }

//   change = 0;
//   for ( i = lens - 1; 0 <= i; i-- )
//   {
//     c = (*filename)[i];
//     if ( '0' <= c && c <= '9' )
//     {
//       change = change + 1;

//       if ( c == '9' )
//       {
//         c = '0';
//         (*filename)[i] = c;
//       }
//       else
//       {
//         c = c + 1;
//         (*filename)[i] = c;
//         return;
//       }
//     }
//   }
// //
// //  No digits were found.  Return blank.
// //
//   if ( change == 0 )
//   {
//     for ( i = lens - 1; 0 <= i; i-- )
//     {
//       (*filename)[i] = ' ';
//     }
//   }

//   return;
// }
//****************************************************************************80

void grid_t6 ( int nx, int ny, int nnodes, int element_num, int element_node[] )

//****************************************************************************80
//
//  Purpose:
//
//    GRID_T6 produces a grid of pairs of 6 node triangles.
//
//  Example:
//
//    Input:
//
//      NX = 4, NY = 3
//
//    Output:
//
//      ELEMENT_NODE =
//         1,  3, 15,  2,  9,  8;
//        17, 15,  3, 16,  9, 10;
//         3,  5, 17,  4, 11, 10;
//        19, 17,  5, 18, 11, 12;
//         5,  7, 19,  6, 13, 12;
//        21, 19,  7, 20, 13, 14;
//        15, 17, 29, 16, 23, 22;
//        31, 29, 17, 30, 23, 24;
//        17, 19, 31, 18, 25, 24;
//        33, 31, 19, 32, 25, 26;
//        19, 21, 33, 20, 27, 26;
//        35, 33, 21, 34, 27, 28.
//
//  Diagram:
//
//   29-30-31-32-33-34-35
//    |\ 8  |\10  |\12  |
//    | \   | \   | \   |
//   22 23 24 25 26 27 28
//    |   \ |   \ |   \ |
//    |  7 \|  9 \| 11 \|
//   15-16-17-18-19-20-21
//    |\ 2  |\ 4  |\ 6  |
//    | \   | \   | \   |
//    8  9 10 11 12 13 14
//    |   \ |   \ |   \ |
//    |  1 \|  3 \|  5 \|
//    1--2--3--4--5--6--7
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NX, NY, controls the number of elements along the
//    X and Y directions.  The number of elements will be
//    2 * ( NX - 1 ) * ( NY - 1 ).
//
//    Input, int NNODES, the number of local nodes per element.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Output, int ELEMENT_NODE[NNODES*ELEMENT_NUM];
//    ELEMENT_NODE(I,J) is the index of the I-th node of the J-th element.
//
{
  int c;
  int e;
  int element;
  int i;
  int j;
  int n;
  int ne;
  int nw;
  int s;
  int se;
  int sw;
  int w;

  element = 0;

  for ( j = 1; j <= ny - 1; j++ )
  {
    for ( i = 1; i <= nx - 1; i++ )
    {
      sw = ( j - 1 ) * 2 * ( 2 * nx - 1 ) + 2 * i - 2;
      w  = sw + 1;
      nw = sw + 2;

      s  = sw + 2 * nx - 1;
      c  = s + 1;
      n  = s + 2;

      se = s  + 2 * nx - 1;
      e  = se + 1;
      ne = se + 2;

      element_node[0+element*nnodes] = sw;
      element_node[1+element*nnodes] = se;
      element_node[2+element*nnodes] = nw;
      element_node[3+element*nnodes] = s;
      element_node[4+element*nnodes] = c;
      element_node[5+element*nnodes] = w;
      element = element + 1;

      element_node[0+element*nnodes] = ne;
      element_node[1+element*nnodes] = nw;
      element_node[2+element*nnodes] = se;
      element_node[3+element*nnodes] = n;
      element_node[4+element*nnodes] = c;
      element_node[5+element*nnodes] = e;
      element = element + 1;
    }
  }

  return;
}
//****************************************************************************80

int i4_max ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MAX returns the maximum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, are two ints to be compared.
//
//    Output, int I4_MAX, the larger of I1 and I2.
//
{
  if ( i2 < i1 )
  {
    return i1;
  }
  else
  {
    return i2;
  }

}
//****************************************************************************80

int i4_min ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MIN returns the smaller of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, two ints to be compared.
//
//    Output, int I4_MIN, the smaller of I1 and I2.
//
{
  if ( i1 < i2 )
  {
    return i1;
  }
  else
  {
    return i2;
  }

}
//****************************************************************************80

void i4vec_print_some ( int n, int a[], int max_print, std::string title )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_PRINT_SOME prints "some" of an I4VEC.
//
//  Discussion:
//
//    The user specifies MAX_PRINT, the maximum number of lines to print.
//
//    If N, the size of the vector, is no more than MAX_PRINT, then
//    the entire vector is printed, one entry per line.
//
//    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
//    followed by a line of periods suggesting an omission,
//    and the last entry.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 November 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries of the vector.
//
//    Input, int A[N], the vector to be printed.
//
//    Input, int MAX_PRINT, the maximum number of lines to print.
//
//    Input, std::string TITLE, an optional title.
//
{
  int i;

  if ( max_print <= 0 )
  {
    return;
  }

  if ( n <= 0 )
  {
    return;
  }

  if ( 0 < s_len_trim ( title ) )
  {
    std::cout << "\n";
    std::cout << title << "\n";
    std::cout << "\n";
  }

  if ( n <= max_print )
  {
    for ( i = 0; i < n; i++ )
    {
      std::cout << std::setw(6)  << i + 1 << "  "
           << std::setw(10) << a[i] << "\n";
    }
  }
  else if ( 3 <= max_print )
  {
    for ( i = 0; i < max_print-2; i++ )
    {
      std::cout << std::setw(6)  << i + 1 << "  "
           << std::setw(10) << a[i]  << "\n";
    }
    std::cout << "......  ..............\n";
    i = n - 1;
    std::cout << std::setw(6)  << i + 1 << "  "
         << std::setw(10) << a[i]  << "\n";
  }
  else
  {
    for ( i = 0; i < max_print-1; i++ )
    {
      std::cout << std::setw(6)  << i + 1 << "  "
           << std::setw(10) << a[i]  << "\n";
    }
    i = max_print - 1;
    std::cout << std::setw(6)  << i + 1 << "  "
         << std::setw(10) << a[i]  << "...more entries...\n";
  }

  return;
}
//****************************************************************************80

int *node_boundary_set ( int nx, int ny, int node_num )

//****************************************************************************80
//
//  Purpose:
//
//    NODE_BOUNDARY_SET assigns an unknown value index at each node.
//
//  Discussion:
//
//    Every node is assigned a value which indicates whether it is
//    an interior node, or a boundary node.
//
//  Example:
//
//    On a simple 5 by 5 grid, where the nodes are numbered starting
//    at the lower left, and increasing in X first, we would have the
//    following values of NODE_BOUNDARY:
//
//       1  1  1  1  1
//       1  0  0  0  1
//       1  0  0  0  1
//       1  0  0  0  1
//       1  1  1  1  1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 April 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Parameters:
//
//    Input, int NX, NY, the number of elements in the X and Y directions.
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Output, int NODE_BOUNDARY[NODE_NUM], is
//    0, if a node is an interior node;
//    1, if a node is a Dirichlet boundary node.
//
{
  int i;
  int j;
  int node;
  int *node_boundary;

  node_boundary = new int[node_num];

  node = 0;

  for ( j = 1; j <= 2 * ny - 1; j++ )
  {
    for ( i = 1; i <= 2 * nx - 1; i++ )
    {
      if ( j == 1 ||
           j == 2 * ny - 1)
           //i == 1 ||
           //i == 2 * nx - 1 )
      {
        node_boundary[node] = 1;
      }
      else
      {
        node_boundary[node] = 0;
      }

      node = node + 1;
    }
  }
  return node_boundary;
}
//****************************************************************************80

void nodes_plot ( std::string file_name, int node_num, double node_xy[],
  bool node_label )

//****************************************************************************80
//
//  Purpose:
//
//    NODES_PLOT plots a pointset.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 September 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, std::string FILE_NAME, the name of the file to create.
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double NODE_XY[2*NODE_NUM], the nodes.
//
//    Input, bool NODE_LABEL, is TRUE if the nodes are to be labeled.
//
{
  int circle_size;
  int delta;
  std::ofstream file_unit;
  int node;
  double x_max;
  double x_min;
  int x_ps;
  int x_ps_max = 576;
  int x_ps_max_clip = 594;
  int x_ps_min = 36;
  int x_ps_min_clip = 18;
  double x_scale;
  double y_max;
  double y_min;
  int y_ps;
  int y_ps_max = 666;
  int y_ps_max_clip = 684;
  int y_ps_min = 126;
  int y_ps_min_clip = 108;
  double y_scale;
//
//  We need to do some figuring here, so that we can determine
//  the range of the data, and hence the height and width
//  of the piece of paper.
//
  x_max = -r8_huge ( );
  for ( node = 0; node < node_num; node++ )
  {
     if ( x_max < node_xy[0+node*2] )
     {
       x_max = node_xy[0+node*2];
     }
  }
  x_min = r8_huge ( );
  for ( node = 0; node < node_num; node++ )
  {
     if ( node_xy[0+node*2] < x_min )
     {
       x_min = node_xy[0+node*2];
     }
  }
  x_scale = x_max - x_min;

  x_max = x_max + 0.05 * x_scale;
  x_min = x_min - 0.05 * x_scale;
  x_scale = x_max - x_min;

  y_max = -r8_huge ( );
  for ( node = 0; node < node_num; node++ )
  {
     if ( y_max < node_xy[1+node*2] )
     {
       y_max = node_xy[1+node*2];
     }
  }
  y_min = r8_huge ( );
  for ( node = 0; node < node_num; node++ )
  {
     if ( node_xy[1+node*2] < y_min )
     {
       y_min = node_xy[1+node*2];
     }
  }
  y_scale = y_max - y_min;

  y_max = y_max + 0.05 * y_scale;
  y_min = y_min - 0.05 * y_scale;
  y_scale = y_max - y_min;

  if ( x_scale < y_scale )
  {
    delta = r8_nint ( ( double ) ( x_ps_max - x_ps_min )
      * ( y_scale - x_scale ) / ( 2.0 * y_scale ) );

    x_ps_max = x_ps_max - delta;
    x_ps_min = x_ps_min + delta;

    x_ps_max_clip = x_ps_max_clip - delta;
    x_ps_min_clip = x_ps_min_clip + delta;

    x_scale = y_scale;
  }
  else if ( y_scale < x_scale )
  {
    delta = r8_nint ( ( double ) ( y_ps_max - y_ps_min )
      * ( x_scale - y_scale ) / ( 2.0 * x_scale ) );

    y_ps_max = y_ps_max - delta;
    y_ps_min = y_ps_min + delta;

    y_ps_max_clip = y_ps_max_clip - delta;
    y_ps_min_clip = y_ps_min_clip + delta;

    y_scale = x_scale;
  }

  file_unit.open ( file_name.c_str ( ) );

  if ( !file_unit )
  {
    std::cout << "\n";
    std::cout << "POINTS_PLOT - Fatal error!\n";
    std::cout << "  Could not open the output EPS file.\n";
    exit ( 1 );
  }

  file_unit << "%!PS-Adobe-3.0 EPSF-3.0\n";
  file_unit << "%%Creator: nodes_plot.cpp\n";
  file_unit << "%%Title: " << file_name << "\n";
  file_unit << "%%Pages: 1\n";
  file_unit << "%%BoundingBox:  "
    << x_ps_min << "  "
    << y_ps_min << "  "
    << x_ps_max << "  "
    << y_ps_max << "\n";
  file_unit << "%%Document-Fonts: Times-Roman\n";
  file_unit << "%%LanguageLevel: 1\n";
  file_unit << "%%EndComments\n";
  file_unit << "%%BeginProlog\n";
  file_unit << "/inch {72 mul} def\n";
  file_unit << "%%EndProlog\n";
  file_unit << "%%Page:      1     1\n";
  file_unit << "save\n";
  file_unit << "%\n";
  file_unit << "% Set the RGB line color to very light gray.\n";
  file_unit << "%\n";
  file_unit << " 0.9000 0.9000 0.9000 setrgbcolor\n";
  file_unit << "%\n";
  file_unit << "% Draw a gray border around the page.\n";
  file_unit << "%\n";
  file_unit << "newpath\n";
  file_unit << x_ps_min << "  "
            << y_ps_min << "  moveto\n";
  file_unit << x_ps_max << "  "
            << y_ps_min << "  lineto\n";
  file_unit << x_ps_max << "  "
            << y_ps_max << "  lineto\n";
  file_unit << x_ps_min << "  "
            << y_ps_max << "  lineto\n";
  file_unit << x_ps_min << "  "
            << y_ps_min << "  lineto\n";
  file_unit << "stroke\n";
  file_unit << "%\n";
  file_unit << "% Set RGB line color to black.\n";
  file_unit << "%\n";
  file_unit << " 0.0000 0.0000 0.0000 setrgbcolor\n";
  file_unit << "%\n";
  file_unit << "%  Set the font and its size:\n";
  file_unit << "%\n";
  file_unit << "/Times-Roman findfont\n";
  file_unit << "0.50 inch scalefont\n";
  file_unit << "setfont\n";
  file_unit << "%\n";
  file_unit << "%  Print a title:\n";
  file_unit << "%\n";
  file_unit << "%  210  702 moveto\n";
  file_unit << "%(Pointset) show\n";
  file_unit << "%\n";
  file_unit << "% Define a clipping polygon\n";
  file_unit << "%\n";
  file_unit << "newpath\n";
  file_unit << x_ps_min_clip << "  "
            << y_ps_min_clip << "  moveto\n";
  file_unit << x_ps_max_clip << "  "
            << y_ps_min_clip << "  lineto\n";
  file_unit << x_ps_max_clip << "  "
            << y_ps_max_clip << "  lineto\n";
  file_unit << x_ps_min_clip << "  "
            << y_ps_max_clip << "  lineto\n";
  file_unit << x_ps_min_clip << "  "
            << y_ps_min_clip << "  lineto\n";
  file_unit << "clip newpath\n";
//
//  Draw the nodes.
//
  if ( node_num <= 200 )
  {
    circle_size = 5;
  }
  else if ( node_num <= 500 )
  {
    circle_size = 4;
  }
  else if ( node_num <= 1000 )
  {
    circle_size = 3;
  }
  else if ( node_num <= 5000 )
  {
    circle_size = 2;
  }
  else
  {
    circle_size = 1;
  }

  file_unit << "%\n";
  file_unit << "%  Draw filled dots at each node:\n";
  file_unit << "%\n";
  file_unit << "%  Set the color to blue:\n";
  file_unit << "%\n";
  file_unit << "0.000  0.150  0.750  setrgbcolor\n";
  file_unit << "%\n";

  for ( node = 0; node < node_num; node++ )
  {
    x_ps = ( int ) (
      ( ( x_max - node_xy[0+node*2]         ) * ( double ) ( x_ps_min )
      + (       + node_xy[0+node*2] - x_min ) * ( double ) ( x_ps_max ) )
      / ( x_max                     - x_min ) );

    y_ps = ( int ) (
      ( ( y_max - node_xy[1+node*2]         ) * ( double ) ( y_ps_min )
      + (         node_xy[1+node*2] - y_min ) * ( double ) ( y_ps_max ) )
      / ( y_max                     - y_min ) );

    file_unit << "newpath  "
      << x_ps << "  "
      << y_ps << "  "
      << circle_size << " 0 360 arc closepath fill\n";
  }
//
//  Label the nodes.
//
  file_unit << "%\n";
  file_unit << "%  Label the nodes:\n";
  file_unit << "%\n";
  file_unit << "%  Set the color to darker blue:\n";
  file_unit << "%\n";
  file_unit << "0.000  0.250  0.850  setrgbcolor\n";
  file_unit << "/Times-Roman findfont\n";
  file_unit << "0.20 inch scalefont\n";
  file_unit << "setfont\n";

  file_unit << "%\n";

  for ( node = 0; node < node_num; node++ )
  {
    x_ps = ( int ) (
      ( ( x_max - node_xy[0+node*2]         ) * ( double ) ( x_ps_min )
      + (       + node_xy[0+node*2] - x_min ) * ( double ) ( x_ps_max ) )
      / ( x_max                     - x_min ) );

    y_ps = ( int ) (
      ( ( y_max - node_xy[1+node*2]         ) * ( double ) ( y_ps_min )
      + (         node_xy[1+node*2] - y_min ) * ( double ) ( y_ps_max ) )
      / ( y_max                     - y_min ) );

    file_unit << "newpath  "
      << x_ps     << "  "
      << y_ps + 5 << "  moveto ("
      << node     << ") show\n";
  }

  file_unit << "%\n";
  file_unit << "restore showpage\n";
  file_unit << "%\n";
  file_unit << "% End of page\n";
  file_unit << "%\n";
  file_unit << "%%Trailer\n";
  file_unit << "%%EOF\n";

  file_unit.close ( );

  return;
}
//****************************************************************************80*

void nodes_write ( int node_num, double node_xy[], std::string output_filename )

//****************************************************************************80*
//
//  Purpose:
//
//    NODES_WRITE writes the nodes to a file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double NODE_XY[2*NODE_NUM], the X and Y coordinates of nodes.
//
//    Input, std::string OUTPUT_FILENAME, the name of the file
//    in which the data should be stored.
//
{
  int node;
  std::ofstream output;
  double x;
  double y;

  output.open ( output_filename.c_str ( ) );

  if ( !output )
  {
    std::cout << "\n";
    std::cout << "NODES_WRITE - Warning!\n";
    std::cout << "  Could not write the node file.\n";
    return;
  }

  for ( node = 0; node < node_num; node++ )
  {
    x = node_xy[0+node*2];
    y = node_xy[1+node*2];

    output << std::setw(8)  << x << "  "
           << std::setw(8)  << y << "\n";
  }

  output.close ( );

  return;
}
//****************************************************************************80

void qbf ( double x, double y, int element, int inode, double node_xy[],
  int element_node[], int element_num, int nnodes,
  int node_num, double *b, double *dbdx, double *dbdy )

//****************************************************************************80
//
//  Purpose:
//
//    QBF evaluates the quadratic basis functions.
//
//  Discussion:
//
//    This routine assumes that the "midpoint" nodes are, in fact,
//    exactly the average of the two extreme nodes.  This is NOT true
//    for a general quadratic triangular element.
//
//    Assuming this property of the midpoint nodes makes it easy to
//    determine the values of (R,S) in the reference element that
//    correspond to (X,Y) in the physical element.
//
//    Once we know the (R,S) coordinates, it's easy to evaluate the
//    basis functions and derivatives.
//
//  The physical element T6:
//
//    In this picture, we don't mean to suggest that the bottom of
//    the physical triangle is horizontal.  However, we do assume that
//    each of the sides is a straight line, and that the intermediate
//    points are exactly halfway on each side.
//
//    |
//    |
//    |        3
//    |       . .
//    |      .   .
//    Y     6     5
//    |    .       .
//    |   .         .
//    |  1-----4-----2
//    |
//    +--------X-------->
//
//  Reference element T6:
//
//    In this picture of the reference element, we really do assume
//    that one side is vertical, one horizontal, of length 1.
//
//    |
//    |
//    1  3
//    |  |.
//    |  | .
//    S  6  5
//    |  |   .
//    |  |    .
//    0  1--4--2
//    |
//    +--0--R--1-------->
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the (global) coordinates of the point
//    at which the basis function is to be evaluated.
//
//    Input, int ELEMENT, the index of the element which contains the point.
//
//    Input, int INODE, the local index, between 0 and 5, that
//    specifies which basis function is to be evaluated.
//
//    Input, double NODE_XY[2*NODE_NUM], the nodes.
//
//    Input, int ELEMENT_NODE[NNODES*ELEMENT_NUM];
//    ELEMENT_NODE(I,J) is the global index of local node I in element J.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input, int NNODES, the number of nodes used to form one element.
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Output, double *B, *DBDX, *DBDY, the value of the basis function
//    and its X and Y derivatives at (X,Y).
//
{
  double dbdr;
  double dbds;
  double det;
  double drdx;
  double drdy;
  double dsdx;
  double dsdy;
  int i;
  double r;
  double s;
  double xn[6];
  double yn[6];

  for ( i = 0; i < 6; i++ )
  {
    xn[i] = node_xy[0+element_node[i+element*nnodes]*2];
    yn[i] = node_xy[1+element_node[i+element*nnodes]*2];
  }
//
//  Determine the (R,S) coordinates corresponding to (X,Y).
//
//  What is happening here is that we are solving the linear system:
//
//    ( X2-X1  X3-X1 ) * ( R ) = ( X - X1 )
//    ( Y2-Y1  Y3-Y1 )   ( S )   ( Y - Y1 )
//
//  by computing the inverse of the coefficient matrix and multiplying
//  it by the right hand side to get R and S.
//
//  The values of dRdX, dRdY, dSdX and dSdY are easily from the formulas
//  for R and S.
//
  det =   ( xn[1] - xn[0] ) * ( yn[2] - yn[0] )
        - ( xn[2] - xn[0] ) * ( yn[1] - yn[0] );

  r = ( ( yn[2] - yn[0] ) * ( x     - xn[0] )
      + ( xn[0] - xn[2] ) * ( y     - yn[0] ) ) / det;

  drdx = ( yn[2] - yn[0] ) / det;
  drdy = ( xn[0] - xn[2] ) / det;

  s = ( ( yn[0] - yn[1] ) * ( x     - xn[0] )
      + ( xn[1] - xn[0] ) * ( y     - yn[0] ) ) / det;

  dsdx = ( yn[0] - yn[1] ) / det;
  dsdy = ( xn[1] - xn[0] ) / det;
//
//  The basis functions can now be evaluated in terms of the
//  reference coordinates R and S.  It's also easy to determine
//  the values of the derivatives with respect to R and S.
//
  if ( inode == 0 )
  {
    *b   =   2.0 *     ( 1.0 - r - s ) * ( 0.5 - r - s );
    dbdr = - 3.0 + 4.0 * r + 4.0 * s;
    dbds = - 3.0 + 4.0 * r + 4.0 * s;
  }
  else if ( inode == 1 )
  {
    *b   =   2.0 * r * ( r - 0.5 );
    dbdr = - 1.0 + 4.0 * r;
    dbds =   0.0;
  }
  else if ( inode == 2 )
  {
    *b   =   2.0 * s * ( s - 0.5 );
    dbdr =   0.0;
    dbds = - 1.0               + 4.0 * s;
  }
  else if ( inode == 3 )
  {
    *b   =   4.0 * r * ( 1.0 - r - s );
    dbdr =   4.0 - 8.0 * r - 4.0 * s;
    dbds =           - 4.0 * r;
  }
  else if ( inode == 4 )
  {
    *b   =   4.0 * r * s;
    dbdr =                           4.0 * s;
    dbds =             4.0 * r;
  }
  else if ( inode == 5 )
  {
    *b   =   4.0 * s * ( 1.0 - r - s );
    dbdr =                         - 4.0 * s;
    dbds =   4.0 - 4.0 * r - 8.0 * s;
  }
  else
  {
    std::cout << "\n";
    std::cout << "QBF - Fatal error!\n";
    std::cout << "  Request for local basis function INODE = " << inode << "\n";
    exit ( 1 );
  }
//
//  We need to convert the derivative information from (R(X,Y),S(X,Y))
//  to (X,Y) using the chain rule.
//
  *dbdx = dbdr * drdx + dbds * dsdx;
  *dbdy = dbdr * drdy + dbds * dsdy;

  return;
}
//****************************************************************************80

void quad_a ( double node_xy[], int element_node[],
  int element_num, int node_num, int nnodes, double wq[], double xq[],
  double yq[] )

//****************************************************************************80
//
//  Purpose:
//
//    QUAD_A sets the quadrature rule for assembly.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 April 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Parameters:
//
//    Input, double NODE_XY[2*NODE_NUM], the nodes.
//
//    Input, int ELEMENT_NODE[NNODES*ELEMENT_NUM];
//    ELEMENT_NODE(I,J) is the global index of local node I in element J.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int NNODES, the number of nodes used to form one element.
//
//    Output, double WQ[3], quadrature weights.
//
//    Output, double XQ[3*ELEMENT_NUM], YQ[3*ELEMENT_NUM], the
//    coordinates of the quadrature points in each element.
//
{
  int element;
  int ip1;
  int ip2;
  int ip3;
  double x1;
  double x2;
  double x3;
  double y1;
  double y2;
  double y3;

  wq[0] = 1.0 / 3.0;
  wq[1] = wq[0];
  wq[2] = wq[0];

  for ( element = 0; element < element_num; element++ )
  {
    ip1 = element_node[0+element*nnodes];
    ip2 = element_node[1+element*nnodes];
    ip3 = element_node[2+element*nnodes];

    x1 = node_xy[0+ip1*2];
    x2 = node_xy[0+ip2*2];
    x3 = node_xy[0+ip3*2];

    y1 = node_xy[1+ip1*2];
    y2 = node_xy[1+ip2*2];
    y3 = node_xy[1+ip3*2];

    xq[0+element*3] = 0.5 * ( x1 + x2 );
    xq[1+element*3] = 0.5 * ( x2 + x3 );
    xq[2+element*3] = 0.5 * ( x1 + x3 );

    yq[0+element*3] = 0.5 * ( y1 + y2 );
    yq[1+element*3] = 0.5 * ( y2 + y3 );
    yq[2+element*3] = 0.5 * ( y1 + y3 );
  }
  return;
}
//****************************************************************************80*

void quad_e ( double node_xy[], int element_node[],
  int element, int element_num, int nnodes, int node_num, int nqe,
  double wqe[], double xqe[], double yqe[] )

//****************************************************************************80*
//
//  Purpose:
//
//    QUAD_E sets a quadrature rule for the error calculation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 April 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Parameters:
//
//    Input, double NODE_XY[2*NODE_NUM], the X and Y coordinates of nodes.
//
//    Input, int ELEMENT_NODE[NNODES*ELEMENT_NUM]; ELEMENT_NODE(I,J) is the global
//    index of local node I in element J.
//
//    Input, int ELEMENT, the index of the element for which the quadrature
//    points are to be computed.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input, int NNODES, the number of nodes used to form one element.
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int NQE, the number of points in the quadrature rule.
//    This is actually fixed at 13.
//
//    Output, double WQE[NQE], the quadrature weights.
//
//    Output, double XQE[NQE], YQE[NQE], the X and Y coordinates
//    of the quadrature points.
//
{
  int i;
  int ii;
  int iii;
  int ip1;
  int ip2;
  int ip3;
  double x1;
  double x2;
  double x3;
  double y1;
  double y2;
  double y3;
  double z1;
  double z2;
  double z3;
  double z4;
  double z5;
  double z6;
  double z7;

  for ( i = 0; i < 3; i++ )
  {
    wqe[i] = 0.175615257433204;
    ii = i + 3;
    wqe[ii] = 0.053347235608839;
    ii = i + 6;
    iii = ii + 3;
    wqe[ii] = 0.077113760890257;
    wqe[iii] = wqe[ii];
  }

  wqe[12] = -0.14957004446767;

  z1 = 0.479308067841923;
  z2 = 0.260345966079038;
  z3 = 0.869739794195568;
  z4 = 0.065130102902216;
  z5 = 0.638444188569809;
  z6 = 0.312865496004875;
  z7 = 0.048690315425316;

  ip1 = element_node[0+element*nnodes];
  ip2 = element_node[1+element*nnodes];
  ip3 = element_node[2+element*nnodes];

  x1 = node_xy[0+ip1*2];
  x2 = node_xy[0+ip2*2];
  x3 = node_xy[0+ip3*2];

  y1 = node_xy[1+ip1*2];
  y2 = node_xy[1+ip2*2];
  y3 = node_xy[1+ip3*2];

  xqe[ 0] = z1 * x1 + z2 * x2 + z2 * x3;
  yqe[ 0] = z1 * y1 + z2 * y2 + z2 * y3;
  xqe[ 1] = z2 * x1 + z1 * x2 + z2 * x3;
  yqe[ 1] = z2 * y1 + z1 * y2 + z2 * y3;
  xqe[ 2] = z2 * x1 + z2 * x2 + z1 * x3;
  yqe[ 2] = z2 * y1 + z2 * y2 + z1 * y3;
  xqe[ 3] = z3 * x1 + z4 * x2 + z4 * x3;
  yqe[ 3] = z3 * y1 + z4 * y2 + z4 * y3;
  xqe[ 4] = z4 * x1 + z3 * x2 + z4 * x3;
  yqe[ 4] = z4 * y1 + z3 * y2 + z4 * y3;
  xqe[ 5] = z4 * x1 + z4 * x2 + z3 * x3;
  yqe[ 5] = z4 * y1 + z4 * y2 + z3 * y3;
  xqe[ 6] = z5 * x1 + z6 * x2 + z7 * x3;
  yqe[ 6] = z5 * y1 + z6 * y2 + z7 * y3;
  xqe[ 7] = z5 * x1 + z7 * x2 + z6 * x3;
  yqe[ 7] = z5 * y1 + z7 * y2 + z6 * y3;
  xqe[ 8] = z6 * x1 + z5 * x2 + z7 * x3;
  yqe[ 8] = z6 * y1 + z5 * y2 + z7 * y3;
  xqe[ 9] = z6 * x1 + z7 * x2 + z5 * x3;
  yqe[ 9] = z6 * y1 + z7 * y2 + z5 * y3;
  xqe[10] = z7 * x1 + z5 * x2 + z6 * x3;
  yqe[10] = z7 * y1 + z5 * y2 + z6 * y3;
  xqe[11] = z7 * x1 + z6 * x2 + z5 * x3;
  yqe[11] = z7 * y1 + z6 * y2 + z5 * y3;
  xqe[12] = ( x1 + x2 + x3 ) / 3.0;
  yqe[12] = ( y1 + y2 + y3 ) / 3.0;

  return;
}
//****************************************************************************80

double r8_huge ( void )

//****************************************************************************80
//
//  Purpose:
//
//    R8_HUGE returns a "huge" R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_HUGE, a "huge" value.
//
{
  return ( double ) HUGE_VAL;
}
//****************************************************************************80

double r8_max ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MAX returns the maximum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 January 2002
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MAX, the maximum of X and Y.
//
{
  if ( y < x )
  {
    return x;
  }
  else
  {
    return y;
  }
}
//****************************************************************************80

double r8_min ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MIN returns the minimum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input double X, Y, the quantities to compare.
//
//    Output, double R8_MIN, the minimum of X and Y.
//
{
  if ( y < x )
  {
    return y;
  }
  else
  {
    return x;
  }
}
//****************************************************************************80

int r8_nint ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_NINT returns the nearest integer to an R8.
//
//  Example:
//
//        X        R8_NINT
//
//      1.3         1
//      1.4         1
//      1.5         1 or 2
//      1.6         2
//      0.0         0
//     -0.7        -1
//     -1.1        -1
//     -1.6        -2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the value.
//
//    Output, int R8_NINT, the nearest integer to X.
//
{
  int s;

  if ( x < 0.0 )
  {
    s = -1;
  }
  else
  {
    s = 1;
  }

  return ( s * ( int ) ( fabs ( x ) + 0.5 ) );
}
//****************************************************************************80

void r8vec_print_some ( int n, double a[], int i_lo, int i_hi, std::string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_PRINT_SOME prints "some" of an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8 values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries of the vector.
//
//    Input, double A[N], the vector to be printed.
//
//    Input, integer I_LO, I_HI, the first and last indices to print.
//    The routine expects 1 <= I_LO <= I_HI <= N.
//
//    Input, std::string TITLE, an optional title.
//
{
  int i;

  if ( 0 < s_len_trim ( title ) )
  {
    std::cout << "\n";
    std::cout << title << "\n";
  }

  std::cout << "\n";
  for ( i = i4_max ( 1, i_lo ); i <= i4_min ( n, i_hi ); i++ )
  {
    std::cout << "  " << std::setw(8)  << i       << "  "
         << "  " << std::setw(14) << a[i-1]  << "\n";
  }

  return;
}
//****************************************************************************80

double rhs ( double x, double y, double time, double N_old, int T, int H,
 float lambda, double ALPHA )

//****************************************************************************80
//
//  Purpose:
//
//    RHS gives the right-hand side of the differential equation.
//
//  Discussion:
//
//    The function specified here depends on the problem being
//    solved.  This is one of the routines that a user will
//    normally want to change.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 April 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the coordinates of a point
//    in the region, at which the right hand side of the
//    differential equation is to be evaluated.
//
//    Input, double TIME, the current time.
//
//    Output, double RHS, the value of the right
//    hand side of the differential equation at (X,Y).
//
{
  double value;
  value = -pow(ALPHA,2)*N_old*(lambda*T + H);  //INTRODUCIR CEL CITOTOXICAS
  return value;
}
//****************************************************************************80

int s_len_trim ( std::string s )

//****************************************************************************80
//
//  Purpose:
//
//    S_LEN_TRIM returns the length of a std::string to the last nonblank.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 October 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, std::string S, a std::string.
//
//    Output, int S_LEN_TRIM, the length of the std::string to the last nonblank.
//    If S_LEN_TRIM is 0, then the std::string is entirely blank.
//
{
  int n;

  n = s.length ( );

  while ( 0 < n )
  {
    if ( s[n-1] != ' ' && s[n-1] != '\n' )
    {
      return n;
    }
    n = n - 1;
  }

  return n;
}
//****************************************************************************80*

void solution_write ( int node_num, double u[], std::string u_file_name )

//****************************************************************************80*
//
//  Purpose:
//
//    SOLUTION_WRITE writes the solution to a file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double U[NODE_NUM], the coefficients of the solution.
//
//    Input, std::string U_FILE_NAME, the name of the file
//    in which the data should be stored.
//
{
  int node;
  std::ofstream u_file;

  u_file.open ( u_file_name.c_str ( ) );

  if ( !u_file )
  {
    std::cout << "\n";
    std::cout << "SOLUTION_WRITE - Warning!\n";
    std::cout << "  Could not write the solution file \""
         << u_file_name << "\".\n";
    return;
  }

  for ( node = 0; node < node_num; node++ )
  {
    u_file << std::setw(14) << u[node] << "\n";
  }

  u_file.close ( );

  return;
}
//****************************************************************************80

void timestamp ( void )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    May 31 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  std::cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
//****************************************************************************80

void xy_set ( int nx, int ny, int node_num, double xl, double xr, double yb,
  double yt, double node_xy[] )

//****************************************************************************80
//
//  Purpose:
//
//    XY_SET sets the XY coordinates of the nodes.
//
//  Discussion:
//
//    The nodes are laid out in an evenly spaced grid, in the unit square.
//
//    The first node is at the origin.  More nodes are created to the
//    right until the value of X = 1 is reached, at which point
//    the next layer is generated starting back at X = 0, and an
//    increased value of Y.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 April 2004
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Parameters:
//
//    Input, int NX, NY, the number of elements in the X and
//    Y direction.
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double XL, XR, YB, YT, the X coordinates of
//    the left and right sides of the rectangle, and the Y coordinates
//    of the bottom and top of the rectangle.
//
//    Output, double NODE_XY[2*NODE_NUM], the nodes.
//
{
  int i;
  int j;

  for ( j = 0; j < 2*ny - 1; j++ )
  {
    for ( i = 0; i < 2*nx - 1; i++ )
    {
      node_xy[0+(i+j*(2*nx-1))*2] =
        ( double ( 2 * nx - i - 2 ) * xl
        + double (          i     ) * xr )
        / double ( 2 * nx     - 2 );

      node_xy[1+(i+j*(2*nx-1))*2] =
        ( double ( 2 * ny - j - 2 ) * yb
        + double (          j     ) * yt )
        / double ( 2 * ny     - 2 );
    }
  }
  return;
}