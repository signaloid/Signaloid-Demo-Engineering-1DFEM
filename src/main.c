/*
  Licensing: This code is distributed under the GNU LGPL license.
  Author: John Burkardt
  Modified by: TODO

*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <uncertain.h>

int main(void);
void assemble(double diag_coeffs[], double left_coeffs[], double right_coeffs[], double f[],
    double x_grid_sizes[], int idx_of_unknown[], int n_basis, int node_indices[], int n_unknowns,
    int n_quadrature, int n_elements, double u_start, double u_end, double x_grid_nodes[], double x_quadrature[]);
double F_function(double x);
void build_mesh(double x_grid_sizes[], int ibc_mode, int idx_of_unknown[], int n_basis, int node_indices[],
    int n_elements, int* n_unknowns, double x_start, double x_grid_nodes[], double x_quadrature[], double x_end);
void init_problem(int* ibc_mode, int* n_quadrature, double* u_start, double* u_end, double* x_start, double* x_end);

void linear_basis_fn(int il, double x, double* phi_i, double* phi_ix, double x_lower, double x_upper);
void print_system(double diag_coeffs[], double left_coeffs[], double right_coeffs[], double f[], int n_unknowns);
void solve(double diag_coeffs[], double left_coeffs[], double right_coeffs[], double f[], int n_unknowns);
double P_function(double x);
double Q_function(double x);

void print_fem_solution(double f[], int ibc_mode, int idx_of_unknown[], int n_elements, int n_unknowns, double u_start,
    double u_end, double x_grid_nodes[], double u_fem[]);
void print_problem_summary(int ibc_mode, int n_quadrature, double u_start, double u_end, double x_start, double x_end);
void print_exact_solution(double x_start, double x_end, double x_grid_nodes[], double u_exact[]);
int main()
{
#define N_ELEMENTS 5
#define N_BASIS 2

    double diag_coeffs[N_ELEMENTS + 1];
    double left_coeffs[N_ELEMENTS + 1];
    double right_coeffs[N_ELEMENTS + 1];
    double f[N_ELEMENTS + 1];
    double x_grid_sizes[N_ELEMENTS];
    int ibc_mode;
    int idx_of_unknown[N_ELEMENTS + 1];
    int node_indices[N_BASIS * N_ELEMENTS];
    int n_quadrature;
    int n_unknowns;
    double u_start;
    double u_end;
    double x_start;
    double x_grid_nodes[N_ELEMENTS + 1];
    double x_quadrature[N_ELEMENTS];
    double x_end;

    double u_exact[N_ELEMENTS + 1];
    double u_fem[N_ELEMENTS + 1];

    init_problem(&ibc_mode, &n_quadrature, &u_start, &u_end, &x_start, &x_end);
    /*  Compute the geometric quantities.*/
    build_mesh(x_grid_sizes, ibc_mode, idx_of_unknown, N_BASIS, node_indices, N_ELEMENTS, &n_unknowns, x_start, x_grid_nodes, x_quadrature, x_end);

    /*  Assemble the linear system.*/
    assemble(diag_coeffs, left_coeffs, right_coeffs, f, x_grid_sizes, idx_of_unknown, N_BASIS, node_indices, n_unknowns, n_quadrature, N_ELEMENTS, u_start, u_end,
        x_grid_nodes, x_quadrature);

    /*  Print out the linear system.*/
    print_system(diag_coeffs, left_coeffs, right_coeffs, f, n_unknowns);

    /* Solve the linear system.*/
    solve(diag_coeffs, left_coeffs, right_coeffs, f, n_unknowns);

    /* Print out the solution.*/
    print_fem_solution(f, ibc_mode, idx_of_unknown, N_ELEMENTS, n_unknowns, u_start, u_end, x_grid_nodes, u_fem);
    print_exact_solution(x_start, x_end, x_grid_nodes, u_exact);

    /* Terminate.*/
    printf("\n");
    printf("FEM1D:");
    printf("  Normal end of execution.\n");

    printf("\n");

    return 0;
}
/******************************************************************************/

void init_problem(int* ibc_mode, int* n_quadrature, double* u_start, double* u_end, double* x_start, double* x_end)
{
    /* IBC declares what the boundary conditions are. */
    *ibc_mode = 1;
    /* NQUAD is the number of quadrature points per subinterval. The program as currently written cannot handle any value for NQUAD except 1. */
    *n_quadrature = 1;
    /* Set the values of U or U' at the endpoints. */
    *u_start = 0.0;
    *u_end = 0.0;
    /* Define the location of the endpoints of the interval. */
    *x_start = 0.0;
    *x_end = 1.0;
    return;
}
/******************************************************************************/
//+ d/dX (P.dU/dX) + F = Q.U

double P_function(double x)
{
    // 30mm^2, Rubber (0.01-0.1 GPa) Nylon (2-4 GPa) Steel (200 GPa)
    //
    double value;
    // double E = (0.1e9); //16.6mm
    // double E = (0.05e9); //33.33mm
    // double E = (0.01e9); //166.6mm
    double E = libUncertainDoubleUniformDist(0.01e9, 0.1e9);
    double A = (30.0e-6);
    value = A * E;
    return value;
}

double F_function(double x)
{

    double value;
    value = 100.0;
    return value;
}

double Q_function(double x)
{
    double value;
    value = 0;
    return value;
}

double exact_solution(double x, double x_end)
{
    return (F_function(x) / (2 * P_function(x))) * (2 * x * x_end - pow(x, 2));
}

/******************************************************************************/
void build_mesh(double x_grid_sizes[], int ibc_mode, int idx_of_unknown[], int n_basis, int node_indices[], int n_elements,
    int* n_unknowns, double x_start, double x_grid_nodes[], double x_quadrature[], double x_end)
{
    int i;
    /* Set the value of XN, the locations of the nodes. */
    printf("\n");
    printf(" Node      Location\n");
    printf("\n");
    for (i = 0; i <= n_elements; i++) {
        x_grid_nodes[i] = ((double)(n_elements - i) * x_start + (double)i * x_end) / (double)(n_elements);
        printf("  %8d  %14f \n", i, x_grid_nodes[i]);
    }
    /* Set the lengths of each subinterval. */
    printf("\n");
    printf("Element    Length\n");
    printf("\n");
    for (i = 0; i < n_elements; i++) {
        x_grid_sizes[i] = x_grid_nodes[i + 1] - x_grid_nodes[i];
        printf("  %8d  %14f\n", i + 1, x_grid_sizes[i]);
    }
    /* Set the quadrature points, each of which is the midpoint of its subinterval.*/
    printf("\n");
    printf("Element    Quadrature point\n");
    printf("\n");
    for (i = 0; i < n_elements; i++) {
        x_quadrature[i] = 0.5 * (x_grid_nodes[i] + x_grid_nodes[i + 1]);
        printf("  %8d  %14f\n", i + 1, x_quadrature[i]);
    }
    /* Set the value of NODE, which records, for each interval, the node_indices numbers at the left and right.*/
    printf("\n");
    printf("Element  Left Node  Right Node\n");
    printf("\n");
    for (i = 0; i < n_elements; i++) {
        node_indices[0 + i * 2] = i;
        node_indices[1 + i * 2] = i + 1;
        printf("  %8d  %8d  %8d\n", i + 1, node_indices[0 + i * 2], node_indices[1 + i * 2]);
    }
    /* Starting with node_indices 0, see if an unknown is associated with the node_indices.  If so, give it an index.*/
    *n_unknowns = 0;
    /*  Handle first node_indices.*/
    i = 0;
    if (ibc_mode == 1 || ibc_mode == 3) {
        idx_of_unknown[i] = -1;
    } else {
        *n_unknowns = *n_unknowns + 1;
        idx_of_unknown[i] = *n_unknowns;
    }
    /* Handle nodes 1 through n_elements-1 */
    for (i = 1; i < n_elements; i++) {
        *n_unknowns = *n_unknowns + 1;
        idx_of_unknown[i] = *n_unknowns;
    }
    /* Handle the last node_indices. */
    i = n_elements;

    if (ibc_mode == 2 || ibc_mode == 3) {
        idx_of_unknown[i] = -1;
    } else {
        *n_unknowns = *n_unknowns + 1;
        idx_of_unknown[i] = *n_unknowns;
    }
    printf("\n");
    printf("  Number of unknowns = %4d\n", *n_unknowns);
    printf("  Node  Unknown\n");
    printf("\n");
    for (i = 0; i <= n_elements; i++) {
        printf("  %8d  %8d\n", i, idx_of_unknown[i]);
    }

    return;
}
/******************************************************************************/

void assemble(double diag_coeffs[], double left_coeffs[], double right_coeffs[], double f[],
    double x_grid_sizes[], int idx_of_unknown[], int n_basis, int node_indices[], int n_unknowns, int n_quadrature,
    int n_elements, double u_start, double u_end, double x_grid_nodes[], double x_quadrature[])
{
    double aij;
    double he;
    int i;
    int ie;
    int ig;
    int il;
    int iq;
    int iu;
    int jg;
    int jl;
    int ju;
    double phi_i;
    double phi_ix;
    double phi_j;
    double phi_jx;
    double x;
    double x_lower;
    double x_quad;
    double x_upper;
    /*
Zero out the arrays that hold the coefficients of the matrix
and the right hand side.
*/
    for (i = 0; i < n_unknowns; i++) {
        f[i] = 0.0;
        diag_coeffs[i] = 0.0;
        left_coeffs[i] = 0.0;
        right_coeffs[i] = 0.0;
    }

    /*  For interval number IE,*/
    for (ie = 0; ie < n_elements; ie++) {
        he = x_grid_sizes[ie];
        x_lower = x_grid_nodes[node_indices[0 + ie * 2]];
        x_upper = x_grid_nodes[node_indices[1 + ie * 2]];

        /*  consider each quadrature point IQ,*/
        for (iq = 0; iq < n_quadrature; iq++) {
            x_quad = x_quadrature[ie];
            /*  and evaluate the integrals associated with the basis functions for the
       left, and for the right nodes. */
            for (il = 1; il <= n_basis; il++) {
                ig = node_indices[il - 1 + ie * 2];
                iu = idx_of_unknown[ig] - 1;

                if (0 <= iu) {
                    linear_basis_fn(il, x_quad, &phi_i, &phi_ix, x_lower, x_upper);
                    f[iu] = f[iu] + he * F_function(x_quad) * phi_i;
                    /*  Take care of boundary nodes at which U' was specified.*/
                    if (ig == 0) {
                        x = 0.0;
                        f[iu] = f[iu] - P_function(x) * u_start;
                    } else if (ig == n_elements) {
                        x = 1.0;
                        f[iu] = f[iu] + P_function(x) * u_end;
                    }
                    /*  Evaluate the integrals that take a product of the basis function
           times itself, or times the other basis function that is nonzero in
           this interval.*/
                    for (jl = 1; jl <= n_basis; jl++) {
                        jg = node_indices[jl - 1 + ie * 2];
                        ju = idx_of_unknown[jg] - 1;

                        linear_basis_fn(jl, x_quad, &phi_j, &phi_jx, x_lower, x_upper);

                        aij = he * (P_function(x_quad) * phi_ix * phi_jx + Q_function(x_quad) * phi_i * phi_j);
                        /* If there is no variable associated with the node_indices, then it's a
             specified boundary value, so we multiply the coefficient times
             the specified boundary value and subtract it from the right hand
             side. */
                        if (ju < 0) {
                            if (jg == 0) {
                                f[iu] = f[iu] - aij * u_start;
                            } else if (jg == n_elements) {
                                f[iu] = f[iu] - aij * u_end;
                            }
                        }
                        /* Otherwise, we add the coefficient we've just computed to the
               diagonal, or left or right entries of row IU of the matrix.*/
                        else {
                            if (iu == ju) {
                                diag_coeffs[iu] = diag_coeffs[iu] + aij;
                            } else if (ju < iu) {
                                left_coeffs[iu] = left_coeffs[iu] + aij;
                            } else {
                                right_coeffs[iu] = right_coeffs[iu] + aij;
                            }
                        }
                    }
                }
            }
        }
    }

    return;
}

void assemble_noQuadrature(double diag_coeffs[], double left_coeffs[], double right_coeffs[], double f[],
    double x_grid_sizes[], int idx_of_unknown[], int n_basis, int node_indices[], int n_unknowns, int n_quadrature,
    int n_elements, double u_start, double u_end, double x_grid_nodes[], double x_quadrature[])
{
    double aij;
    double he;
    int i;
    int ie;
    int ig;
    int il;
    int iq;
    int iu;
    int jg;
    int jl;
    int ju;
    double phi_i;
    double phi_ix;
    double phi_j;
    double phi_jx;
    double x;
    double x_lower;
    double x_quad;
    double x_upper;
    /*
Zero out the arrays that hold the coefficients of the matrix
and the right hand side.
*/
    for (i = 0; i < n_unknowns; i++) {
        f[i] = 0.0;
        diag_coeffs[i] = 0.0;
        left_coeffs[i] = 0.0;
        right_coeffs[i] = 0.0;
    }

    /*  For interval number IE,*/
    for (ie = 0; ie < n_elements; ie++) {
        he = x_grid_sizes[ie];
        x_lower = x_grid_nodes[node_indices[0 + ie * 2]];
        x_upper = x_grid_nodes[node_indices[1 + ie * 2]];

        /*  consider each quadrature point IQ,*/
        for (iq = 0; iq < n_quadrature; iq++) {
            x_quad = x_quadrature[ie];
            /*  and evaluate the integrals associated with the basis functions for the
       left, and for the right nodes. */
            for (il = 1; il <= n_basis; il++) {
                ig = node_indices[il - 1 + ie * 2];
                iu = idx_of_unknown[ig] - 1;

                if (0 <= iu) {
                    linear_basis_fn(il, x_quad, &phi_i, &phi_ix, x_lower, x_upper);
                    f[iu] = f[iu] + he * F_function(x_quad) * phi_i;
                    /*  Take care of boundary nodes at which U' was specified.*/
                    if (ig == 0) {
                        x = 0.0;
                        f[iu] = f[iu] - P_function(x) * u_start;
                    } else if (ig == n_elements) {
                        x = 1.0;
                        f[iu] = f[iu] + P_function(x) * u_end;
                    }
                    /*  Evaluate the integrals that take a product of the basis function
           times itself, or times the other basis function that is nonzero in
           this interval.*/
                    for (jl = 1; jl <= n_basis; jl++) {
                        jg = node_indices[jl - 1 + ie * 2];
                        ju = idx_of_unknown[jg] - 1;

                        linear_basis_fn(jl, x_quad, &phi_j, &phi_jx, x_lower, x_upper);

                        aij = he * (P_function(x_quad) * phi_ix * phi_jx + Q_function(x_quad) * phi_i * phi_j);
                        /* If there is no variable associated with the node_indices, then it's a
             specified boundary value, so we multiply the coefficient times
             the specified boundary value and subtract it from the right hand
             side. */
                        if (ju < 0) {
                            if (jg == 0) {
                                f[iu] = f[iu] - aij * u_start;
                            } else if (jg == n_elements) {
                                f[iu] = f[iu] - aij * u_end;
                            }
                        }
                        /* Otherwise, we add the coefficient we've just computed to the
               diagonal, or left or right entries of row IU of the matrix.*/
                        else {
                            if (iu == ju) {
                                diag_coeffs[iu] = diag_coeffs[iu] + aij;
                            } else if (ju < iu) {
                                left_coeffs[iu] = left_coeffs[iu] + aij;
                            } else {
                                right_coeffs[iu] = right_coeffs[iu] + aij;
                            }
                        }
                    }
                }
            }
        }
    }

    return;
}
/******************************************************************************/

void linear_basis_fn(int il, double x, double* phi_i, double* phi_ix, double x_lower,
    double x_upper)
{
    if (x_lower <= x && x <= x_upper) {
        if (il == 1) {
            *phi_i = (x_upper - x) / (x_upper - x_lower);
            *phi_ix = -1.0 / (x_upper - x_lower);
        } else {
            *phi_i = (x - x_lower) / (x_upper - x_lower);
            *phi_ix = 1.0 / (x_upper - x_lower);
        }
    }
    /* If X is outside of the interval, just set everything to 0. */
    else {
        *phi_i = 0.0;
        *phi_ix = 0.0;
    }

    return;
}
/******************************************************************************/

void solve(double diag_coeffs[], double left_coeffs[], double right_coeffs[], double f[], int n_unknowns)
{
    int i;
    /*  Carry out Gauss elimination on the matrix, saving information needed for the backsolve. */
    right_coeffs[0] = right_coeffs[0] / diag_coeffs[0];

    for (i = 1; i < n_unknowns - 1; i++) {
        diag_coeffs[i] = diag_coeffs[i] - left_coeffs[i] * right_coeffs[i - 1];
        right_coeffs[i] = right_coeffs[i] / diag_coeffs[i];
    }
    diag_coeffs[n_unknowns - 1] = diag_coeffs[n_unknowns - 1] - left_coeffs[n_unknowns - 1] * right_coeffs[n_unknowns - 2];
    /*  Carry out the same elimination steps on F that were done to the matrix. */
    f[0] = f[0] / diag_coeffs[0];
    for (i = 1; i < n_unknowns; i++) {
        f[i] = (f[i] - left_coeffs[i] * f[i - 1]) / diag_coeffs[i];
    }
    /* And now carry out the steps of "back substitution". */
    for (i = n_unknowns - 2; 0 <= i; i--) {
        f[i] = f[i] - right_coeffs[i] * f[i + 1];
    }

    return;
}
/******************************************************************************/

void print_problem_summary(int ibc_mode, int n_quadrature, double u_start, double u_end, double x_start, double x_end)
{
    printf("\n");
    printf("Modified FEM1D, by CTS, after John Burkardt.\n");

    printf("  Based on the code to solve the general two-point boundary value problem\n");
    printf("  - d/dX (P dU/dX) + Q U  =  F\n");

    printf("  P, Q and F specified to Solve dynamics of a fixed cantilever under constant axial force \n");
    printf("  given by the differential equation: d/dX (AE dU(x)/dx) + a = 0\n");
    printf("  i.e. P = AE, F = a, Q = 0 \n");
    printf("\n");
    printf("  on the interval [s_start,x_end], specifying the value of U or U' at each end.\n");
    printf("  The interval [s_start,x_end] is broken into N_ELEMENTS = %d elements\n", N_ELEMENTS);
    printf("  Number of basis functions per element is N_BASIS = %d\n", N_BASIS);

    printf("\n");
    printf("  The equation is to be solved for\n");
    printf("  X greater than x_start = %f and less than x_end = %f\n", x_start, x_end);
    printf("\n");
    printf("  The boundary conditions are:\n");

    if (ibc_mode == 1 || ibc_mode == 3) {
        printf("  At X = x_start, U = %f (from u_start)\n", u_start);
    } else {
        printf("  At X = x_start, U' = %f (from u_start)\n", u_start);
    }

    if (ibc_mode == 2 || ibc_mode == 3) {
        printf("  At X = x_end, U = %f (from u_end)\n", u_end);
    } else {
        printf("  At X = x_end, U' = %f (from u_end)\n", u_end);
    }

    printf("\n");
    printf("  Number of quadrature points per element is %d\n", n_quadrature);

    return;
}

void print_system(double diag_coeffs[], double left_coeffs[], double right_coeffs[], double f[], int n_unknowns)
{
    int i;

    printf("\n");
    printf("Printout of tridiagonal linear system:\n");
    printf("\n");
    printf("Equation  ALEFT  ADIAG  right_coeffs  RHS\n");
    printf("\n");

    for (i = 0; i < n_unknowns; i++) {
        printf("  %8d  %14f  %14f  %14f  %14f\n", i + 1, left_coeffs[i], diag_coeffs[i],
            right_coeffs[i], f[i]);
    }

    return;
}

void print_fem_solution(double f[], int ibc_mode, int idx_of_unknown[], int n_elements, int n_unknowns, double u_start,
    double u_end, double x_grid_nodes[], double u_fem[])
{
    int i;
    double u;

    printf("\n");
    printf("  ----- FEM solution -----\n");
    printf("  Node    X(I)        U(X(I)) [m]    U(X(I)) [mm]\n");

    for (i = 0; i <= n_elements; i++) {
        /* If we're at the first node_indices, check the boundary condition. */
        if (i == 0) {
            if (ibc_mode == 1 || ibc_mode == 3) {
                u = u_start;
            } else {
                u = f[idx_of_unknown[i] - 1];
            }
        }
        /* If we're at the last node_indices, check the boundary condition. */
        else if (i == n_elements) {
            if (ibc_mode == 2 || ibc_mode == 3) {
                u = u_end;
            } else {
                u = f[idx_of_unknown[i] - 1];
            }
        }
        /* Any other node_indices, we're sure the value is stored in F. */
        else {
            u = f[idx_of_unknown[i] - 1];
        }
        u_fem[i] = u;
        printf("  %3d     %8f      %8f      %8f\n", i, x_grid_nodes[i], u, u * 1000);
    }

    return;
}

void print_exact_solution(double x_start, double x_end, double x_grid_nodes[], double u_exact[])
{
    int i = 0;

    printf("\n");
    printf("  ---- Exact solution ----\n");
    printf("  Node    X(I)        U(X(I)) [m]    U(X(I)) [mm]\n");

    for (i = 0; i <= N_ELEMENTS; i++) {
        u_exact[i] = exact_solution(x_grid_nodes[i], x_end);
        printf("  %3d     %8f      %8f      %8f\n", i, x_grid_nodes[i], u_exact[i], u_exact[i] * 1000);
    }

    return;
}

