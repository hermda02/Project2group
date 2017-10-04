//------Starting Project 2-----------------------------------------
//------first the non-interacting case-----------------------------
#include <iostream>
#include <new>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <armadillo>
#include <cstring>
#include <string>
#include <math.h>
#include <vector>
#include <fstream>
#include <time.h>
#include <ctime>

// #define CATCH_CONFIG_MAIN
// #include "catch.hpp"

using namespace std;
using namespace arma;

void offdiamax(mat A, int* p, int* q, int n){
    double maxm = 0;
    for (int i=0; i<n; i++){
        for (int j = i+1; j<n; j++){
            double aij = fabs(A(i,j));
            if (aij > maxm){
                maxm = aij; *p=i; *q=j;
                // cout << "offdiagmax:" << maxm << endl;
            }
        }
    }
//------------------------- Test to verify that the above code finds the non-diagonal element with the greatest distance from 0 -----------------------------------


        // intialize a new matrix with 0's on the diagonal so that the Armadillo function .max doesn't return a diagonal value
        // mat B = A;
        // for (int i=0; i<n; i++){
        //     B(i,i) = 0;
        // }

        // if ((maxm == B.max()) || (maxm = fabs(B.min()))){
        //      cout << "Maximum off diagonal successfully found!" << endl;
        // }
        // else {
        //     cout << "Failure!" << endl;
        // }
}

void jacobi(mat A, mat R, int k, int l, int n){
    double s, c;
    if (A(k,l) != 0.0){
        double t, tau;
        tau = (A(l,l)-A(k,k))/(2*A(k,l));
        if (tau >= 0) {
            t = 1.0/(tau + sqrt(1.0+tau*tau));
        } else {
            t = -1.0/(-tau + sqrt(1.0+tau*tau));
        }

        c = 1.0/sqrt(1+t*t);
        s = c*t;
    }
    else {
        c = 1.0;
        s = 0.0;
    }

    double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
    a_kk = A(k,k);
    a_ll = A(l,l);
    A(k,k) = c*c*a_kk - 2.0*c*s*A(k,l) + s*s*a_ll;
    A(l,l) = s*s*a_kk + 2.0*c*s*A(k,l) + c*c*a_ll;
    A(k,l) = 0.0;                                   // hard-coding non-diagonal elements by hand
    A(l,k) = 0.0;                                   // same here

    for ( int i = 0; i < n; i++ ) {
      if ( i != k && i != l ) {
        a_ik = A(i,k);
        a_il = A(i,l);
        A(i,k) = c*a_ik - s*a_il;
        A(k,i) = A(i,k);
        A(i,l) = c*a_il + s*a_ik;
        A(l,i) = A(i,l);
      }
    //  storing the new eigenvectors
    r_ik = R(i,k);
    r_il = R(i,l);

    R(i,k) = c*r_ik - s*r_il;
    R(i,l) = c*r_il + s*r_ik;
      }
}

int main(){

    // the essentials for our specific case [-d^2/d\rho^2 + wr^2*\rho^2 + 1/\rho]u(\rho) = \lambda u(\rho)
    int n = 400;                        // number of steps
    double rhomax = 10;                 // maximum radial distance
    double h = rhomax/(n+1);            // step size
    double hh = h*h;                    // step size squared
    double d = 2./(hh);                 // diagonal elements
    double e = -(1./(hh));              // off-diagonal elements
    double wr = 1;                      // the oscillator potential for the interacting case: 0.01, 0.5, 1, 5
    double wr2 = wr*wr;                 // w_r^2 as used in the potential


    // initializing the matricies and vectors we need

    vec rho(n);                         // the dimensionless radial distance rho_i
    vec v(n);                           // the potential v_i = omegar^2 + 1/rho_i

    mat A(n,n, fill::zeros);            // matrix to be solved
    mat R(n,n, fill::eye);              // the eigenvector matrix, initialized as the identity matrix


    // initializing the tridiagonal matrix for the differential equation discretization

    for (int i=0; i<n; i++) {
        rho(i) = (i+1)*h;
        v(i) = wr2*rho(i)*rho(i) + 1/rho(i);
        for (int j=0; j<n; j++){
            if (i == j){
                A(i,j) = d + v(i);
            }
            else if (i == j+1){
                A(i,j) = e;
            }
            else if (i == j-1){
                A(i,j) = e;
            };

        };
    };

    // cout << "Matrix A:" << A << endl;

    //-------------- Armadillo solver-------------------

    vec eigval;
    mat eigvec;

    eig_sym(eigval,eigvec,A);

    // cout << "Eigenvectors:" << endl << eigvec << endl;

    // cout << "Eigenvalues(arma):" << endl << eigval << endl;

    // the Armadillo solver will be used to verify that the eigenvalues of the matrix stays the same as it evolves

    //--------------------------------------------------


    // while loop to run the jacobi rotation until the eigenvalues are found

    double tolerance = 1.0E-10;
    int iterations = 0;
    int maxiter = 100;
    double maxm = 2.0E-10;
    while ( maxm > tolerance && iterations <= maxiter)
    {
        // find the maximum non-diagonal element which we will focus on for the rotation (removing the largest non-diagonal first)
        int p,q;
        offdiamax(A, &p, &q, n);
        jacobi(A,R,p,q,n);

        // do the matrix rotation focusing on the element discovered in the previous step:

        iterations++;

//------------------------------------------------------ Testing to show that the eigenvalues are preserved through rotation --------------------------------------

        // vec eval;
        // mat evec;
        // eig_sym(eval,evec,A);

        // cout << "eigval" << endl << eigval << endl;
        // cout << "post rotation" << endl << eval << endl;

        // int check = 0;

        // for ( int l = 0; l < n; l++) {
        //     if (eval(l) == eigval(l)){
        //         check++;
        //     }
        // }

        // if (check = n){
        //     cout << "Eigenvalues preserved!!" << endl;
        // }
        // else {
        //     cout << "BOOOOO!!!!" << endl;
        // }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------

    }

    vec S = sort(diagvec(A));
    cout << S << endl;

    ofstream myfile;
    myfile.open("eigenvalues.txt");
    myfile << "N = " << n << ", rho_max = " << rhomax << ", iterations = " << maxiter << endl;
    myfile << "Armadillo:" << setw(15) << "Jacobi:" << endl;

    for (int j = 0; j <=2; j++){
        myfile << eigval(j) << setw(15) << S(j) << endl;
    }

    myfile.close();
}

