//------Starting Project 2-----------------------------------------
//------the interacting case-----------------------------
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

#include "catch.hpp"

using namespace std;
using namespace arma;

int main(){

    // the essentials for our specific case [-d^2/d\rho^2 + wr^2*\rho^2 + 1/\rho]u(\rho) = \lambda u(\rho)
    int n = 10;                         // number of steps
    double rhomax = 10;                 // maximum radial distance
    double h = rhomax/(n+1);            // step size
    double hh = h*h;                    // step size squared
    double d = 2./(hh);                 // diagonal elements
    double e = -(1./(hh));              // off-diagonal elements
    double wr = 0.5;                    // the oscillator potential for the interacting case: 0.01, 0.5, 1, 5
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


    //-------------- Armadillo solver-------------------

    vec eigval;
    mat eigvec;

    eig_sym(eigval,eigvec,A);

    // cout << "Eigenvectors:" << endl << eigvec << endl;

    cout << "Eigenvalues:" << endl << eigval << endl;

    //--------------------------------------------------


    // while loop to run the jacobi rotation until the eigenvalues are found
    double tolerance = 1.0E-10;
    int iterations = 0;
    int maxiter = 5500;
    double maxm = 1;
    while ( maxm > tolerance && iterations <= maxiter)
    {

        // find the maximum non-diagonal element which we will focus on for the rotation (removing the largest non-diagonal first)
        int p,q;
        double maxm = 0;
        for (int i=0; i<n; i++){
            for (int j = i+1; j<n; j++){
                double aij = fabs(A(i,j));
                if (aij > maxm){
                    maxm = aij; p=i; q=j;
                }
            }
        }

        // do the matrix rotation focusing on the element discovered in the previous step:

        int k=p, l=q;

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
        } else {
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

        iterations++;
    }
    // cout << "The eigenvector matrix:" << endl << R << endl;
    // cout << "Eigenvalues (Jacobi):" << endl;
    // for (int i = 0; i<n; i++){
    //   cout << A(i,i) << endl;
    // }
}
