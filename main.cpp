#include <iostream>
#include <string>
#include <vector>
#include <cassert>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <stdio.h>

#include "solver.hpp"

int main(int argc, char** argv)
{
    std::cout << "Integral Observations Problem" << std::endl;
    if (argc < 8) 
    {
        std::cout << "Program needs M, TN, maxIters, mu, T, tau, alpha" << std::endl;
        exit(0);
    }

    int M = atoi(argv[1]); 
    int TN = atoi(argv[2]);
    int maxIters = atoi(argv[3]);
    double mu = atof(argv[4]);
    double T = atof(argv[5]);
    double tau = atof(argv[6]);
    double alpha = atof(argv[7]);
    
    SolverFDM solver(mu, M, T, TN, tau, alpha);

    int k = 0;
    std::vector<double> phi((M+1)*(TN+1));
    std::vector<double> q((M+1)*(TN+1));
    std::vector<double> u_c(TN+1);
    
    solver.setInitialGuess(u_c);
    while(k < maxIters)
    {
        solver.solveForwardProblem(phi, u_c);
        solver.solveBackwardProblem(q, phi);
        solver.solveNextIteration(u_c, q);    

        ++k;
    }
    
    double phi_fullerr_L2, phi_err_L2, u_c_fullerr_L2;
    phi_fullerr_L2 = solver.phi_full_errorL2(phi);
    phi_err_L2 = solver.phi_errorL2(phi, TN);
    u_c_fullerr_L2 = solver.u_c_full_errorL2(u_c);

    std::cout << "Iterations: " << k << std::endl;
    std::cout << "phi full error (L2) : " << phi_fullerr_L2 << std::endl;
    std::cout << "phi error (L2) : " << phi_err_L2 << std::endl;
    std::cout << "u_c full error (L2) : " << u_c_fullerr_L2 << std::endl;

    /*
    std::ofstream results;
    results.open("results.txt");
    results << X;
    results.close();
    */

    return 0;
}