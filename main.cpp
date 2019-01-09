#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>
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
    
    const int num_params = 11;
    int k = 0;
    std::fstream myfile("config");
    std::string line;
    std::vector<std::string> params(num_params);
    while(k < num_params)
    {
        // Read line
        std::getline(myfile, line);
        std::cout << line << std::endl;
        // Read from line
        std::istringstream iss(line);
        // Split line into tokens
        std::vector<std::string> tokens;    
        std::copy(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(), back_inserter(tokens));
        // Get token
        params[k] = tokens[0];

        ++k;
    }
    // Read params
    int M = std::stoi(params[0]);
    int TN = std::stoi(params[1]);
    double T = std::stof(params[2]);
    int maxIters = std::stoi(params[3]);
    double eps = std::stof(params[4]);
    double mu = std::stof(params[5]);
    double b = std::stof(params[6]);
    double tau = std::stof(params[7]);
    double alpha = std::stof(params[8]);
    double coeff1 = std::stof(params[9]);
    double coeff2 = std::stof(params[10]);

    if (tau == 0)
    {
        tau = 2.0 / (2.0 * alpha + 1);
        std::cout << "Adaptive tau, tau = 2 / (2a + 1) = " << tau << std::endl;
    }

    if (tau < 0)
    {
        tau = fabs(tau) / pow(alpha, 0.2);
        std::cout << "Adaptive tau, tau = C / (a**0.2) = " << tau << std::endl;
    }
  
    SolverFDM solver(mu, b, M, T, TN, tau, alpha);

    k = 0;
    double err = 1.0;
    std::vector<double> phi((M+1)*(TN+1));
    std::vector<double> q((M+1)*(TN+1));
    std::vector<double> u_c(TN+1);
    
    solver.setInitialGuess(u_c, coeff1, coeff2);
    std::cout << "eps = " << eps << " max iters = " << maxIters << std::endl;
    while(err > eps && k < maxIters)
    {
        solver.solveForwardProblem(phi, u_c);
        solver.solveBackwardProblem(q, phi);
        err = solver.solveNextIteration(u_c, q);

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

    // Output
    std::ofstream results;

    results.open("u_c.txt");
    results << TN-1 << std::endl;
    for(int j = 1; j <= TN-1; ++j)
    {
        results << u_c[j] << std::endl;
    }
    results.close();

    results.open("u_c_solution.txt");
    results << TN-1 << std::endl;
    for(int j = 1; j <= TN-1; ++j)
    {
        results << solver.u_c_solution(j) << std::endl;
    }
    results.close();

    results.open("phi.txt");
    results << M+1 << std::endl; 
    results << TN+1 << std::endl;
    for(int j = 0; j <= TN; ++j)
        for(int i = 0; i <= M; ++i)
        {
            results << phi[solver.index(i, j)] << std::endl;
        }
    results.close();

    results.open("phi_solution.txt");
    results << M+1 << std::endl; 
    results << TN+1 << std::endl;
    for(int j = 0; j <= TN; ++j)
        for(int i = 0; i <= M; ++i)
        {
            results << solver.phi_solution(i, j) << std::endl;
        }
    results.close();

    results.open("final_phi.txt");
    results << M+1 << std::endl; 
    {
        int j = TN;
        for(int i = 0; i <= M; ++i)
        {
            results << phi[solver.index(i, j)] << std::endl;
        }
    }
    results.close();

    results.open("final_phi_solution.txt");
    results << M+1 << std::endl; 
    {
        int j = TN;
        for(int i = 0; i <= M; ++i)
        {
            results << solver.phi_solution(i, j) << std::endl;
        }
    }
    results.close();

    return 0;
}