#include <iostream>
#include <string>
#include <vector>
#include <cassert>
#include <cmath> 

#ifndef SOLVER_H
#define SOLVER_H

class SolverFDM
{
private:
    // Physics
    double mu;
    double b;
    
    // Mesh config
    double A, B;
    double T;
    int M, TN;
    double h;
    double dt;

    // Solver options
    double tau;
    double alpha;
    
    // No copy
    SolverFDM(const SolverFDM& m);
    SolverFDM& operator=(const SolverFDM& m);

    // Mesh functions
    double w_c(int i) const
    {
        return 1.0;
    }

    double g_obs(int i) const
    {
        return 1.0;
    }

    double f(int i, int j) const
    {
        return 0.0;
    }

    double phi_obs(int j) const
    {
        double tj = j*dt;

        return (2.0/3.0)*(exp(-tj) - 3.0);
    }

    double phi_0(int i) const
    {
        double xi = A + i*h;
        return xi*xi - 1.0;
    }

    double q_TN(int i) const
    {
        return 0.0;
    }

    // Solutions
    double phi_solution(int i, int j) const
    {
        double xi = A + i*h;
        double tj = j*dt;

        return exp(-tj)*pow(xi, 2) - 1.0;
    }

    double u_c_solution(int j) const
    {
        double tj = j*dt;

        return -(1.0 + 2.0*mu*exp(-tj));
    }

public:
    SolverFDM(double mu, int M, double T, int TN, double tau, double alpha) : 
        mu(mu), M(M), T(T), TN(TN), tau(tau), alpha(alpha)
    {
        A = -1.0; B = 1.0;
        b = 1.0;

        h = (B - A) / M;
        dt = T / TN;

        std::cout << "Created solver. Options: A = " << A << " B = " << B << " T = " << T 
                    << " M = " << M << " TN = " << TN << std::endl;
        
        std::cout << "mu = " << mu << " b = " << b << std::endl;
        
        std::cout << "h = " << h << " dt = " << dt
            << " tau = " << tau << " alpha = " << alpha << std::endl;
    }

    // Index in vector for 2D arrays
    int index(int i, int j) const
    {
        return i + M*j;
    }

    int setInitialGuess(std::vector<double>& u_c) const;

    // Trapezodial Rule, integrate v(xi, tj)*g_obs(xi)*dx for fixed tj
    double integrate_g_obs(const std::vector<double>& v, int j) const;
    // Trapezodial Rule, integrate v(xi, tj)*w_c(xi)*dx for fixed tj
    double integrate_w_c(const std::vector<double>& v, int j) const;

    // Solve: L phi = f + B u_c
    int solveForwardProblem(std::vector<double>& phi, const std::vector<double>& u_c) const;
    
    // Solve: L* q = C*(C phi - phi_obs)
    int solveBackwardProblem(std::vector<double>& q, const std::vector<double>& phi) const;

    // Solve: alpha u_c + B* q = 0 
    int solveNextIteration(std::vector<double>& u_c, const std::vector<double>& q) const;

    // Errors
    double phi_full_errorL2(const std::vector<double>& phi) const;
    double phi_errorL2(const std::vector<double>& phi, int j) const;
    double u_c_full_errorL2(const std::vector<double>& u_c) const;
};

#endif