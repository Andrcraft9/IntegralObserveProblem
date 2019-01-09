#include <iostream>
#include <string>
#include <vector>
#include <cassert>
#include <cmath> 
#include <random>

#ifndef SOLVER_H
#define SOLVER_H

#define GAMMA 100.0

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
    
    // tmp
    //std::vector<double> rand_nums;

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
        double xi = A + i*h;
        double tj = j*dt;

        //return 0.0;
        //return rand_nums[index(i, j)] + (sin(M_PI * xi) * (GAMMA + mu * GAMMA * tj * pow(M_PI, 2) + b * GAMMA * tj) - GAMMA * tj);
        return (sin(M_PI * xi) * (GAMMA + mu * GAMMA * tj * pow(M_PI, 2) + b * GAMMA * tj) - GAMMA * tj);
    }

    double phi_obs(int j) const
    {
        double tj = j*dt;

        //return -GAMMA*(4.0/3.0)*exp(-tj*b);
        return 2 * GAMMA * tj / M_PI;
    }

    double phi_0(int i) const
    {
        double xi = A + i*h;
        
        //return GAMMA*(xi*xi - 1.0);
        return 0.0;
    }

    double q_TN(int i) const
    {
        return 0.0;
    }

public:
    SolverFDM(double mu, double b, int M, double T, int TN, double tau, double alpha) : 
        mu(mu), b(b), M(M), T(T), TN(TN), tau(tau), alpha(alpha)
    {
        /*
        rand_nums.resize((M+1)*(TN+1));
        std::uniform_real_distribution<double> unif(-1.0, 1.0);
        std::default_random_engine re;
        for(int j = 0; j <= TN; ++j)
            for(int i = 0; i <= M; ++i)
                rand_nums[index(i, j)] = (1000.0) * unif(re);
        */

        //A = -1.0; B = 1.0;
        //b = 1.0;
        A = 0.0; B = 1.0;
        //b = 0.5;

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

    int setInitialGuess(std::vector<double>& u_c, double coeff1, double coeff2) const;

    // Trapezodial Rule, integrate v(xi, tj)*g_obs(xi)*dx for fixed tj
    double integrate_g_obs(const std::vector<double>& v, int j) const;
    // Trapezodial Rule, integrate v(xi, tj)*w_c(xi)*dx for fixed tj
    double integrate_w_c(const std::vector<double>& v, int j) const;

    // Solve: L phi = f + B u_c
    int solveForwardProblem(std::vector<double>& phi, const std::vector<double>& u_c) const;
    
    // Solve: L* q = C*(C phi - phi_obs)
    int solveBackwardProblem(std::vector<double>& q, const std::vector<double>& phi) const;

    // Solve: alpha u_c + B* q = 0 
    double solveNextIteration(std::vector<double>& u_c, const std::vector<double>& q) const;

    // Errors
    double phi_full_errorL2(const std::vector<double>& phi) const;
    double phi_errorL2(const std::vector<double>& phi, int j) const;
    double u_c_full_errorL2(const std::vector<double>& u_c) const;

    // Solutions
    double phi_solution(int i, int j) const
    {
        double xi = A + i*h;
        double tj = j*dt;

        //return GAMMA*exp(-tj*b)*(pow(xi, 2) - 1.0);
        return GAMMA * tj * sin(xi * M_PI);
    }

    double u_c_solution(int j) const
    {
        double tj = j*dt;

        //return -GAMMA*2.0*mu*exp(-tj*b);
        return GAMMA * tj;
    }
};

#endif