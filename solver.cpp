#include "solver.hpp"

int SolverFDM::setInitialGuess(std::vector<double>& u_c, double coeff1, double coeff2) const
{
    assert(u_c.size() == TN+1);
    std::uniform_real_distribution<double> unif(-coeff2, coeff2);
    std::default_random_engine re;
    unif(re);

    for(int j = 0; j <= TN; ++j)
        u_c[j] = coeff1*u_c_solution(j) + unif(re);

    return 0;
}

int SolverFDM::solveForwardProblem(std::vector<double>& phi, const std::vector<double>& u_c) const
{
    int i, j;
    assert(phi.size() == (M+1)*(TN+1));
    assert(u_c.size() == TN+1);

    // Boundary conditions
    i = 0;
    for(j = 0; j <= TN; ++j)
        phi[index(i, j)] = 0.0;
    i = M;
    for(j = 0; j <= TN; ++j)
        phi[index(i, j)] = 0.0;

    // Initial conditions
    j = 0;
    for(i = 0; i <= M; ++i)
        phi[index(i, j)] = phi_0(i);

    // Inner area
    for(j = 1; j <= TN; ++j)
        for(i = 1; i <= M-1; ++i)
        {
            phi[index(i, j)] = phi[index(i, j-1)] + (mu*dt/h/h) * (phi[index(i+1, j-1)] - 2.0*phi[index(i, j-1)] + phi[index(i-1, j-1)]) 
                                - dt*b * phi[index(i, j-1)] + dt * w_c(i) * u_c[j-1] + dt * f(i, j-1);
        }
        
    return 0.0;
}

int SolverFDM::solveBackwardProblem(std::vector<double>& q, const std::vector<double>& phi) const
{
    int i, j;
    assert(q.size() == (M+1)*(TN+1));
    assert(phi.size() == (M+1)*(TN+1));

    // Boundary conditions
    i = 0;
    for(j = 0; j <= TN; ++j)
        q[index(i, j)] = 0.0;
    i = M;
    for(j = 0; j <= TN; ++j)
        q[index(i, j)] = 0.0;

    // Initial conditions
    j = TN;
    for(i = 0; i <= M; ++i)
        q[index(i, j)] = q_TN(i);
        //q[index(i, j)] = integrate_g_obs(phi, j) - phi_obs(j);

    // Inner area
    for(j = TN-1; j >= 0; --j)
        for(i = 1; i <= M-1; ++i)
        {
            double rhs = g_obs(i) * (integrate_g_obs(phi, j+1) - phi_obs(j+1));

            q[index(i, j)] = q[index(i, j+1)] + (mu*dt/h/h) * (q[index(i+1, j+1)] - 2.0*q[index(i, j+1)] + q[index(i-1, j+1)]) 
                                - dt*b * q[index(i, j+1)] + dt * rhs;
        }

    return 0;
}

double SolverFDM::solveNextIteration(std::vector<double>& u_c, const std::vector<double>& q) const
{
    int i, j;
    double err = 0.0;
    assert(q.size() == (M+1)*(TN+1));
    assert(u_c.size() == TN+1);

    for(j = 0; j <= TN; ++j)
    {
        double s = tau*(alpha * u_c[j] + integrate_w_c(q, j));
        err = err + dt*s*s;
        u_c[j] = u_c[j] - s;
    }

    return sqrt(err);
}

double SolverFDM::integrate_g_obs(const std::vector<double>& v, int j) const
{
    assert(v.size() == (M+1)*(TN+1));
    double sum = 0.0;    

    for(int i = 0; i <= M-1; ++i)
        sum = sum + g_obs(i) * h * 0.5 * (v[index(i, j)] + v[index(i + 1, j)]);

    return sum;
}

double SolverFDM::integrate_w_c(const std::vector<double>& v, int j) const
{
    assert(v.size() == (M+1)*(TN+1));
    double sum = 0.0;

    for(int i = 0; i <= M-1; ++i)
        sum = sum + w_c(i) * h * 0.5 * (v[index(i, j)] + v[index(i + 1, j)]);

    return sum;
}

double SolverFDM::phi_full_errorL2(const std::vector<double>& phi) const
{
    assert(phi.size() == (M+1)*(TN+1));
    int j;
    double err = 0.0;

    for(j = 0; j <= TN-1; ++j)
    {
        double s1 = phi_errorL2(phi, j);
        double s2 = phi_errorL2(phi, j+1);
        err = err + dt*0.5*(s1*s1 + s2*s2);
    }

    return sqrt(err);
}

double SolverFDM::phi_errorL2(const std::vector<double>& phi, int j) const
{
    assert(phi.size() == (M+1)*(TN+1));
    double err = 0.0;

    for(int i = 0; i <= M-1; ++i)
    {
        double s1 = phi[index(i, j)] - phi_solution(i, j);
        double s2 = phi[index(i+1, j)] - phi_solution(i+1, j);
        err = err + h*0.5*(s1*s1 + s2*s2);
    }

    return sqrt(err);
}

double SolverFDM::u_c_full_errorL2(const std::vector<double>& u_c) const
{
    assert(u_c.size() == TN+1);
    int i, j;
    double err = 0.0;

    for(j = 0; j <= TN-1; ++j)
    {
        double s1 = u_c[j] - u_c_solution(j);
        double s2 = u_c[j+1] - u_c_solution(j+1);
        err = err + dt*0.5*(s1*s1 + s2*s2);
    }

    return sqrt(err);
}