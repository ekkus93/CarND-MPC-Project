#ifndef MPC_H
#define MPC_H

#include <vector>
#include <cppad/cppad.hpp>
#include "Eigen-3.3/Eigen/Core"

using namespace std;
using CppAD::AD;

typedef CPPAD_TESTVECTOR(double) Dvector;

class FG_eval {
  public:
    // Fitted polynomial coefficients
    Eigen::VectorXd coeffs;
    FG_eval(Eigen::VectorXd coeffs);

    typedef CPPAD_TESTVECTOR(AD<double>) ADvector;

    void SetReferenceStateCost(ADvector &fg, const ADvector &vars);
    void SetupConstraints(ADvector &fg, const ADvector &vars);
    void operator()(ADvector& fg, const ADvector& vars);
  private:
    int cost_cte_factor_ = 3000;
    int cost_epsi_factor_ = 500; 
    int cost_v_factor_ = 1;
    int cost_current_delta_factor_ = 1;
    int cost_diff_delta_factor_ = 200;
    int cost_current_a_factor_ = 1;
    int cost_diff_a_factor_ = 1;      

    double ref_cte_ = 0;
    double ref_epsi_ = 0;
    double ref_v_ = 40;  

    const double Lf_ = 2.67;    
};

class MPC
{
  public:
    MPC();

    virtual ~MPC();

    void InitState(double x, double y, double psi,
                      double v, double cte,
                      double epsi, Dvector &vars);

    // Solve the model given an initial state and polynomial coefficients.
    // Return the first actuatotions.
    vector<double> Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs);
  private:
    size_t n_vars_ = 0;
    size_t n_constraints_ = 0;
};

#endif /* MPC_H */
