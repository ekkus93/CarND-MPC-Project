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
    double Lf; 
    double dt;

    // Fitted polynomial coefficients
    Eigen::VectorXd coeffs;
    FG_eval(Eigen::VectorXd coeffs, double Lf, double dt);

    typedef CPPAD_TESTVECTOR(AD<double>) ADvector;

    void SetReferenceStateCost(ADvector &fg, const ADvector &vars);
    void SetupConstraints(ADvector &fg, const ADvector &vars);
    void operator()(ADvector& fg, const ADvector& vars);
  private:
    int cost_cte_factor_ = 3000;
    int cost_epsi_factor_ = 600; 
    int cost_v_factor_ = 1;
    int cost_current_delta_factor_ = 1;
    int cost_diff_delta_factor_ = 200;
    int cost_current_a_factor_ = 1;
    int cost_diff_a_factor_ = 1;      

    double ref_cte_ = 0;     // reference cte
    double ref_epsi_ = 0;    // refenence epsi
    double ref_v_ = 40;      // reference velocity
};

class MPC
{
  public:
    double Lf; 
    double dt;

    MPC();

    virtual ~MPC();

    void InitVarBounds();

    void InitContraintBounds(double x, double y, double psi,
                              double v, double cte,
                              double epsi);

    void InitState(double x, double y, double psi,
                      double v, double cte,
                      double epsi, Dvector &vars);

    bool IpoptSolve(const Eigen::VectorXd coeffs, const Dvector &vars, 
                      vector<double> &result);

    // Solve the model given an initial state and polynomial coefficients.
    // Return the first actuatotions.
    vector<double> Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs);
  private:
    
    size_t n_vars_;          
    size_t n_constraints_;  
    Dvector vars_lowerbound_;
    Dvector vars_upperbound_;    
    Dvector constraints_lowerbound_;
    Dvector constraints_upperbound_;    
};

#endif /* MPC_H */
