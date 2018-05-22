#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

// TODO: Set the timestep length and duration
const size_t N = 10;     // time step length
const double DT = 0.1;   // delta t

const double LF = 2.67;

size_t x_start = 0;
size_t y_start = x_start + N;
size_t psi_start = y_start + N;
size_t v_start = psi_start + N;
size_t cte_start = v_start + N;
size_t epsi_start = cte_start + N;
size_t delta_start = epsi_start + N;
size_t a_start = delta_start + N - 1;

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.

FG_eval::FG_eval(Eigen::VectorXd coeffs, double Lf, double dt)
{  
  this->Lf = Lf;
  this->dt = dt;
  this->coeffs = coeffs; 
}

void FG_eval::SetReferenceStateCost(ADvector &fg, const ADvector &vars)
{
  fg[0] = 0;

  for (unsigned int i = 0; i < N; i++)
  {
    fg[0] += cost_cte_factor_ * CppAD::pow(vars[cte_start + i] - ref_cte_, 2);
    fg[0] += cost_epsi_factor_ * CppAD::pow(vars[epsi_start + i] - ref_epsi_, 2);
    fg[0] += cost_v_factor_ * CppAD::pow(vars[v_start + i] - ref_v_, 2);
  }

  for (unsigned int i = 0; i < N - 1; i++)
  {
    fg[0] += cost_current_delta_factor_ * CppAD::pow(vars[delta_start + i], 2);
    fg[0] += cost_current_a_factor_ * CppAD::pow(vars[a_start + i], 2);
  }

  for (unsigned int i = 0; i < N - 2; i++)
  {
    fg[0] += cost_diff_delta_factor_ * CppAD::pow(vars[delta_start + i] - vars[delta_start + i], 2);
    fg[0] += cost_diff_a_factor_ * CppAD::pow(vars[a_start + i + 1] - vars[a_start + i], 2);
  }
}

void FG_eval::SetupConstraints(ADvector &fg, const ADvector &vars)
{
  // Initial constraints
  //
  // We add 1 to each of the starting indices due to cost being located
  // index at 0 of `fg`
  // This bumps up the position of all the other values.
  fg[1 + x_start] = vars[x_start];
  fg[1 + y_start] = vars[y_start];
  fg[1 + psi_start] = vars[psi_start];
  fg[1 + v_start] = vars[v_start];
  fg[1 + cte_start] = vars[cte_start];
  fg[1 + epsi_start] = vars[epsi_start];

  // The rest of the constraints
  for (unsigned int i = 0; i < N - 1; i++)
  {
    // t + 1
    AD<double> x1 = vars[x_start + i + 1];
    AD<double> y1 = vars[y_start + i + 1];
    AD<double> psi1 = vars[psi_start + i + 1];
    AD<double> v1 = vars[v_start + i + 1];
    AD<double> cte1 = vars[cte_start + i + 1];
    AD<double> epsi1 = vars[epsi_start + i + 1];

    // t
    AD<double> x0 = vars[x_start + i];
    AD<double> y0 = vars[y_start + i];
    AD<double> psi0 = vars[psi_start + i];
    AD<double> v0 = vars[v_start + i];
    AD<double> cte0 = vars[cte_start + i];
    AD<double> epsi0 = vars[epsi_start + i];

    AD<double> delta0 = vars[delta_start + i];
    AD<double> a0 = vars[a_start + i];

    AD<double> f0 = coeffs[0] + coeffs[1] * x0 + coeffs[2] * pow(x0, 2) + coeffs[3] * pow(x0, 3);
    AD<double> psides0 = CppAD::atan(coeffs[1] + (2 * coeffs[2] * x0) + (3 * coeffs[3] * pow(x0, 2)));

    // Here's `x` to get you started.
    // The idea here is to contraint this value to be 0.
    //
    // TODO: Setup the rest of the model constraints
    fg[2 + x_start + i] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
    fg[2 + y_start + i] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
    fg[2 + psi_start + i] = psi1 - (psi0 - v0 * delta0 / Lf * dt);
    fg[2 + v_start + i] = v1 - (v0 + a0 * dt);
    fg[2 + cte_start + i] = cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * dt));
    fg[2 + epsi_start + i] = epsi1 - ((psi0 - psides0) - v0 * delta0 / Lf * dt);
  }
}

void FG_eval::operator()(ADvector &fg, const ADvector &vars)
{
  // TODO: implement MPC
  // `fg` a vector of the cost constraints, `vars` is a vector of variable values (state & actuators)
  // NOTE: You'll probably go back and forth between this function and
  // the Solver function below.


  // Reference State Cost
  // TODO: Define the cost related the reference state and
  // anything you think may be beneficial.
  SetReferenceStateCost(fg, vars);

  // Setup Constraints
  //
  // NOTE: In this section you'll setup the model constraints.
  SetupConstraints(fg, vars);
}

//
// MPC class definition implementation.
//
MPC::MPC() {
  Lf = LF; 
  dt = DT;

  // TODO: Set the number of model variables (includes both states and inputs).
  // For example: If the state is a 4 element vector, the actuators is a 2
  // element vector and there are 10 timesteps. The number of variables is:
  //
  // 4 * 10 + 2 * 9
  n_vars_ = N * 6 + 2*(N - 1);
  // TODO: Set the number of constraints
  n_constraints_ = N * 6;

  // TODO: Set lower and upper limits for variables.
  InitVarBounds();
}
MPC::~MPC() {}

void MPC::InitVarBounds()
{
  vars_lowerbound_ = Dvector(n_vars_);
  vars_upperbound_ = Dvector(n_vars_);   

  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  for (unsigned int i = 0; i < n_constraints_; i++) {
    vars_lowerbound_[i] = -1.0e19;
    vars_upperbound_[i] = 1.0e19;
  }

  // The upper and lower limits of delta are set to -25 and 25
  // degrees (values in radians)
  // NOTE: Feel free to change this to something else.
  for (unsigned int i = delta_start; i < a_start; i++) {
    vars_upperbound_[i] = M_PI/8; 
    vars_lowerbound_[i] = -M_PI/8;
  }

  // Acceleration/deceleration upper and lower limits 
  // NOTE: Feel free to change this to something else.
  for (unsigned int i = a_start; i < n_vars_; i++) {
    vars_lowerbound_[i] = -1.0;
    vars_upperbound_[i] = 1.0;
  }
}

void MPC::InitContraintBounds(double x, double y, double psi,
                              double v, double cte,
                              double epsi)
{
  constraints_lowerbound_ = Dvector(n_constraints_);
  constraints_upperbound_ = Dvector(n_constraints_);   

  // Lower and upper limits for the constraints
  // Should be 0 besides initial state
  for (unsigned int i = 0; i < n_constraints_; i++) {
    constraints_lowerbound_[i] = 0;
    constraints_upperbound_[i] = 0;
  }

  constraints_lowerbound_[x_start] = x;
  constraints_lowerbound_[y_start] = y;
  constraints_lowerbound_[psi_start] = psi;
  constraints_lowerbound_[v_start] = v;
  constraints_lowerbound_[cte_start] = cte;
  constraints_lowerbound_[epsi_start] = epsi; 

  constraints_upperbound_[x_start] = x;
  constraints_upperbound_[y_start] = y;
  constraints_upperbound_[psi_start] = psi;
  constraints_upperbound_[v_start] = v;
  constraints_upperbound_[cte_start] = cte;
  constraints_upperbound_[epsi_start] = epsi;
}

void MPC::InitState(double x, double y, double psi,
                      double v, double cte,
                      double epsi, Dvector &vars)
{
  for (unsigned int i = 0; i < n_vars_; i++) {
    vars[i] = 0;
  }
  
  vars[x_start] = x;
  vars[y_start] = y;
  vars[psi_start] = psi;
  vars[v_start] = v;
  vars[cte_start] = cte;
  vars[epsi_start] = epsi;  
}

bool MPC::IpoptSolve(const Eigen::VectorXd coeffs, const Dvector &vars, 
                      vector<double> &result)
{
  bool ok = true;

  // object that computes objective and constraints
  FG_eval fg_eval(coeffs, Lf, dt);

  //
  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";
  
  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;
  
  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
                                        options, vars, vars_lowerbound_, vars_upperbound_, 
                                        constraints_lowerbound_, constraints_upperbound_, 
                                        fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;
  
  if (!ok) {
    // TODO: Maybe this should do more than just print an error message if not success.
    cout << "ERROR: solution status was not success (" << solution.status << ")" << endl;
  }

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

  // TODO: Return the first actuator values. The variables can be accessed with
  // `solution.x[i]`.
  //
  // {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
  // creates a 2 element double vector.

  result.push_back(solution.x[delta_start]);
  result.push_back(solution.x[a_start]);

  for (unsigned int i = 0; i < N-1; i++)
  {
    /*
    cout << "###" << i << "\n";
    cout << "###x_start * i + 1: " << (x_start + i + 1) << "\n";
    cout << "### x: " << solution.x[x_start + i + 1] << "\n";
    cout << "###y_start * i + 1: " << (y_start + i + 1) << "\n";
    cout << "### y: " << solution.x[y_start + i + 1] << "\n";
    */
    result.push_back(solution.x[x_start + i + 1]);
    result.push_back(solution.x[y_start + i + 1]);
  }

  return ok;
}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  // unpack state variables
  double x = state[0];
  double y = state[1];
  double psi = state[2];
  double v = state[3];
  double cte = state[4];
  double epsi = state[5];  

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars_);
  InitState(x, y, psi, v, cte, epsi, vars);

  //cout << "###vars: " << vars << endl;

  InitContraintBounds(x, y, psi, v, cte, epsi);

  vector<double> result;

  bool ok = MPC::IpoptSolve(coeffs, vars, result);

  if (!ok) {
    // TODO: Maybe this should do more than just print an error message if not success.
    cout << "ERROR: something bad happened" << endl;
  }
 
  return result;
}



                
