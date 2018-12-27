/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef BROYDEN_H
#define BROYDEN_H

#include <vector>

#include <Eigen/Dense>

namespace thermalfist {

  /**
  * Abstract class which implements the non-linear equations to be solved by the Broyden's method.
  */
  class BroydenEquations
  {
  public:
    /**
    * Default constructror.
    */
    BroydenEquations(void) { }

    /// Destructor.
    virtual ~BroydenEquations(void) { }

    /// Dimensions
    virtual int Dimension() const { return m_N; }

    void SetDimension(int dim) { m_N = dim; }

    /// Equation
    virtual std::vector<double> Equations(const std::vector<double> &x) = 0;

  protected:
    int m_N;           // Number of equations
  };


  /**
  * Class which implements calculation of the Jacobian needed for the Broyden's method.
  * Here it is done numerically using the BroydenEquations instance and the finite difference method.
  * A derived class may implement an analytic calculation.
  */
  class BroydenJacobian
  {
  public:
    static const double EPS; /// Finite difference for calculating J numerically

    /**
    * Default constructror.
    */
    BroydenJacobian(BroydenEquations *eqs = NULL) : m_Equations(eqs), m_dx(EPS) { }

    /// Destructor.
    virtual ~BroydenJacobian(void) { }

    /// Jacobian
    virtual Eigen::MatrixXd Jacobian(const std::vector<double> &x);

    void SetDx(double dx) { m_dx = dx; }
    double CurrentDx() const { return m_dx; }

  private:
    BroydenEquations *m_Equations;  // Pointer to BroydenEquations object
    double           m_dx;
  };

  /**
  * Implementation of the Broyden method to solve a system of non-linear equations.
  */
  class Broyden
  {
  public:
    class BroydenSolutionCriterium {
    public:
      BroydenSolutionCriterium(double relative_error = TOL) { m_RelativeError = relative_error; }
      virtual ~BroydenSolutionCriterium() { }
      virtual bool IsSolved(const std::vector<double>& x, const std::vector<double>& f, const std::vector<double>& xdelta = std::vector<double>()) const;
    protected:
      double m_RelativeError;
    };
  public:
    static const double TOL;       /// Default desired solution accuracy
    static const int    MAX_ITERS; /// Maximum number of Broyden iterations before terminating

    /**
    * Default constructror.
    */
    Broyden(BroydenEquations *eqs = NULL, BroydenJacobian *jaco = NULL) : m_Equations(eqs), m_Jacobian(jaco), m_Iterations(0), m_MaxDifference(0.), m_UseNewton(false) { }
    
    /// Destructor.
    virtual ~Broyden(void) { }

    virtual std::vector<double> Solve(const std::vector<double> &x0, BroydenSolutionCriterium *solcrit = NULL, int max_iterations = MAX_ITERS);
    //virtual std::vector<double> NewtonSolve(const std::vector<double> &x0, double relative_error = TOL, int max_iterations = MAX_ITERS);

    int Iterations() const { return m_Iterations; }
    int MaxIterations() const { return m_MaxIterations; }
    int MaxDifference() const { return m_MaxDifference; }

    void UseNewton(bool flag) { m_UseNewton = flag; }
    bool UseNewton() const { return m_UseNewton; }

  protected:
    BroydenEquations *m_Equations;  // Pointer to BroydenEquations object
    BroydenJacobian  *m_Jacobian;   // Pointer to BroydenJacobian object
    int m_Iterations;               // Current number of iterations
    int m_MaxIterations;            // Maximum number of iterations set
    double m_MaxDifference;         // Value of the maximum difference of equations from zero for the solution found
    bool m_UseNewton;
  };

} // namespace thermalfist

#endif
