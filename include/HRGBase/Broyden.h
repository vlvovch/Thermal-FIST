/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2018-2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */

/**
 * \file Broyden.h
 * \brief Implementation of the generic Broyden's method routines
 * 
 */
#ifndef BROYDEN_H
#define BROYDEN_H

#include <cstddef>
#include <vector>

namespace thermalfist {

  /**
   * \brief Abstract class which defines the system of non-linear 
   *        equations to be solved by the Broyden's method.
   * 
   * It is presumed that all equations
   * are of the form 'l.h.s. == 0'.
   * Actual equations are to be specified
   * in a derived class.
   */
  class BroydenEquations
  {
  public:
    /// Default constructor. Does nothing.
    BroydenEquations(void) { }

    /// Destructor.
    virtual ~BroydenEquations(void) { }

    /// Number of equations.
    virtual int Dimension() const { return m_N; }

    /**
     * Set the number of equations.
     * 
     * \param dim The number of equations.
     */
    void SetDimension(int dim) { m_N = dim; }

    /**
     * Evaluates the l.h.s. of all the equations
     * being solved for the specified values of the variables. 
     * It is presumed that all equations
     * are of the form 'l.h.s. == 0'.
     * Pure virtual function which should be
     * implemented in a derived class.
     * 
     * \param x Vector of the variables' values.
     * \return std::vector<double> Vector of l.h.s.
     * computed.
     */
    virtual std::vector<double> Equations(const std::vector<double> &x) = 0;

  protected:
    /// The number of equations
    int m_N;           
  };


  /**
   * \brief Class which implements calculation of the Jacobian needed for the Broyden's method.
   * 
   * Here it is done numerically using the BroydenEquations instance and the finite difference method.
   * A derived class may implement an analytic calculation.
   */
  class BroydenJacobian
  {
  public:
    /// The default finite variable difference value
    /// used for calculations the Jacobian numerically.
    static const double EPS;
    /**
     * \brief Construct a new BroydenJacobian object
     * 
     * \param eqs A pointer to BroydenEquations object
     * which implements the equations to be solved by
     * the Broyden's method.
     */
    BroydenJacobian(BroydenEquations *eqs = NULL) : m_Equations(eqs), m_dx(EPS) { }

    /// Destructor.
    virtual ~BroydenJacobian(void) { }

    /**
     * \brief Evaluates the Jacobian for given values of the variables.
     * 
     * \param x Vector of the variables' values.
     * \return std::vector<double> The computed Jacobian matrix. 
     *         Returns a vector of matrix elements stored in a RowMajor ordering.
     * 
     * Here finite differences are used to approximates
     * the derivatives. This method can be overriden in
     * a derived class e.g. to implement analytic
     * calculations for a specific system of equations. 
     */
    virtual std::vector<double> Jacobian(const std::vector<double> &x);

    /**
     * \brief Set the finite variable difference value 
     *        used for calculating the Jacobian numerically.
     * 
     * \param dx The finite difference value.
     */
    void SetDx(double dx) { m_dx = dx; }

    /**
     * \return double Current finite variable difference value.
     */
    double CurrentDx() const { return m_dx; }

  private:
    /// Pointer to BroydenEquations object.
    BroydenEquations *m_Equations;  
    /// Finite variable difference value.
    double           m_dx;
  };

  /**
   * \brief Class implementing the Broyden method to solve a system of non-linear equations.
   */
  class Broyden
  {
  public:
    /**
     * \brief Sub-class where it is determined whether
     *        the required accuracy is achieved
     *        in the Broyden's method.
     * 
     * By default, the desired accuracy is achieved
     * if all of the equations deviate from zero by less
     * than the specified accuracy.
     * A derived class may implement a different criterium.
     */
    class BroydenSolutionCriterium {
    public:
      /**
       * Construct a new BroydenSolutionCriterium object
       * 
       * \param maximum_error The maximum deviation from zero
       * of any of the l.h.s. of all the equations for the desired
       * accuracy to be considered achieved.
       */
      BroydenSolutionCriterium(double maximum_error = TOL) { m_MaximumError = maximum_error; }
      
      /**
       * \brief Destroy the BroydenSolutionCriterium object
       * 
       */
      virtual ~BroydenSolutionCriterium() { }

      /**
       * Determines whether the solution with the
       * desired accuracy has been achieved.
       * 
       * \param x Values of the variables.
       * \param f Values of the l.h.s. of all the equations
       * \param xdelta Difference in variables' values between
       * current and previous iterations.
       * \return true The desired accuracy is achieved.
       * \return false The desired accuracy not achieved.
       */
      virtual bool IsSolved(const std::vector<double>& x, const std::vector<double>& f, const std::vector<double>& xdelta = std::vector<double>()) const;
    protected:
      /**
       * The maximum deviation from zero
       * of any of the l.h.s. of all the equations for the desired
       * accuracy to be considered achieved.
       */
      double m_MaximumError;
    };
  public:
    /// Default desired solution accuracy.
    static const double TOL;       
    /// Maximum number of Broyden iterations before terminating.
    static const int    MAX_ITERS; 

    /**
     * \brief Construct a new Broyden object
     * 
     * \param eqs A pointer to BroydenEquations object
     * which specifies the equations to be solved.
     * \param jaco A pointer to BroydenJacobian object
     * which calulates the Jacobian matrix. If nullptr is passed
     * then the Jacobian will be computed using finite
     * differences.
     */
    Broyden(BroydenEquations *eqs = NULL, BroydenJacobian *jaco = NULL) : m_Equations(eqs), m_Jacobian(jaco), m_Iterations(0), m_MaxDifference(0.), m_UseNewton(false) { }
    
    /**
     * \brief Destroy the Broyden object
     * 
     */
    virtual ~Broyden(void) { }

    /**
     * Solves the system of equations.
     * Will use at most max_iterations iterations.
     * Uses the Broyden's method by default, see [https://en.wikipedia.org/wiki/Broyden%27s_method](https://en.wikipedia.org/wiki/Broyden%27s_method)
     * Will use the Newton's method [https://en.wikipedia.org/wiki/Newton%27s_method](https://en.wikipedia.org/wiki/Newton%27s_method) 
     * if the corresponding flag was set with the UseNewton(bool) method.
     * 
     * \param x0 A vector of starting values of the variables.
     * \param solcrit Criterium used to determine whether the desired accuracy is achieved. 
     * \param max_iterations Maximum number of iterations before the Broyden's method terminates.
     * \return std::vector<double> A vector of the variables' values which wolve the equations.
     */
    virtual std::vector<double> Solve(const std::vector<double> &x0, BroydenSolutionCriterium *solcrit = NULL, int max_iterations = MAX_ITERS);

    /**
     * Returns number of Broyden/Newton iterations
     * used to obtain the solution.
     * Only meaningful to call this after the Solve()
     * method was called.
     * If Iterations() equals MaxIterations(),
     * then the convergence was not achieved.
     * 
     * \return Number of Broyden/Newton iterations.
     */
    int Iterations() const { return m_Iterations; }

    /**
     * \return Maximum number of Broyden/Newton iterations.
     */
    int MaxIterations() const { return m_MaxIterations; }

    /**
     * \return Maximum deviation from zero of any of the equations.
     */
    double MaxDifference() const { return m_MaxDifference; }

    /**
     * Specify whether to use Newton's method instead of
     * the Broyden's method.
     * In the Newton's method the Jacobian will be re-evaluated
     * from scrath at each iteration.
     * 
     * \param flag Use Newton's method if true, use Broyden's method otherwise.
     */
    void UseNewton(bool flag) { m_UseNewton = flag; }

    /**
     * \return true Newton's method is to be used
     * \return false Broyden's method is to be used
     */
    bool UseNewton() const { return m_UseNewton; }

  protected:
    BroydenEquations *m_Equations; 
    BroydenJacobian  *m_Jacobian;   
    int m_Iterations;              
    int m_MaxIterations;      
    double m_MaxDifference;         
    bool m_UseNewton;
  };

} // namespace thermalfist

#endif
