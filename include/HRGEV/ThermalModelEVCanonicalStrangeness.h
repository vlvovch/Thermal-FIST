/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2018-2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef THERMALMODELEVCANONICALSTRANGENESS_H
#define THERMALMODELEVCANONICALSTRANGENESS_H

#include <map>
#include <iostream>

#include "HRGBase/ThermalModelCanonicalStrangeness.h"
#include "HRGEV/ThermalModelEVDiagonal.h"

namespace thermalfist {

  /**
   * \brief Class implementing the diagonal
   *        excluded-volume model in the strangeness-canonical ensemble.
   * 
   * ---
   * **NOTE**
   * 
   * The calculations are approximate and assume
   * that strange particles form a small part of the
   * total system, i.e. their contribution to the total
   * density/pressure is close to negligible.
   * Calculations may not be accurate if this condition
   * is not fulfilled.
   * 
   * ---
   * 
   */
  class ThermalModelEVCanonicalStrangeness : public ThermalModelCanonicalStrangeness
  {
  public:
    /**
     * \brief Construct a new Thermal ModelEVCanonicalStrangeness object
     * 
     * \param TPS A pointer to the ThermalParticleSystem object containing the particle list
     * \param params ThermalModelParameters object with current thermal parameters
     */
    ThermalModelEVCanonicalStrangeness(ThermalParticleSystem *TPS, const ThermalModelParameters& params = ThermalModelParameters());

    /**
     * \brief Destroy the ThermalModelEVCanonicalStrangeness object
     * 
     */
    virtual ~ThermalModelEVCanonicalStrangeness(void);

    /// \copydoc thermalfist::ThermalModelEVDiagonal::FillVirialEV()
    void FillVirialEV(const std::vector<double> & vi = std::vector<double>(0));

    /// \copydoc thermalfist::ThermalModelEVDiagonal::ExcludedVolume(int)
    double ExcludedVolume(int i) const;

    virtual void CalculateEnergyDensitiesGCE();

    virtual void CalculatePressuresGCE();

    // Override functions begin

    void SetRadius(double rad);
    
    void SetRadius(int i, double rad);

    /// \copydoc thermalfist::ThermalModelEVDiagonal::FillVirial(const std::vector<double> &)
    void FillVirial(const std::vector<double> & ri = std::vector<double>(0));

    /// \copydoc thermalfist::ThermalModelEVDiagonal::ReadInteractionParameters(const std::string &)
    virtual void ReadInteractionParameters(const std::string &filename);

    /// \copydoc thermalfist::ThermalModelEVDiagonal::WriteInteractionParameters(const std::string &)
    virtual void WriteInteractionParameters(const std::string &filename);

    virtual double CalculateEigenvolumeFraction();

    virtual void CalculateDensitiesGCE();

    virtual void CalculatePrimordialDensities();

    virtual double CalculateEnergyDensity();

    virtual double CalculateEntropyDensity();

    virtual double CalculatePressure();

    virtual bool IsConservedChargeCanonical(ConservedCharge::Name charge) const;

    // Override functions end

  protected:
    /// Calculates the necessary auxiliary quantities
    /// in the Diagonal EV model in the GCE consisting
    /// of non-strange particles
    void PrepareModelEV(); 

    /// Clears m_modelEV
    void ClearModelEV();    

    /// \copydoc thermalfist::ThermalModelEVDiagonal::MuShift()
    virtual double MuShift(int id) const;

    ThermalModelEVDiagonal *m_modelEV; /**< Pointer to the diagonal EV model in the GCE with non-strange particles only */
    std::vector<double> m_v;   /**< Vector of eigenvolumes of all hadrons */
    double m_PNS;              /**< Pressure of all non-strange hadrons */
    double m_Suppression;      /**< Common suppression factor, from non-strange hadrons */
    double m_EVNS;             /**< Total eigenvolume of all non-strange hadrons */
    double m_EVS;              /**< Total eigenvolume of all strange hadrons */
  };

} // namespace thermalfist

#endif // THERMALMODELEVCANONICALSTRANGENESS_H
