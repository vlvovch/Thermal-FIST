/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2018-2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef THERMALMODELVDWCANONICALSTRANGENESS_H
#define THERMALMODELVDWCANONICALSTRANGENESS_H

#include <map>
#include <iostream>

#include "HRGBase/ThermalModelCanonicalStrangeness.h"
#include "HRGVDW/ThermalModelVDW.h"

namespace thermalfist {

  /**
   * \brief Class implementing the quantum van der Waals
   *        model in the strangeness-canonical ensemble.
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
  class ThermalModelVDWCanonicalStrangeness : public ThermalModelCanonicalStrangeness
  {
  public:
    /**
     * \brief Construct a new Thermal ModelEVCanonicalStrangeness object
     * 
     * \param TPS A pointer to the ThermalParticleSystem object containing the particle list
     * \param params ThermalModelParameters object with current thermal parameters
     */
    ThermalModelVDWCanonicalStrangeness(ThermalParticleSystem *TPS, const ThermalModelParameters& params = ThermalModelParameters());

    /**
     * \brief Destroy the ThermalModelEVCanonicalStrangeness object
     * 
     */
    virtual ~ThermalModelVDWCanonicalStrangeness(void);
    
    /// \copydoc thermalfist::ThermalModelVDW::FillVirialEV()
    void FillVirialEV(const std::vector< std::vector<double> > & bij = std::vector< std::vector<double> >(0));

    // Override functions begin

    void FillVirial(const std::vector<double> & ri = std::vector<double>(0));

    void FillAttraction(const std::vector< std::vector<double> > & aij = std::vector< std::vector<double> >(0));
    
    void SetVirial(int i, int j, double b) { if (i >= 0 && i < static_cast<int>(m_Virial.size()) && j >= 0 && j < static_cast<int>(m_Virial[i].size())) m_Virial[i][j] = b; }
    
    void SetAttraction(int i, int j, double a) { if (i >= 0 && i < static_cast<int>(m_Attr.size()) && j >= 0 && j < static_cast<int>(m_Attr[i].size()))     m_Attr[i][j] = a; }

    double VirialCoefficient(int i, int j) const;

    double AttractionCoefficient(int i, int j) const;

    virtual void ReadInteractionParameters(const std::string &filename);

    virtual void WriteInteractionParameters(const std::string &filename);

    virtual void CalculateDensitiesGCE();

    virtual void CalculateEnergyDensitiesGCE();

    virtual void CalculatePressuresGCE();

    virtual void CalculatePrimordialDensities();

    virtual double CalculateEnergyDensity();

    virtual double CalculateEntropyDensity();

    virtual double CalculatePressure();

    virtual bool IsConservedChargeCanonical(ConservedCharge::Name charge) const;
    
    // Override functions end

  protected:
    /// Calculates the necessary auxiliary quantities
    /// in the QvdW model in the GCE consisting
    /// of non-strange particles
    void PrepareModelVDW();
    
    /// Clears m_modelVDW
    void ClearModelVDW();		
  
    /**
     * \brief Calculates the necessary 
     *        strangeness-canonical partition functions
     * 
     * Calculations are performed in a reduced strangeness correlation volume
     * due to the eigenvolume corrections from non-strange particles 
     * 
     * \param Vcs A vector of effective correlation volumes which is seen
     *            by each particle species
     */
    virtual void CalculateSums(const std::vector<double> &  Vcs);

    /// \copydoc thermalfist::ThermalModelVDW::MuShift()
    virtual double MuShift(int id) const;

    ThermalModelVDWFull *m_modelVDW; /**< Pointer to the QvdW model in the GCE with non-strange particles only */
    std::vector< std::vector<double> > m_Virial; /**Matrix of the excluded volume coefficients \f$ \tilde{b}_{ij} \f$ */
    std::vector< std::vector<double> > m_Attr;   /**Matrix of the attractive QvdW coefficients \f$ a_{ij} \f$ */
    double m_PNS;		                      /**< Pressure of all non-strange hadrons */
    std::vector<double> m_MuStar;         /**< Vector of the shifted chemical potentials */
    std::vector<double> m_Suppression;		/**< Common suppression factor, from non-strange hadrons */
  };

} // namespace thermalfist

#endif // THERMALMODELEVCANONICALSTRANGENESS_H
