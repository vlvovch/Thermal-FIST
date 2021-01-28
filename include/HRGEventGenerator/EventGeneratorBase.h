/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef EVENTGENERATORBASE_H
#define EVENTGENERATORBASE_H


#include <sstream>

#include "HRGEventGenerator/SimpleEvent.h"
#include "HRGEventGenerator/Acceptance.h"
#include "HRGEventGenerator/RandomGenerators.h"
#include "HRGBase/xMath.h"
#include "HRGBase/ThermalModelBase.h"

namespace thermalfist {

  template <typename T>
  std::string to_string_fix(T value)
  {
    //create an output string stream
    std::ostringstream os;

    //throw the value into the string stream
    os << value;

    //convert the string stream into a string and return
    return os.str();
  }

  /// \brief Structure containing the thermal event generator configuration.
  struct EventGeneratorConfiguration {
    /// Enumerates the statistical ensembles 
    enum Ensemble { 
      GCE, ///< Grand-canonical
      CE,  ///< Canonical
      SCE, ///< Strangeness-canonical
      CCE  ///< Charm-canonical
    };

    /// Enumerates the different interaction models
    enum ModelType { 
      PointParticle, ///< Ideal gas
      DiagonalEV,    ///< Diagonal excluded-volume
      CrosstermsEV,  ///< Crossterms excluded-volume
      MeanFieldEV,   ///< Excluded-volume in the thermodynamic mean field approach (currently not used)
      QvdW           ///< Quantum van der Waals
    };
    
    /// The statistical ensemble used
    Ensemble fEnsemble;
    
    /// The type of interaction model
    ModelType fModelType;
    
    /// The chemical freeze-out parameters
    ThermalModelParameters CFOParameters;
    
    /// The total values of conserved charges in the CE
    int B, Q, S, C;

    /// Mixed-canonical configuration (full canonical by default)
    bool CanonicalB;
    bool CanonicalQ;
    bool CanonicalS;
    bool CanonicalC;

    /// The matrix of excluded volume coefficients \f$ \tilde{b}_{ij} \f$
    std::vector< std::vector<double> > bij;

    /// The matrix of van der Waals attraction coefficients \f$ a_{ij} \f$
    std::vector< std::vector<double> > aij;

    /// Whether partial chemical equilibrium (PCE) is used
    bool fUsePCE;

    /// PCE chemical potentials
    std::vector<double> fPCEChems;

    /// Default configuration
    EventGeneratorConfiguration();
  };

  /// \brief Base class for generating events with the Thermal Event Generator
  class EventGeneratorBase
  {
  public:
    /// Constructor
    EventGeneratorBase() { m_THM = NULL; fCEAccepted = fCETotal = 0; }

    /// Destructor
    virtual ~EventGeneratorBase();

    /// Clears the momentum generators for all particles
    void ClearMomentumGenerators();

    /// Sets the projectile laboratory kinetic energy per nucleon of the collision
    void SetCollisionKineticEnergy(double ekin) {
      SetCollisionCMSEnergy(sqrt(2.*xMath::mnucleon()*(ekin + 2. * xMath::mnucleon())));
    }

    /// Sets the projectile laboratory energy per nucleon of the collision
    void SetCollisionLabEnergy(double elab) {
      SetCollisionCMSEnergy(sqrt(2.*xMath::mnucleon()*(elab + xMath::mnucleon())));
    }

    /// Sets the center of mass energy \f$ \sqrt{s_{_{NN}}} \f$ of the collision
    void SetCollisionCMSEnergy(double ssqrt) {
      m_ssqrt = ssqrt;
      m_ekin = m_ssqrt * m_ssqrt / 2. / xMath::mnucleon() - 2. * xMath::mnucleon();
      m_elab = xMath::mnucleon() + m_ekin;
      double plab = sqrt(m_elab*m_elab - xMath::mnucleon() * xMath::mnucleon());
      m_ycm = 0.5 * log((m_elab + xMath::mnucleon() + plab) / (m_elab + xMath::mnucleon() - plab));
    }

    /// The y-pT acceptance map (not used by default).
    //std::vector<Acceptance::AcceptanceFunction>& GetAcceptance() { return m_acc; }

    /// Read the acceptance map from file.
    //virtual void ReadAcceptance(std::string accfolder);

    /// The center-of-mass longitudinal rapidity relative to the lab frame.
    double getYcm() const { return m_ycm; }

    /**
     * \brief Samples the primordial yields for each particle species.
     *
     * \return pair< std::vector<int>, double > The sampled yields. 
     *                                          The first element is a vector of the sampled yields.
     *                                          The second element is the weight.
     */
    virtual std::pair< std::vector<int>, double > SampleYields() const;

    /**
     * \brief Samples the momenta of the particles and returns the sampled list of particles as an event.
     *
     * The sampled SimpleEvent is assigned the weight of unity. 
     * This weight should be overriden if importance sampling is used.
     *
     * \param  yields       Vector of yields for each particle species for the given event.
     *                      Make sure the indices match the particle list pointed to by \ref m_THM.
     * \return SimpleEvent  The generated event containing the primordial particles.
     */
    virtual SimpleEvent SampleMomenta(const std::vector<int>& yields) const;

    /**
     * \brief Generates a single event.
     * 
     * \param PerformDecays If set to true, the decays of all particles 
     *                      marked unstable are performed until
     *                      only stable particles remain.
     *                      Otherwise only primordial particles are
     *                      generated and appear in the output
     * \return SimpleEvent  The generated event
     */
    virtual SimpleEvent GetEvent(bool PerformDecays = true) const;

    /**
     * \brief Performs decays of all unstable particles until only stable ones left.
     *
     * \param evtin An event structure contains the list of all the primordial particles.
     * \param TPS   Pointer to the particle list instance that contains all the decay properties.
     * \return      A SimpleEvent instance containing all particles after resonance decays.
     */
    static SimpleEvent PerformDecays(const SimpleEvent& evtin, ThermalParticleSystem* TPS);

    /**
     * \brief The grand-canonical mean yields.
     *
     * \return std::vector<double> The computed grand-canonical mean yields.
     */
    virtual std::vector<double> GCEMeanYields() const;

    /// Helper variable to monitor the Acceptance rate of the rejection
    /// sampling used for canonical ensemble and/or eigenvolumes.
    static int fCEAccepted, fCETotal;

    /**
     * \brief Set system volume.
     * 
     * Can be used to include volume fluctuations
     */
    void SetVolume(double V);

    /**
     * \brief Rescale the precalculated GCE means. 
     * 
     * Called when the system volume is changed
     */
    void RescaleCEMeans(double Vmod);

    /**
    * \brief Pointer to an underlying GCE Thermal Model.
    */
    ThermalModelBase* ThermalModel() { return m_THM; }


    double ComputeWeight(const std::vector<int>& totals) const;
    double ComputeWeightNew(const std::vector<int>& totals) const;

  protected:
    /**
     * \brief Sets the event generator configuration.
     * 
     * Must be called before generating any events.
     * 
     * \param TPS       Pointer to a particle list object
     * \param config    Event generator configuration
     */
    void SetConfiguration(ThermalParticleSystem *TPS,
      const EventGeneratorConfiguration& config);

    /// Prepares the parameters of multinomial distribution used
    /// for sampling the yields in the canonical ensemble
    void PrepareMultinomials();
    
    /// Samples the multiplicities of all the
    /// particle species from the given statistical ensemble
    /// \return A vector of the sampled multiplicities
    std::vector<int> GenerateTotals() const;

    /// Samples the multiplicities of all the
    /// particle species from the grand canonical ensemble
    /// \return A vector of the sampled multiplicities
    std::vector<int> GenerateTotalsGCE() const;

    /// Samples the multiplicities of all the
    /// particle species from the canonical ensemble
    ///
    /// Uses rejection sampling, and the multi-step
    /// procedure from F. Becattini, L. Ferroni, Eur. Phys. J. **C38**, 225 (2004) [hep-ph/0407117]
    /// \return A vector of the sampled multiplicities
    std::vector<int> GenerateTotalsCE() const;

    /// Samples the multiplicities of all the
    /// particle species from the strangeness-canonical ensemble
    ///
    /// Takes into account the case when the strangeness correlation volume
    /// is different from the total volume 
    /// \return A vector of the sampled multiplicities
    std::vector<int> GenerateTotalsSCE() const;

    /// Samples the multiplicities of all the
    /// particle species from the strangeness-canonical ensemble
    /// with the specified (sub)system volume
    ///
    /// \param VolumeSC The system volume
    /// \return A vector of the sampled multiplicities
    std::vector<int> GenerateTotalsSCESubVolume(double VolumeSC) const;

    /// Samples the multiplicities of all the
    /// particle species from the charm-canonical ensemble
    ///
    /// Takes into account the case when the strangeness correlation volume
    /// is different from the total volume 
    /// \return A vector of the sampled multiplicities
    std::vector<int> GenerateTotalsCCE() const;

    /// Samples the multiplicities of all the
    /// particle species from the charm-canonical ensemble
    /// with the specified (sub)system volume
    ///
    /// \param VolumeSC The (canonical) system volume
    /// \return A vector of the sampled multiplicities
    std::vector<int> GenerateTotalsCCESubVolume(double VolumeSC) const;

    EventGeneratorConfiguration m_Config;
    ThermalModelBase *m_THM;

    /// Ideal gas densities used for sampling an interacting HRG
    std::vector<double> m_DensitiesIdeal;

    /// Vector of momentum generators for each particle species
    std::vector<RandomGenerators::ParticleMomentumGenerator*>    m_MomentumGens;

    /// Vector of particle mass generators for each particle species
    /// Used if finite resonance widths are considered
    std::vector<RandomGenerators::ThermalBreitWignerGenerator*>  m_BWGens;

  private:

    /// Currently not used
    //static SimpleEvent PerformDecaysAlternativeWay(const SimpleEvent& evtin, ThermalParticleSystem* TPS);

    double m_ekin, m_ycm, m_ssqrt, m_elab;
    // Acceptance discontinued
    //std::vector<Acceptance::AcceptanceFunction> m_acc;

    //@{
    /// Indices and multinomial probabilities for an efficient CE sampling
    std::vector< std::pair<double, int> > m_Baryons;
    std::vector< std::pair<double, int> > m_AntiBaryons;
    std::vector< std::pair<double, int> > m_StrangeMesons;
    std::vector< std::pair<double, int> > m_AntiStrangeMesons;
    std::vector< std::pair<double, int> > m_ChargeMesons;
    std::vector< std::pair<double, int> > m_AntiChargeMesons;
    std::vector< std::pair<double, int> > m_CharmMesons;
    std::vector< std::pair<double, int> > m_AntiCharmMesons;
    std::vector< std::pair<double, int> > m_CharmAll;
    std::vector< std::pair<double, int> > m_AntiCharmAll;

    std::vector<double> m_BaryonsProbs;
    std::vector<double> m_AntiBaryonsProbs;
    std::vector<double> m_StrangeMesonsProbs;
    std::vector<double> m_AntiStrangeMesonsProbs;
    std::vector<double> m_ChargeMesonsProbs;
    std::vector<double> m_AntiChargeMesonsProbs;
    std::vector<double> m_CharmMesonsProbs;
    std::vector<double> m_AntiCharmMesonsProbs;
    std::vector<double> m_CharmAllProbs;
    std::vector<double> m_AntiCharmAllProbs;
    //@}

    double m_MeanB, m_MeanAB;
    double m_MeanSM, m_MeanASM;
    double m_MeanCM, m_MeanACM; 
    double m_MeanCHRMM, m_MeanACHRMM;
    double m_MeanCHRM, m_MeanACHRM;

    static double m_LastWeight;
    static double m_LastLogWeight;
    static double m_LastNormWeight;
  };

} // namespace thermalfist

#endif
