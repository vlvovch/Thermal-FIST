#ifndef HYPERSURFACESAMPLER_H
#define HYPERSURFACESAMPLER_H

/*
 * Thermal-FIST package
 *
 * Copyright (c) 2021 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */

#include "HRGBase/SplineFunction.h"
#include "HRGEventGenerator/MomentumDistribution.h"
#include "HRGEventGenerator/RandomGenerators.h"
#include "HRGEventGenerator/SphericalBlastWaveEventGenerator.h"
#include "HRGEventGenerator/CylindricalBlastWaveEventGenerator.h"

namespace thermalfist {

  struct ParticlizationHypersurfaceElement {
    double tau, x, y, eta;
    double dsigma[4];
    double u[4];
    double T, muB, muQ, muS;
    double edens, rhoB;
  };

  typedef std::vector<ParticlizationHypersurfaceElement> ParticlizationHypersurface;

  namespace RandomGenerators {

    /**
      * \brief Sample the volume element on a hypersurface from a multinomial distribution
      *
      */
    class VolumeElementSampler {
      std::vector<double> m_CumulativeProbabilities;
    public:
      VolumeElementSampler(const ParticlizationHypersurface* Hypersurface = NULL);
      VolumeElementSampler(const std::vector<double>& Weights);
      void FillProbabilities(const ParticlizationHypersurface* Hypersurface);
      void FillProbabilities(const std::vector<double>& Weights);
      void SetProbabilities(const std::vector<double>& CumulativeProbabilities) { m_CumulativeProbabilities = CumulativeProbabilities; }
      int SampleVolumeElement(MTRand& rangen = RandomGenerators::randgenMT) const;
    };


    /**
     * \brief Class for generating momentum of a particle from a hypersurface.
     *
     *  First chooses the hypersurface element multinomially, then sample the momentum from it.
     *  Stricly valid only for the ideal gas equation of state, but may be a sufficient approximation for
     *  interacting models in some cases as well.
     *
     */
    class HypersurfaceMomentumGenerator
      : public ParticleMomentumGenerator
    {
    public:

      /**
       * \brief Samples the Cartesian phase-space coordinates of a particle emmited from a hypersurface element.
       *
       * \param hypersurface    Pointer to a ParticlizationHypersurface object.
       * \param particle        Pointer to a ThermalParticle object representing the particle to sample.
       * \param mass            Particle mass in GeV. If negative, the pole/vacuum mass is used.
       * \param etasmear        The smear in longitudinal rapidity
       *
       * \return                A vector of 7 elements, the first 3 elements are the three-momentum (px,py,pz) in GeV,
       *                        the remaining four elements is the space-time coordinate (r0,rx,ry,rz) in fm/c
       */
      static std::vector<double> SamplePhaseSpaceCoordinateFromElement(
        const ParticlizationHypersurfaceElement* elem,
        const ThermalParticle* particle,
        const double& mass = -1.,
        const double& etasmear = 0.
      );

      /**
       * \brief Construct a new BoostInvariantMomentumGenerator object
       *
       * \param hypersurface    Pointer to a ParticlizationHypersurface object. Not deleted on destruction!
       * \param particle        Pointer to a ThermalParticle object. Not deleted on destruction!
       * \param positionsampler Pointer to a VolumeElementSampler object. Not deleted on destruction!
       * \param etasmear        The smear in longitudinal rapidity
       */
      HypersurfaceMomentumGenerator(
        const ParticlizationHypersurface* hypersurface = NULL,
        const ThermalParticle* particle = NULL,
        const VolumeElementSampler* positionsampler = NULL,
        double etasmear = 0.0);

      /**
       * \brief BoostInvariantMomentumGenerator desctructor.
       *
       * Will free the memory used by the object pointed to by m_FreezeoutModel
       */
      virtual ~HypersurfaceMomentumGenerator() {}

      double EtaSmear() const { return m_EtaSmear; }
      double Mass() const { return m_Particle->Mass(); }

      // Override functions begin

      virtual std::vector<double> GetMomentum(double mass = -1.) const;

      // Override functions end


    protected:

    private:
      const ParticlizationHypersurface* m_ParticlizationHypersurface;
      const ThermalParticle* m_Particle;
      const VolumeElementSampler* m_VolumeElementSampler;
      //ThermalMomentumGenerator m_Generator;
      //double m_Tkin;
      double m_EtaSmear;
    };


    /**
     * \brief Class for generating momentum of a particle
     *              from a longitudinally boost-invariant (2+1)-D hypersurface.
     *
     */
    class BoostInvariantHypersurfaceMomentumGenerator
      : public ParticleMomentumGenerator
    {
    public:


      /**
       * \brief Construct a new BoostInvariantMomentumGenerator object
       *
       * \param hypersurface Pointer to a ParticlizationHypersurface object. Not deleted on destruction!
       * \param hypersurface Pointer to a VolumeElementSampler object. Not deleted on destruction!
       * \param Tkin       The kinetic temperature (in GeV)
       * \param etamax     The longitudinal space-time rapidity cut-off
       * \param mass       Particle mass (in GeV)
       * \param statistics Quantum statistics (default: Maxwell-Boltzmann)
       * \param mu         Chemical potential (in GeV). Only matters for quantum statistics
       */
      BoostInvariantHypersurfaceMomentumGenerator(
        const ParticlizationHypersurface* hypersurface = NULL,
        VolumeElementSampler* positionsampler = NULL,
        double Tkin = 0.100, double etamax = 3.0, double mass = 0.938, int statistics = 0, double mu = 0);

      /**
       * \brief BoostInvariantMomentumGenerator desctructor.
       *
       * Will free the memory used by the object pointed to by m_FreezeoutModel
       */
      virtual ~BoostInvariantHypersurfaceMomentumGenerator() {}

      double EtaMax() const { return m_EtaMax; }
      double Mass() const { return m_Mass; }

      // Override functions begin

      virtual std::vector<double> GetMomentum(double mass = -1.) const;

      // Override functions end


    protected:

    private:
      const ParticlizationHypersurface* m_ParticlizationHypersurface;
      VolumeElementSampler* m_VolumeElementSampler;
      ThermalMomentumGenerator m_Generator;
      double m_Tkin;
      double m_EtaMax;
      double m_Mass;
    };

  }

  

  ///// \brief Class implementing the sampling of momenta of HRG particles
  /////        from an arbitrary hydro hypersurface
  /////
  /////        Calculates the multinomial weights for sampling the hypersurface element.
  /////        Samples the momenta of all the hadrons given the multiplcites in the current event externally.
  /////        
  class HypersurfaceEventGenerator : public EventGeneratorBase
  {
  public:

    /**
     * \brief Construct a new HypersurfaceEventGenerator object
     * 
     * \param hypersurface   A pointer to the particlization hypersurface. Not deleted at destruction!
     * \param model          A pointer to the thermal model object for calculating the densities at each hypersurface element.  Not deleted at destruction!
     * \param etasmear       Smearing in rapidity
     */
    HypersurfaceEventGenerator(
      const ParticlizationHypersurface* hypersurface = NULL,
      ThermalModelBase* model = NULL,
      double etasmear = 0.0) : EventGeneratorBase()
    {
      SetHypersurface(hypersurface);
      SetEtaSmear(etasmear);
      SetRescaleTmu();
      m_THM = model;
      //SetParameters(hypersurface, model, etasmear);
    }

    /**
     * \brief Construct a new HypersurfaceEventGenerator object
     *
     * \param TPS    A pointer to the particle list.  Not deleted at destruction!
     * \param config Event generator configuration
     * \param hypersurface   A pointer to the particlization hypersurface.  Not deleted at destruction!
     * \param etasmear       Smearing in rapidity
     */
    HypersurfaceEventGenerator(
      ThermalParticleSystem* TPS,
      const EventGeneratorConfiguration& config = EventGeneratorConfiguration(),
      const ParticlizationHypersurface* hypersurface = NULL,
      double etasmear = 0.0);

    virtual ~HypersurfaceEventGenerator() {}

    /// Override

    //virtual std::pair< std::vector<int>, double > SampleYields() const { exit(1); return std::pair< std::vector<int>, double >(); }
    //virtual SimpleEvent GetEvent(bool PerformDecays = true) const { exit(1); return SimpleEvent(); }
    virtual std::vector<double> GCEMeanYields() const;
    virtual std::vector<double>& GCEMeanYields();

    /// End override

    

    void SetModel(ThermalModelBase* model) { m_THM = model; m_ParametersSet = false; }

    void SetHypersurface(const ParticlizationHypersurface* hypersurface) { m_ParticlizationHypersurface = hypersurface; m_ParametersSet = false; }

    void SetEtaSmear(double etaSmear) { m_EtaSmear = etaSmear; m_ParametersSet = false; }
    double GetEtaSmear() const { return m_EtaSmear; }

    void SetRescaleTmu(bool rescale = false, double edens = 0.26);

    /// Sets the hypersurface parameters
    //void SetParameters(const ParticlizationHypersurface* hypersurface, ThermalModelBase* model, double etasmear = 0.0);
    //virtual void SetParameters();

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


    /// Sets the hypersurface parameters
    //virtual void CheckSetParameters() { if (!m_ParametersSet) SetParameters(); }

  protected:
    /// Sets up the random generators of particle momenta
    /// and resonances masses
    void SetMomentumGenerators();

    /// Processes the volume elements to calculate the multinomial volume element sampling probabilities and the full-space yields
    void ProcessVolumeElements();

    /// Calculates the (T,muB,muS,muQ) values as function of baryon density at fixed constant energy density
    static std::vector<std::vector<double>> CalculateTMuMap(ThermalModelBase* model, double edens, double rhomin = 0.0, double rhomax = 0.27, double drho = 0.001);

    /// Sets the hypersurface parameters
    //void SetParameters(const ParticlizationHypersurface* hypersurface, ThermalModelBase* model, double etasmear = 0.0);
    virtual void SetParameters();

    /// The computed grand-canonical yields in 4pi
    const std::vector<double>& FullSpaceYields() const { return m_FullSpaceYields; }
    
    //bool m_ParametersSet;

  private:
    const ParticlizationHypersurface* m_ParticlizationHypersurface;
    std::vector<RandomGenerators::VolumeElementSampler> m_VolumeElementSamplers;
    std::vector<double> m_FullSpaceYields;
    double m_EtaSmear;
    double m_Tav;
    std::vector<double> m_Musav;
    bool m_RescaleTmu;
    double m_edens;
    std::vector<SplineFunction> m_SplinesTMu;

    // Find T and muB = 0 to match energy density
    class BroydenEquationsTen : public BroydenEquations
    {
    public:
      BroydenEquationsTen(ThermalModelBase* model, double edens = 0.5) : BroydenEquations(),
        m_THM(model), m_edens(edens) {
        m_N = 1;
      }

      std::vector<double> Equations(const std::vector<double>& x) {
        const double& T = x[0];

        m_THM->SetTemperature(T);
        m_THM->SetBaryonChemicalPotential(0.);
        m_THM->SetElectricChemicalPotential(0.);
        m_THM->SetStrangenessChemicalPotential(0.);
        m_THM->CalculatePrimordialDensities();

        double en = m_THM->EnergyDensity();

        std::vector<double> ret(1, 0);
        ret[0] = (en / m_edens - 1.);

        return ret;
      }
    private:
      ThermalModelBase* m_THM;
      double m_edens;
    };

    // Find T,muB to match energy and baryon densities
    class BroydenEquationsTmuB : public BroydenEquations
    {
    public:
      BroydenEquationsTmuB(ThermalModelBase* model, double edens = 0.5, double rhoB = 0.16, int zeroMuBmode = 0) : BroydenEquations(),
        m_THM(model), m_edens(edens), m_rhoB(rhoB), m_zeroMuBmode(zeroMuBmode)
      {
        //m_N = 2;
        m_N = 4;
      }
      std::vector<double> Equations(const std::vector<double>& x) {
        const double& T = x[0];
        const double& muB = x[1];
        const double& muS = x[2];
        const double& muQ = x[3];

        //m_THM->SetQoverB(0.4);
        //m_THM->ConstrainMuQ(true);
        //m_THM->ConstrainMuS(true);

        m_THM->SetTemperature(T);
        if (!m_zeroMuBmode) {
          //m_THM->SetBaryonChemicalPotential(muB);
          //m_THM->ConstrainChemicalPotentials(false);
          m_THM->SetBaryonChemicalPotential(muB);
          m_THM->SetStrangenessChemicalPotential(muS);
          m_THM->SetElectricChemicalPotential(muQ);
          m_THM->CalculatePrimordialDensities();
        }
        else {
          m_THM->SetBaryonChemicalPotential(0.);
          m_THM->SetElectricChemicalPotential(0.);
          m_THM->SetStrangenessChemicalPotential(0.);
          m_THM->CalculatePrimordialDensities();
        }

        double en = m_THM->EnergyDensity();
        double rhoB = m_THM->BaryonDensity();

        std::vector<double> ret(4, 0);
        ret[0] = (en / m_edens - 1.);

        if (!m_zeroMuBmode) {
          ret[1] = rhoB / m_rhoB - 1.;
          ret[2] = m_THM->StrangenessDensity() / m_THM->AbsoluteStrangenessDensity();
          ret[3] = m_THM->ElectricChargeDensity() / m_THM->BaryonDensity() / 0.4 - 1.;
        }
        else {
          ret[1] = 0.;
          ret[2] = 0.;
          ret[3] = 0.;
        }
        return ret;
      }
    private:
      ThermalModelBase* m_THM;
      double m_edens, m_rhoB;
      int m_zeroMuBmode;
    };

    // returns (T,muB,muS,muQ) from given energy and baryon densities, assuming Q/B = 0.4, S = 0
    static std::vector<double> MatchEnergyBaryonDensities(ThermalModelBase* model, double edens, double rhoB) {
      if (edens == 0.0)
        return { 0.,0.,0.,0. };

      int zeromuBmode = 0;
      if (rhoB == 0.0)
        zeromuBmode = 1;

      if (!zeromuBmode) {
        BroydenEquationsTmuB eqs(model, edens, rhoB, zeromuBmode);
        Broyden broydn(&eqs);
        if (model->Parameters().muB == 0.0)
          model->SetBaryonChemicalPotential(0.010);
        //broydn.Solve({ model->Parameters().T, model->Parameters().muB });
        broydn.Solve({ model->Parameters().T, model->Parameters().muB, model->Parameters().muS, model->Parameters().muQ });
      }
      else {
        BroydenEquationsTen eqs(model, edens);
        Broyden broydn(&eqs);
        broydn.Solve({ model->Parameters().T });
      }

      return {
        model->Parameters().T,
        model->Parameters().muB,
        model->Parameters().muS,
        model->Parameters().muQ
      };
    }
  };

  ///// \brief Class implementing the sampling of momenta of HRG particles
  /////        from an arbitrary hydro hypersurface with excluded volume effect for baryons.
  /////
  /////        Corresponds to the EV-HRG model from https://arxiv.org/abs/1708.02852 and https://arxiv.org/abs/2107.00588
  /////        abd optimized accordingly
  /////        
  class HypersurfaceEventGeneratorEVHRG : public HypersurfaceEventGenerator
  {
  public:
    HypersurfaceEventGeneratorEVHRG(
      ThermalParticleSystem* TPS,
      const EventGeneratorConfiguration& config = EventGeneratorConfiguration(),
      const ParticlizationHypersurface* hypersurface = NULL,
      double etasmear = 0.0) : HypersurfaceEventGenerator(TPS, config, hypersurface, etasmear) {
      m_b = 0.0;
      m_rad = 0.0;
      m_MeanB = m_MeanAB = m_VEff = 0.0;
    }

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
    virtual SimpleEvent GetEvent(bool DoDecays = true) const;

    void SetExcludedVolume(double b) { m_b = b; m_ParametersSet = false; }
    void SetBaryonRadius(double r) { m_rad = r; m_ParametersSet = false; }
  protected:
    static double EVHRGWeight(int sampledN, double meanN, double V, double b);
    virtual void SetParameters();

    virtual SimpleEvent SampleParticles(const std::vector<int>& yields) const;

  private:

    std::pair<int, int> ComputeNBNBbar(const std::vector<int>& yields) const;

    double m_b, m_rad;
    double m_MeanB, m_MeanAB, m_VEff;
  };

  /// \brief Class implementing the Thermal Event Generator for
  ///        the boost-invariant (2+1)-d hydro hypersurface
  class BoostInvariantHypersurfaceEventGenerator : public EventGeneratorBase
  {
  public:

    /**
     * \brief Construct a new BoostInvariantHypersurfaceEventGenerator object
     *
     * \param TPS    A pointer to the particle list
     * \param config Event generator configuration
     * \param hypersurface   A pointer to the particlization hypersurface
     */
    BoostInvariantHypersurfaceEventGenerator(ThermalParticleSystem* TPS = NULL,
      const EventGeneratorConfiguration& config = EventGeneratorConfiguration(),
      double etamax = 0.5,
      const ParticlizationHypersurface* hypersurface = NULL);

    virtual ~BoostInvariantHypersurfaceEventGenerator();

    /// Sets the hypersurface parameters
    void SetParameters(double etamax, const ParticlizationHypersurface* hypersurface);

    double GetEtaMax() const { return m_EtaMax; }

  protected:
    /// Sets up the random generators of particle momenta
    /// and resonances masses
    void SetMomentumGenerators();

  private:
    double m_EtaMax;
    const ParticlizationHypersurface* m_ParticlizationHypersurface;
    RandomGenerators::VolumeElementSampler* m_VolumeElementSampler;
  };

} // namespace thermalfist

#endif