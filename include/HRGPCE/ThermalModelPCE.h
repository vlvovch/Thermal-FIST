#ifndef THERMALMODELPCE_H
#define THERMALMODELPCE_H
#include <map>

#include "HRGBase/ThermalModelBase.h"
#include "HRGBase/Broyden.h"

namespace thermalfist {

  /**
   * \brief Class implementing HRG in partial chemical equilibrium.
   * 
   * Partial chemical equilbiirum (PCE) describes the hadronic phase dynamics.
   *
   * The basic ideas about PCE can be found in
   * H. Bebie, P. Gerber, J.L. Goity, H. Leutwyler,
   * Nucl. Phys. B **378**, 95 (1992)
   * 
   * The present implementation was used (and described) in a paper
   * V. Vovchenko, K. Gallmeister, J. Schaffner-Bielich, C. Greiner,
   * Phys. Lett. B **800**, 135131 (2020),
   * [arXiv:1903.10024](https://arxiv.org/abs/1903.10024)
   * 
   * A typical usage pattern is the following:
   * 1. Create a ThermalModelBase() instance implementing a particular grand-canonical HRG model.
   * 2. Create a ThermalModelPCE() instance using the ThermalModelBase() instance from the first step.
   * Configure PCE options like which hadrons should be considered stable and whether the Saha
   * equation should be used for light nuclei.
   * 3. Set the chemical freeze-out conditions via SetChemicalFreezeout()
   * 4. Evaluate the PCE conditions at a given temperature or volume in the hadronic phase via CalculatePCE()
   * 5. Use the pointer ThermalModel() to the HRG model to access the various physical properties at the given
   * PCE point, like various hadron yields or the equation of state.
   * 
   */
  class ThermalModelPCE
  {
  public:
    /**
     * \brief Whether partial chemical equilibrium should be calculated
     * at a fixed value of the temperature or a fixed value of the volume.
     * 
     * Used in CalculatePCE()
     *
     */
    enum PCEMode {
      /// Partial chemical equilibrium at fixed value of the temperature
      AtFixedTemperature = 0, 
      /// Partial chemical equilibrium at fixed value of the volume
      AtFixedVolume = 1
    };


    /**
     * \brief Construct a new ThermalModelPCE object
     *
     * \param THMbase               A pointer to the ThermalModelBase object containing the HRG model implementation
     * \param FreezeLonglived       Whether long-lived resonance abundances should be frozen at Tch along with the stable hadrons
     * \param LonglivedResoWidthCut Threshold width for the long-lived resonances to be considered stable in the PCE
     */
    ThermalModelPCE(ThermalModelBase *THMbase, bool FreezeLonglived = false, double LonglivedResoWidthCut = 0.015);

    /**
     * \brief Destroy the ThermalModelPCE object
     *
     */
    virtual ~ThermalModelPCE(void) { }

    //@{
      /**
       * \brief Whether the nuclear abundances are evaluated through the Saha equation.
       *
       * If not set, the nuclear abundances are frozed at the chemical freeze-out of stable hadrons.
       *
       * \param flag Flag whether the Saha equation is used.
       */
    void UseSahaForNuclei(bool flag) { m_UseSahaForNuclei = flag; m_StabilityFlagsSet = false; }
    bool UseSahaForNuclei() const { return m_UseSahaForNuclei; }
    //@}

    //@{
      /**
       * \brief Whether long-lived resonances yields should be frozen.
       *
       * If set, the long-lived resonance with width smaller than LonglivedResonanceWidthCut() are frozen at the chemical freeze-out of stable hadrons.
       * Otherwise, their abundances are evaluated through the chemical potentials of their decay products.
       *
       * \param flag Whether long-lived resonances yields should be frozen.
       */
    void FreezeLonglivedResonances(bool flag) { m_FreezeLonglivedResonances = flag; m_StabilityFlagsSet = false; }
    bool FreezeLonglivedResonances() const { return m_FreezeLonglivedResonances; }
    //@}

    //@{
      /**
       * \brief The threshold resonance width value to consider the resonance long-lived and its abundance frozen in the hadronic phase.
       *
       * Note: This function sets FreezeLonglivedResonances() to true.
       *
       * \param width_cut The threshold resonance width [GeV].
       */
    void SetLonglivedResonanceWidthCut(double width_cut) { m_ResoWidthCut = width_cut; FreezeLonglivedResonances(true); }
    double LonglivedResonanceWidthCut() const { return m_ResoWidthCut; }
    //@}

    //@{
      /**
       * \brief Manually set the PCE stability flags for all species.
       *
       * \param StabilityFlags Vector of stability flags. 
       *                       Each element corresponds to a particle specie with the same 0-based index in the particle list.
       *                       Element i equal to zero means particle yield is not frozen in the hadronic phase, otherwise it is frozen.
       */
    virtual void SetStabilityFlags(const std::vector<int>& StabilityFlags);
    const std::vector<int>& StabilityFlags() const { return m_StabilityFlags; }
    //@}

    //@{
      /**
       * \brief Sets the chemical freeze-out conditions to be used as an initial condition for PCE calculations
       *
       * \param params   Thermal parameters at the chemical freeze-out.
       * \param ChemInit Chemical potentials of all species at the chemical freeze-out. 
       *                 If this vector is empty (default value), the chemical potentials are populated using \p params
       */
    void SetChemicalFreezeout(const ThermalModelParameters& params, const std::vector<double>& ChemInit = std::vector<double>(0));
    //@}

    //@{
      /**
       * \brief Sets the entropy density at the chemical freeze-out.
       *
       * \param sinit  Entropy density at the chemical freeze-out.
       */
    void SetEntropyDensityChem(double sinit) { m_EntropyDensityInit = sinit; }
    double EntropyDensityChem() const { return m_EntropyDensityInit; }
    //@}

    /**
     * \brief Pointer to the HRG model used in calculations.
     */
    ThermalModelBase* ThermalModel() const { return m_model; }

    /**
     * \brief Solves the equations of partial chemical equilibrium at a fixed temperature or a fixed volume
     *
     * \param param  Temperature (in GeV) if \p mode correspond to PCE at fixed temperature, 
     * Volume (in fm^3)  if \p mode correspond to PCE at fixed volume.
     * \param mode   Determines whether the PCE calculation
     * is perfored at a fixed temperature (default) or at a fixed volume. 
     */
    virtual void CalculatePCE(double param, PCEMode mode = AtFixedTemperature);

    /**
     * \return Vector of chemical potentials of all particles species, as resulted from the last CalculatePCE() call
     */
    const std::vector<double>& ChemicalPotentials() const { return m_ChemCurrent; }
          std::vector<double>& ChemicalPotentials() { return m_ChemCurrent; }

    /**
     * \return The system volume, as resulted from the last CalculatePCE() call
     */
    double Volume() const { return m_ParametersCurrent.V; }

    /**
     * \brief Fills the "decay" products of light nuclei in accordance with their baryon content
     * 
     * \param TPS  Pointer to the particle list
     */
    static void PrepareNucleiForPCE(ThermalParticleSystem *TPS);

    /**
     * \brief Computes the PCE stability flags based on the provided particle list and a number of parameters
     *
     * \param TPS                    Pointer to the particle list
     * \param SahaEquationForNuclei  Whether the Saha equation is applied to light nuclei
     * \param FreezeLongLived        Whether long-lived resonance yields are frozen at the chemical freeze-out
     * \param WidthCut               The threshold resonance width value to consider the resonance long-lived and its abundance frozen [GeV]
     */
    static std::vector<int> ComputePCEStabilityFlags(
      const ThermalParticleSystem* TPS,
      bool SahaEquationForNuclei = true,
      bool FreezeLongLived = false,
      double WidthCut = 0.015);

    /**
     * \brief Modifies the decay threshold masses of bosonic resonances such that the Bose-Condesation does not occur due to large fugacities.
     *
     *        This is only necessary if energy-dependent Breit-Wigner widths are used.
     */
    void ApplyFixForBoseCondensation();

    const std::vector< std::vector<double> >& EffectiveCharges() const { return m_EffectiveCharges; }
    const std::vector<int>& StableMapTo() const { return m_StableMapTo; }
    

  protected:
    ThermalModelBase *m_model;

    /// Whether nuclear abundances are calculated via the Saha equation
    bool m_UseSahaForNuclei;

    /// Whether long-lived resonances are frozen at Tch
    bool m_FreezeLonglivedResonances;
    /// Resonance width cut for freezeing the resonance abundances
    double m_ResoWidthCut;

    /// Whether the chemical freeze-out "initial" condition has been set
    bool m_ChemicalFreezeoutSet;

    /// Whether PCE has been calculated
    bool m_IsCalculated;

    /// PCE configuration, list of stable species etc.
    bool m_StabilityFlagsSet;
    std::vector<int> m_StabilityFlags;
    int m_StableComponentsNumber;
    std::vector< std::vector<double> > m_EffectiveCharges;
    std::vector<int> m_StableMapTo;

    /// Parameters at the chemical freeze-out
    ThermalModelParameters m_ParametersInit;
    std::vector<double> m_ChemInit;
    std::vector<double> m_DensitiesInit;
    std::vector<double> m_StableDensitiesInit;
    double m_EntropyDensityInit;
    double m_ParticleDensityInit;

    /// The current PCE thermal paratmeres and chemical potentials
    ThermalModelParameters m_ParametersCurrent;
    std::vector<double> m_ChemCurrent;

    class BroydenEquationsPCE : public BroydenEquations
    {
    public:
      BroydenEquationsPCE(ThermalModelPCE *model, int mode = 0) : BroydenEquations(), m_THM(model), m_Mode(mode) { m_N = m_THM->m_StableComponentsNumber + 1; }
      std::vector<double> Equations(const std::vector<double> &x);
    private:
      ThermalModelPCE *m_THM;
      int m_Mode;
    };

  };

} // namespace thermalfist

#endif

