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
   *
   * H. Bebie, P. Gerber, J.L. Goity, H. Leutwyler,
   * Nucl. Phys. B **378**, 95 (1992)
   * 
   * The present implementation was used (and described) in a paper
   * 
   * V. Vovchenko, K. Gallmeister, J. Schaffner-Bielich, C. Greiner,
   * [arXiv:1903.10024](https://arxiv.org/abs/1903.10024)
   * 
   */
  class ThermalModelPCE
  {
  public:

    ThermalModelPCE(ThermalModelBase *THMbase, bool FreezeLongLived = false, double WidthCut = 0.015);

    virtual ~ThermalModelPCE(void) { }

    virtual void SetStabilityFlags(const std::vector<int>& StabilityFlags);
    const std::vector<int>& StabilityFlags() const { return m_StabilityFlags; }

    void SetChemicalFreezeout(const ThermalModelParameters& params, const std::vector<double>& ChemInit);
    double EntropyDensityInit() const { return m_EntropyDensityInit; }
    void SetEntropyDensityInit(double sinit) { m_EntropyDensityInit = sinit; }

    ThermalModelBase* ThermalModel() const { return m_model; }

    /**
     * \brief Solves the equation of partial chemical equilibrium at fixed temperature or volume
     *
     * \param param  Temperature (in GeV) if mode == 0, Volume (in fm^3) if mode == 1
     * \param mode   PCE at fixed temperature for mode == 0, at fixed voluime for mode == 1
     */
    virtual void CalculatePCE(double param, int mode = 0);

    const std::vector<double>& ChemicalPotentials() const { return m_ChemCurrent; }
    double Volume() const { return m_ParametersCurrent.V; }

    /**
     * \brief Fills the "decay" products of light nuclei in accordance with their baryon content
     * 
     * \param TPS  Pointer to the particle list
     */
    static void PrepareNucleiForPCE(ThermalParticleSystem *TPS);

    static std::vector<int> ComputePCEStabilityFlags(
      const ThermalParticleSystem* TPS,
      bool SahaEquationForNuclei = true,
      bool FreezeLongLived = false,
      double WidthCut = 0.015);

    void ApplyFixForBoseCondensation();
    

  protected:
    ThermalModelBase *m_model;
    bool m_ChemicalFreezeoutSet;
    bool m_IsCalculated;
    std::vector<int> m_StabilityFlags;
    ThermalModelParameters m_ParametersInit;
    std::vector<double> m_ChemInit;
    std::vector<double> m_DensitiesInit;
    std::vector<double> m_StableDensitiesInit;
    double m_EntropyDensityInit;
    double m_ParticleDensityInit;
    int m_StableComponentsNumber;
    std::vector< std::vector<double> > m_EffectiveCharges;
    std::vector<int> m_StableMapTo;
    //std::map<int> m_StableMapFrom;

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

