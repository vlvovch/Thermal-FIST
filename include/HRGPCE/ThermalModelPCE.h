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
   * Implementation is still under development.
   * Use the current class on your own risk.
   *
   * Early ideas about PCE can be found in
   *
   * H. Bebie, P. Gerber, J.L. Goity, H. Leutwyler,
   * Nucl. Phys. B **378**, 95 (1992)
   * 
   * The present implementation was used (and described) in the paper
   * 
   * V. Vovchenko, K. Gallmeister, J. Schaffner-Bielich, C. Greiner,
   * arXiv:1904.10024,
   * [http://arxiv.org/pdf/1904.10024.pdf](http://arxiv.org/pdf/1904.10024.pdf) 
   * 
   */
  class ThermalModelPCE
  {
  public:

    ThermalModelPCE(ThermalModelBase *THMbase, bool FreezeLongLived = false);

    virtual ~ThermalModelPCE(void) { }

    virtual void SetStabilityFlags(const std::vector<int>& StabilityFlags);
    const std::vector<int>& StabilityFlags() const { return m_StabilityFlags; }

    void SetChemicalFreezeout(const ThermalModelParameters& params, const std::vector<double>& ChemInit);

    ThermalModelBase* ThermalModel() const { return m_model; }

    virtual void CalculatePCE(double T);

    const std::vector<double>& ChemicalPotentials() const { return m_ChemCurrent; }
    double Volume() const { return m_ParametersCurrent.V; }

    // Fills the "decay" products of light nuclei in the table in accordance with their baryon content
    /**
     * \brief Fills the "decay" products of light nuclei in accordance with their baryon content
     * 
     * \param TPS  Pointer to the particle list
     */
    static void PrepareNucleiForPCE(ThermalParticleSystem *TPS);

  protected:
    ThermalModelBase *m_model;
    bool m_ChemicalFreezeoutSet;
    bool m_IsCalculated;
    std::vector<int> m_StabilityFlags;
    ThermalModelParameters m_ParametersInit;
    std::vector<double> m_ChemInit;
    std::vector<double> m_DensitiesInit;
    double m_EntropyDensityInit;
    int m_StableComponentsNumber;
    std::vector< std::vector<double> > m_EffectiveCharges;
    std::vector<int> m_StableMapTo;
    //std::map<int> m_StableMapFrom;

    ThermalModelParameters m_ParametersCurrent;
    std::vector<double> m_ChemCurrent;

    class BroydenEquationsPCE : public BroydenEquations
    {
    public:
      BroydenEquationsPCE(ThermalModelPCE *model) : BroydenEquations(), m_THM(model) { m_N = m_THM->m_StableComponentsNumber + 1; }
      std::vector<double> Equations(const std::vector<double> &x);
    private:
      ThermalModelPCE *m_THM;
    };

  };

} // namespace thermalfist

#endif

