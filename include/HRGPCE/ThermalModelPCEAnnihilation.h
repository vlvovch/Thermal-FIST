#ifndef THERMALMODELPCEANNIHILATION_H
#define THERMALMODELPCEANNIHILATION_H
#include <map>

#include "HRGPCE/ThermalModelPCE.h"
#include "HRGBase/Broyden.h"

namespace thermalfist {

  /**
  * \brief Class implementing HRG in partial chemical equilibrium with baryon annihilation.
  *
  * Implements equilibrium of N\bar{N} <-> <N_pi} pi reactions
  * N is protons and neutrons by default, but can also include other stable baryons like Lambda's
  *
  * The implementation was used (and described) in a paper
  * V. Vovchenko, V. Koch,
  * Phys. Lett. B **835**, 137577 (2022),
  * [arXiv:2210.15641](https://arxiv.org/abs/2210.15641)
  *
  * A typical usage pattern is similar to ThermalModelPCE() except
  * one can specify the annihilating (anti)baryons by a value 2 through SetStabilityFlags()
  *
  */
  class ThermalModelPCEAnnihilation : public ThermalModelPCE
  {
  public:

    /**
     * \brief Construct a new ThermalModelPCEAnnihilation object
     *
     * \param THMbase               A pointer to the ThermalModelBase object containing the HRG model implementation
     * \param FreezeLonglived       Whether long-lived resonance abundances should be frozen at Tch along with the stable hadrons
     * \param LonglivedResoWidthCut Threshold width for the long-lived resonances to be considered stable in the PCE
     */
    ThermalModelPCEAnnihilation(ThermalModelBase* THMbase, bool FreezeLonglived = false, double LonglivedResoWidthCut = 0.015);

    /**
     * \brief Destroy the ThermalModelPCEAnnihilation object
     *
     */
    virtual ~ThermalModelPCEAnnihilation(void) { }

    std::vector<int> RecalculateStabilityFlags(const std::vector<long long>& annihilationpdgs = {2212, 2112});

    /**
     * \brief Manually set the PCE stability flags for all species.
     *
     * \param StabilityFlags Vector of stability flags.
     *                       Each element corresponds to a particle specie with the same 0-based index in the particle list.
     *                       Element i equal to zero means particle yield is not frozen in the hadronic phase,
     *                       if it is equal to one it means its yield is stable,
     *                       and if it is equal to two it means the corresponding baryon annihilates.
     */
    virtual void SetStabilityFlags(const std::vector<int>& StabilityFlags);

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
     * \brief Set the average number of pion produced in baryon-antibaryon annihilations.
     *
     * \param npi  The average number of pion produced in baryon-antibaryon annihilations.
     */
    void SetPionAnnihilationNumber(double npi) { m_PionAnnihilationNumber = npi; }

  protected:

    /// Returns the PCE-based index of the stable hadron based on its global (particle list) index
    int StableHadronIndexByGlobalId(int globalid);

    /// Chemical potentials of all PCE-based hadrons from the solution to PCE equations
    std::vector<double> StableChemsFromBroydenInput(const std::vector<double>& x);

  private:

    std::vector<int> m_StableNormal;
    std::vector<int> m_StableAnnihilate;
    std::vector<int> m_Pions;
    std::vector<int> m_StableAnnihilateAnti;

    double m_PionAnnihilationNumber;

    class BroydenEquationsPCEAnnihilation : public BroydenEquations
    {
    public:
      BroydenEquationsPCEAnnihilation(ThermalModelPCEAnnihilation *model, PCEMode mode = PCEMode::AtFixedTemperature) : BroydenEquations(), m_THM(model), m_Mode(mode) { m_N = m_THM->m_StableNormal.size() + m_THM->m_StableAnnihilate.size() + m_THM->m_Pions.size() + 1; }
      std::vector<double> Equations(const std::vector<double> &x);
    private:
      ThermalModelPCEAnnihilation *m_THM;
      PCEMode m_Mode;
    };

  };

} // namespace thermalfist

#endif

