/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2022 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef THERMALMODELEVCROSSTERMS_H
#define THERMALMODELEVCROSSTERMS_H

#include "HRGVDW/ThermalModelVDW.h"

namespace thermalfist {

  /**
   * \brief Class implementing the crossterms
   *        excluded-volume model.
   * 
   * The model formulation can be found in
   * 
   * M.I. Gorenstein, A.P. Kostyuk, Y.D. Krivenko,
   * J.Phys. G **25**, L75 (1999),
   * [http://arxiv.org/pdf/nucl-th/9906068.pdf](http://arxiv.org/pdf/nucl-th/9906068.pdf)
   * 
   * and in
   * 
   * V. Vovchenko, H. Stoecker,
   * Phys. Rev. C **95**, 044904 (2017),
   * [http://arxiv.org/pdf/1606.06218.pdf](http://arxiv.org/pdf/1606.06218.pdf) 
   * 
   * This class implements the crossterms excluded-volume model as a partial case of
   * the van der Waals HRG model where all the attraction terms are set to zero.
   * 
   */
  class ThermalModelEVCrossterms : public ThermalModelVDW
  {
  public:
    /**
     * \brief Construct a new ThermalModelEVCrossterms object
     * 
     * \param TPS A pointer to the ThermalParticleSystem object containing the particle list
     * \param params ThermalModelParameters object with current thermal parameters
     */
    ThermalModelEVCrossterms(ThermalParticleSystem *TPS, const ThermalModelParameters& params = ThermalModelParameters()) :
      ThermalModelVDW(TPS, params) {
      m_TAG = "ThermalModelEVCrossterms";
      m_InteractionModel = CrosstermsEV;
    }

    /**
     * \brief Destroy the ThermalModelEVCrossterms object
     * 
     */
    virtual ~ThermalModelEVCrossterms(void) { }

    // Override functions begin

    virtual void ReadInteractionParameters(const std::string &filename);

    virtual void WriteInteractionParameters(const std::string &filename);

    virtual void SetRadius(double rad) { FillVirial(std::vector<double>(m_TPS->Particles().size(), rad)); }

    virtual void SetAttraction(int i, int j, double a);

    /// No need to search for multiple soultions in EV-HRG model
    virtual void SetMultipleSolutionsMode(bool search);

    // Override functions end

    const std::vector< std::vector<int> >& EVComponentIndices() const { return VDWComponentIndices(); }
  };

} // namespace thermalfist

#endif

