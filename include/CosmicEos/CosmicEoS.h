#ifndef COSMICEOS_H
#define COSMICEOS_H
#include <map>

#include "HRGBase/ThermalModelBase.h"
#include "HRGBase/Broyden.h"
#include "CosmicEos/EffectiveMassModel.h"

namespace thermalfist {

  // Isospin = (u - d) / 2, e.g. I_pi = 1,0,-1
  // Using https://en.wikipedia.org/wiki/Gell-Mann%E2%80%93Nishijima_formula
  double IsospinCharge(const ThermalParticle& part);

  /**
   * \brief An auxiliary struct containing the list of conserved lepton flavor charges.
   *
   */
  struct LeptonFlavor {
    /**
     * \brief Set of all conserved charges considered.
     *
     */
    enum Name {
      Electron = 0, ///< Electron
      Muon = 1,     ///< Muon
      Tau = 2,      ///< Tau
    };

    /// Number of lepton flavors
    static const int NumberOfFlavors = 3;

    /// Electron mass
    static double m_e;
    
    /// Muon mass
    static double m_mu;

    /// Tauon mass
    static double m_tau;

  };

  /**
   * \brief Class implementing cosmological equation of state
   *        as a mixture of HRG model, ideal gases of leptons and photons
   *
   */
  class CosmicEoS
  {
  public:

    /**
      * \brief Constructor.
      * 
      * \param THMbase        Pointer to the HRG model object.
      * \param pionsinteract  Specifies whether to incorporate pions interactions through effective mass model.
      */
    CosmicEoS(ThermalModelBase* THMbase, bool pionsinteract = false);

    /// Destructor
    virtual ~CosmicEoS(void) { }

    /// Pointer to the HRG model object
    ThermalModelBase* HRGModel() const { return m_modelHRG; }

    /// Chemical potentials (baryon number, electric charge, three lepton charges)
    const std::vector<double>& ChemicalPotentials() const { return m_ChemCurrent; }

    /**
     * \brief Set the temperature
     *
     * \param T Temperature (GeV)
     */
    virtual void SetTemperature(double T);
    double Temperature() const { return m_T; }

    /**
     * \brief Set the baryon chemical potential
     *
     * \param muB Baryon chemical potential (GeV)
     */
    virtual void SetBaryonChemicalPotential(double muB) { m_ChemCurrent[0] = muB; m_modelHRG->SetBaryonChemicalPotential(muB); }
    double BaryonChemicalPotential() const { return m_ChemCurrent[0]; }

    /**
     * \brief Set the electric chemical potential
     *
     * \param muQ Electric chemical potential (GeV)
     */
    virtual void SetElectricChemicalPotential(double muQ);// { m_ChemCurrent[1] = muQ; m_modelHRG->SetElectricChemicalPotential(muQ); }
    double ElectricChemicalPotential() const { return m_ChemCurrent[1]; }

    /**
     * \brief Set the lepton chemical potential
     *
     * \param muL Lepton chemical potential (GeV)
     * \param flavor Lepton flavor
     */
    virtual void SetLeptonChemicalPotential(LeptonFlavor::Name flavor, double muL) { m_ChemCurrent[2 + static_cast<int>(flavor)] = muL; }
    double LeptonChemicalPotential(LeptonFlavor::Name flavor) const { return m_ChemCurrent[2 + static_cast<int>(flavor)]; }

    /// Set baryon number asymmetry (baryon over entropy density)
    virtual void SetBaryonAsymmetry(double b) { m_Asymmetries[0] = b; }
    double BaryonAsymmetry() const { return m_Asymmetries[0]; }

    /// Set electric charge asymmetry (charge over entropy density)
    virtual void SetChargeAsymmetry(double q) { m_Asymmetries[1] = q; }
    double ChargeAsymmetry() const { return m_Asymmetries[1]; }

    /// Set lepton flavor asymmetry (lepton flavor over entropy density)
    virtual void SetLeptonAsymmetry(LeptonFlavor::Name flavor, double l) { m_Asymmetries[2 + static_cast<int>(flavor)] = l; }
    double LeptonAsymmetry(LeptonFlavor::Name flavor) const { return m_Asymmetries[2 + static_cast<int>(flavor)]; }



    /// Set all asymmetries at once
    void SetAsymmetries(const std::vector<double>& asymmetries) { m_Asymmetries = asymmetries; }

    /// Calculate number densities of all particle species
    void CalculatePrimordialDensities();

    /// Total entropy density
    double EntropyDensity();

    /// Entropy density of the HRG part
    double EntropyDensityHRG();

    /// Total pressure
    double Pressure();

    /// Partial pressure of the HRG part
    double PressureHRG();

    /// Total energy density
    double EnergyDensity();

    /// Energy density of the HRG part
    double EnergyDensityHRG();

    /// Net density of charged lepton flavor
    double NetDensityChargedLepton(int iL);

    /// Partial pressure of charged lepton flavor 
    double PressureChargedLepton(int iL);

    /// Energy density of charged lepton flavor
    double EnergyDensityChargedLepton(int iL);
    
    /// \brief Total baryon density
    /// \param absolute If true, calculates absolute baryon density
    double BaryonDensity(bool absolute = false);

    /// \brief Electric charge density
    /// \param absolute If true, calculates absolute charge density
    double ElectricChargeDensity(bool absolute = false);

    /// \brief Electric charge density of the HRG part
    /// \param absolute If true, calculates absolute charge density
    double ElectricChargeDensityHRG(bool absolute = false);

    /// \brief Isospin charge density
    /// \param absolute If true, calculates absolute isospin density
    double IsospinChargeDensity(bool absolute = false);

    /// \brief Lepton flavor density of the HRG part
    /// \param flavor   Lepton flavor
    /// \param absolute If true, calculates absolute lepton flavor density
    double LeptonFlavorDensity(LeptonFlavor::Name flavor, bool absolute = false);

    /// Calculates the values of the chemical potential (B,Q,{L})
    /// that satisfy the given asymmetry constraints at given temperature
    ///
    /// \param T      Temperature [GeV]
    /// \param muInit Initial guesses of the chemical potentials
    std::vector<double> SolveChemicalPotentials(double T, const std::vector<double>& muInit = std::vector<double>());

    // Pion decay constant for pion interactions a la ChPT
    static double fpi;

    /// Mass of pi+
    /// Used to determine if Bose condensation reached in ideal gas
    double GetPionMass() const;

    /// Whether to include pion interactions via effective mass model
    void SetPionsInteracting(bool pionsinteract = true, double fpiChPT = fpi);
    bool InteractingPions() const { return m_InteractingPions; }

    /// Whether the system has non-zero BEC of pions (as a result of the calculation)
    bool InPionCondensedPhase() const;

    // std::vector<EffectiveMassModel>& EMMPions() { return m_Pions; }
    // double EMMPionCharge(int ipi) const { return m_PionCharges[ipi]; }

    /// ThermalParticle() object instance corresponding to photons
    const ThermalParticle& PhotonParticle() const { return m_Photon; }

    int NumberOfElectroWeakSpecies() const { return 1 + 2 * m_ChargedLeptons.size() + 2 * m_Neutrinos.size(); }

    /// Returns the name of particle species of given id
    std::string GetSpeciesName(int id) const;

    /// @brief Returns the number density for given species
    /// @param id 0 - photon, 1 - e+, 2 - e-,  3 - mu+, ..., 7 - nu_e, 8 - anti-nu_e, ...
    /// @return 
    double GetDensity(int id) const;

  protected:
    /// Pointer to an HRG model object
    ThermalModelBase* m_modelHRG;

    /// Photons
    ThermalParticle m_Photon;

    /// Charged leptons
    std::vector<ThermalParticle> m_ChargedLeptons;
    
    /// Neutrinos
    std::vector<ThermalParticle> m_Neutrinos;

    // std::vector<EffectiveMassModel> m_Pions;
    // std::vector<double> m_PionCharges;

    bool m_InteractingPions;

    bool m_IsCalculated;

    double m_T;
    std::vector<double> m_ChemCurrent;
    std::vector<double> m_Asymmetries;

    void ClearEMMs();

    class BroydenEquationsCosmology : public BroydenEquations
    {
    public:
      BroydenEquationsCosmology(CosmicEoS* model) : BroydenEquations(), m_THM(model) { m_N = 2 + LeptonFlavor::NumberOfFlavors; }
      std::vector<double> Equations(const std::vector<double>& x);
    private:
      CosmicEoS* m_THM;
    };
  };

} // namespace thermalfist

#endif

