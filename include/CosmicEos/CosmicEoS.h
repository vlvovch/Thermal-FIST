#ifndef COSMICEOS_H
#define COSMICEOS_H
#include <map>

#include "HRGBase/ThermalModelBase.h"
#include "HRGBase/Broyden.h"
#include "CosmicEos/EffectiveMassModel.h"

namespace thermalfist {

  /**
   * \brief Calculates the isospin charge of a particle.
   * 
   * Isospin = (u - d) / 2, e.g. I_pi = 1,0,-1
   * Using https://en.wikipedia.org/wiki/Gell-Mann%E2%80%93Nishijima_formula
   * 
   * \param part The particle for which to calculate the isospin charge
   * \return double The isospin charge value
   */
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

    /**
     * \brief Destructor.
     */
    virtual ~CosmicEoS(void) { }

    /**
     * \brief Gets the pointer to the HRG model object.
     * 
     * \return Pointer to the HRG model object.
     */
    ThermalModelBase* HRGModel() const { return m_modelHRG; }

    /**
     * \brief Gets the chemical potentials (baryon number, electric charge, three lepton charges).
     * 
     * \return Vector of chemical potentials.
     */
    const std::vector<double>& ChemicalPotentials() const { return m_ChemCurrent; }

    /**
     * \brief Set the temperature
     *
     * \param T Temperature (GeV)
     */
    virtual void SetTemperature(double T);
    
    /**
     * \brief Gets the current temperature.
     * 
     * \return The temperature value in GeV.
     */
    double Temperature() const { return m_T; }

    /**
     * \brief Set the baryon chemical potential
     *
     * \param muB Baryon chemical potential (GeV)
     */
    virtual void SetBaryonChemicalPotential(double muB) { m_ChemCurrent[0] = muB; m_modelHRG->SetBaryonChemicalPotential(muB); }
    
    /**
     * \brief Gets the baryon chemical potential.
     * 
     * \return The baryon chemical potential value in GeV.
     */
    double BaryonChemicalPotential() const { return m_ChemCurrent[0]; }

    /**
     * \brief Set the electric chemical potential
     *
     * \param muQ Electric chemical potential (GeV)
     */
    virtual void SetElectricChemicalPotential(double muQ);// { m_ChemCurrent[1] = muQ; m_modelHRG->SetElectricChemicalPotential(muQ); }
    
    /**
     * \brief Gets the current electric chemical potential.
     * 
     * \return The electric chemical potential value in GeV.
     */
    double ElectricChemicalPotential() const { return m_ChemCurrent[1]; }

    /**
     * \brief Set the lepton chemical potential
     *
     * \param muL Lepton chemical potential (GeV)
     * \param flavor Lepton flavor
     */
    virtual void SetLeptonChemicalPotential(LeptonFlavor::Name flavor, double muL) { m_ChemCurrent[2 + static_cast<int>(flavor)] = muL; }
    
    /**
     * \brief Gets the current lepton chemical potential for a specific flavor.
     * 
     * \param flavor The lepton flavor.
     * \return The lepton chemical potential value in GeV.
     */
    double LeptonChemicalPotential(LeptonFlavor::Name flavor) const { return m_ChemCurrent[2 + static_cast<int>(flavor)]; }

    /**
     * \brief Sets the baryon number asymmetry (baryon over entropy density).
     * 
     * \param b The baryon asymmetry value.
     */
    virtual void SetBaryonAsymmetry(double b) { m_Asymmetries[0] = b; }
    
    /**
     * \brief Gets the current baryon asymmetry.
     * 
     * \return The baryon asymmetry value.
     */
    double BaryonAsymmetry() const { return m_Asymmetries[0]; }

    /**
     * \brief Sets the electric charge asymmetry (charge over entropy density).
     * 
     * \param q The charge asymmetry value.
     */
    virtual void SetChargeAsymmetry(double q) { m_Asymmetries[1] = q; }
    
    /**
     * \brief Gets the current electric charge asymmetry.
     * 
     * \return The charge asymmetry value.
     */
    double ChargeAsymmetry() const { return m_Asymmetries[1]; }

    /**
     * \brief Sets the lepton flavor asymmetry (lepton flavor over entropy density).
     * 
     * \param flavor The lepton flavor.
     * \param l The lepton asymmetry value.
     */
    virtual void SetLeptonAsymmetry(LeptonFlavor::Name flavor, double l) { m_Asymmetries[2 + static_cast<int>(flavor)] = l; }
    
    /**
     * \brief Gets the current lepton flavor asymmetry.
     * 
     * \param flavor The lepton flavor.
     * \return The lepton asymmetry value.
     */
    double LeptonAsymmetry(LeptonFlavor::Name flavor) const { return m_Asymmetries[2 + static_cast<int>(flavor)]; }

    /**
     * \brief Sets all asymmetries at once.
     * 
     * \param asymmetries Vector of asymmetries (baryon, charge, lepton flavors).
     */
    void SetAsymmetries(const std::vector<double>& asymmetries) { m_Asymmetries = asymmetries; }

    /**
     * \brief Calculates number densities of all particle species.
     */
    void CalculatePrimordialDensities();

    /**
     * \brief Calculates the total entropy density.
     * 
     * \return The total entropy density value.
     */
    double EntropyDensity();

    /**
     * \brief Calculates the entropy density of the HRG part.
     * 
     * \return The HRG entropy density value.
     */
    double EntropyDensityHRG();

    /**
     * \brief Calculates the total pressure.
     * 
     * \return The total pressure value.
     */
    double Pressure();

    /**
     * \brief Calculates the partial pressure of the HRG part.
     * 
     * \return The HRG pressure value.
     */
    double PressureHRG();

    /**
     * \brief Calculates the total energy density.
     * 
     * \return The total energy density value.
     */
    double EnergyDensity();

    /**
     * \brief Calculates the energy density of the HRG part.
     * 
     * \return The HRG energy density value.
     */
    double EnergyDensityHRG();

    /**
     * \brief Calculates the net density of charged lepton flavor.
     * 
     * \param iL The lepton flavor index.
     * \return The net density value.
     */
    double NetDensityChargedLepton(int iL);

    /**
     * \brief Calculates the partial pressure of charged lepton flavor.
     * 
     * \param iL The lepton flavor index.
     * \return The partial pressure value.
     */
    double PressureChargedLepton(int iL);

    /**
     * \brief Calculates the energy density of charged lepton flavor.
     * 
     * \param iL The lepton flavor index.
     * \return The energy density value.
     */
    double EnergyDensityChargedLepton(int iL);
    
    /**
     * \brief Calculates the total baryon density.
     * 
     * \param absolute If true, calculates absolute baryon density.
     * \return The baryon density value.
     */
    double BaryonDensity(bool absolute = false);

    /**
     * \brief Calculates the electric charge density.
     * 
     * \param absolute If true, calculates absolute charge density.
     * \return The electric charge density value.
     */
    double ElectricChargeDensity(bool absolute = false);

    /**
     * \brief Calculates the electric charge density of the HRG part.
     * 
     * \param absolute If true, calculates absolute charge density.
     * \return The HRG electric charge density value.
     */
    double ElectricChargeDensityHRG(bool absolute = false);

    /**
     * \brief Calculates the isospin charge density.
     * 
     * \param absolute If true, calculates absolute isospin density.
     * \return The isospin charge density value.
     */
    double IsospinChargeDensity(bool absolute = false);

    /**
     * \brief Calculates the lepton flavor density.
     * 
     * \param flavor The lepton flavor.
     * \param absolute If true, calculates absolute lepton flavor density.
     * \return The lepton flavor density value.
     */
    double LeptonFlavorDensity(LeptonFlavor::Name flavor, bool absolute = false);

    /**
     * \brief Calculates the values of the chemical potential (B,Q,{L})
     *        that satisfy the given asymmetry constraints at given temperature.
     *
     * \param T Temperature in GeV.
     * \param muInit Initial guesses of the chemical potentials.
     * \return Vector of calculated chemical potentials.
     */
    std::vector<double> SolveChemicalPotentials(double T, const std::vector<double>& muInit = std::vector<double>());

    /**
     * \brief Pion decay constant for pion interactions a la ChPT.
     */
    static double fpi;

    /**
     * \brief Gets the mass of pi+.
     * 
     * Used to determine if Bose condensation reached in ideal gas.
     * 
     * \return The mass of pi+ in GeV.
     */
    double GetPionMass() const;

    /**
     * \brief Sets whether to include pion interactions via effective mass model.
     * 
     * \param pionsinteract If true, pion interactions are included.
     * \param fpiChPT The pion decay constant (default: fpi).
     */
    void SetPionsInteracting(bool pionsinteract = true, double fpiChPT = fpi);
    
    /**
     * \brief Checks if pions are interacting in the model.
     * 
     * \return True if pions are interacting, false otherwise.
     */
    bool InteractingPions() const { return m_InteractingPions; }

    /**
     * \brief Checks if the system has non-zero BEC of pions.
     * 
     * \return True if the system is in a pion condensed phase, false otherwise.
     */
    bool InPionCondensedPhase() const;

    // std::vector<EffectiveMassModel>& EMMPions() { return m_Pions; }
    // double EMMPionCharge(int ipi) const { return m_PionCharges[ipi]; }

    /**
     * \brief Gets the ThermalParticle object instance corresponding to photons.
     * 
     * \return Reference to the photon particle object.
     */
    const ThermalParticle& PhotonParticle() const { return m_Photon; }

    /**
     * \brief Gets the number of electroweak species (charged leptons and neutrinos) in the model.
     * 
     * \return The number of electroweak species.
     */
    int NumberOfElectroWeakSpecies() const { return 1 + 2 * m_ChargedLeptons.size() + 2 * m_Neutrinos.size(); }

    /**
     * \brief Gets the name of particle species of given id.
     * 
     * \param id The particle species id.
     * \return The name of the particle species.
     */
    std::string GetSpeciesName(int id) const;

    /**
     * \brief Gets the number density for given species.
     * 
     * \param id 0 - photon, 1 - e+, 2 - e-,  3 - mu+, ..., 7 - nu_e, 8 - anti-nu_e, ...
     * \return The number density value.
     */
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

    /// Whether to include pion interactions
    bool m_InteractingPions;

    /// Whether the system has been calculated
    bool m_IsCalculated;

    /// Temperature in GeV
    double m_T;
    
    /// Vector of chemical potentials (baryon, charge, lepton flavors)
    std::vector<double> m_ChemCurrent;
    
    /// Vector of asymmetries (baryon, charge, lepton flavors)
    std::vector<double> m_Asymmetries;

    /**
     * \brief Clears the effective mass models.
     */
    void ClearEMMs();

    /**
     * \brief Class implementing the Broyden equations for cosmology.
     * 
     * This class is used to solve for chemical potentials that satisfy
     * the asymmetry constraints.
     */
    class BroydenEquationsCosmology : public BroydenEquations
    {
    public:
      /**
       * \brief Constructor for the BroydenEquationsCosmology class.
       * 
       * \param model Pointer to the CosmicEoS model.
       */
      BroydenEquationsCosmology(CosmicEoS* model) : BroydenEquations(), m_THM(model) { m_N = 2 + LeptonFlavor::NumberOfFlavors; }
      
      /**
       * \brief Implements the equations to be solved.
       * 
       * \param x The vector of variables (chemical potentials).
       * \return The vector of equation values.
       */
      std::vector<double> Equations(const std::vector<double>& x);
    private:
      /// Pointer to the CosmicEoS model
      CosmicEoS* m_THM;
    };
  };

} // namespace thermalfist

#endif
