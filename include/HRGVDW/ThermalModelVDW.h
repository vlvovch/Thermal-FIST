/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef THERMALMODELVDW_H
#define THERMALMODELVDW_H

#include "HRGBase/ThermalModelBase.h"

namespace thermalfist {

  class ThermalModelVDW : public ThermalModelBase
  {
  public:
    ThermalModelVDW(ThermalParticleSystem *TPS_, const ThermalModelParameters& params = ThermalModelParameters(), int mode = 0);

    virtual ~ThermalModelVDW(void);

    void FillChemicalPotentials();
    virtual void SetChemicalPotentials(const std::vector<double> & chem = std::vector<double>(0));

    void FillVirial(const std::vector<double> & ri = std::vector<double>(0));
    void FillVirialEV(const std::vector< std::vector<double> > & bij = std::vector< std::vector<double> >(0));
    void FillAttraction(const std::vector< std::vector<double> > & aij = std::vector< std::vector<double> >(0));

    virtual void ReadInteractionParameters(const std::string &filename);
    virtual void WriteInteractionParameters(const std::string &filename);

    void SetVirial(int i, int j, double b) { if (i >= 0 && i < m_Virial.size() && j >= 0 && j < m_Virial[i].size()) m_Virial[i][j] = b; }
    void SetAttraction(int i, int j, double a) { if (i >= 0 && i < m_Attr.size() && j >= 0 && j < m_Attr[i].size())     m_Attr[i][j] = a; }

    double VirialCoefficient(int i, int j) const;
    double AttractionCoefficient(int i, int j) const;

    void SetVirialdT(int i, int j, double dbdT) { if (i >= 0 && i < m_VirialdT.size() && j >= 0 && j < m_VirialdT[i].size()) m_VirialdT[i][j] = dbdT; }
    void SetAttractiondT(int i, int j, double dadT) { if (i >= 0 && i < m_AttrdT.size() && j >= 0 && j < m_AttrdT[i].size())     m_AttrdT[i][j] = dadT; }

    double VirialCoefficientdT(int i, int j) const;
    double AttractionCoefficientdT(int i, int j) const;

    void SetTemperatureDependentAB(bool Tdep) { m_TemperatureDependentAB = Tdep; }
    bool TemperatureDependentAB() const { return m_TemperatureDependentAB; }

    virtual void SetParameters(double T, double muB, double muS, double muQ, double gammaS, double V, double R);

    virtual void SetParameters(const ThermalModelParameters& params);

    void UpdateParameters();

    virtual void ChangeTPS(ThermalParticleSystem *TPS_);

    std::vector<double> ComputeNp(const std::vector<double>& dmustar);
    std::vector<double> ComputeNp(const std::vector<double>& dmustar, const std::vector<double>& ns);

    virtual std::vector<double> SearchSolutionsSingle(const std::vector<double> & muStarInit);
    std::vector<double> SearchSolutionsMulti(int iters = 300);
    void SearchSolutions();

    virtual void CalculateDensities();

    virtual std::vector<double> CalculateChargeFluctuations(const std::vector<double> &chgs, int order = 4);
    virtual std::vector< std::vector<double> >  CalculateFluctuations(int order);


    void CalculateTwoParticleCorrelations();
    // TODO higher orders
    void CalculateFluctuations();


    virtual double CalculateEnergyDensity();

    virtual double CalculateEntropyDensity();

    // Dummy
    virtual double CalculateBaryonMatterEntropyDensity();

    virtual double CalculateMesonMatterEntropyDensity();

    virtual double CalculatePressure();

    virtual double ParticleScalarDensity(int part);

    void SetMultipleSolutionsMode(bool search) { m_SearchMultipleSolutions = search; }
    bool UseMultipleSolutionsMode() const { return m_SearchMultipleSolutions; }

    double GetMaxDiff() const { return m_MaxDiff; }
    bool   IsLastSolutionOK() const { return m_LastBroydenSuccessFlag && m_LastCalculationSuccessFlag; }

    double DensityId(int index) { return m_DensitiesId[index]; }

    double MuStar(int index) const { return m_MuStar[index]; }
    std::vector<double> GetMuStar() const { return m_MuStar; }
    void SetMuStar(const std::vector<double> & MuStar) { m_MuStar = MuStar; }

  protected:
    // TODO: test
    virtual double MuShift(int id);

    std::vector<double> m_DensitiesId;
    std::vector<double> m_scaldens;

    std::vector< std::vector<double> > m_Virial;
    std::vector< std::vector<double> > m_Attr;
    std::vector< std::vector<double> > m_VirialdT;
    std::vector< std::vector<double> > m_AttrdT;

    std::vector< std::vector<double> > m_chi;
    std::vector<double> m_chiarb;

    bool   m_TemperatureDependentAB;


    virtual void CalculateDensitiesOld();
    virtual void CalculateDensitiesNew();

  protected:
    bool   m_SearchMultipleSolutions;
    bool   m_LastBroydenSuccessFlag;
    double m_MaxDiff;

    std::vector<double> m_MuStar;

    std::vector<int> m_MapTodMuStar;
    std::vector<int> m_MapFromdMuStar;
    std::vector< std::vector<int> > m_dMuStarIndices;

  private:
      class BroydenEquationsVDW : public BroydenEquations
      {
      public:
        BroydenEquationsVDW(ThermalModelVDW *model) : BroydenEquations(), m_THM(model) { m_N = model->m_MapFromdMuStar.size(); }
        std::vector<double> Equations(const std::vector<double> &x);
      private:
        ThermalModelVDW *m_THM;
      };

      class BroydenJacobianVDW : public BroydenJacobian
      {
      public:
        BroydenJacobianVDW(ThermalModelVDW *model) : BroydenJacobian(), m_THM(model) { }
        Eigen::MatrixXd Jacobian(const std::vector<double> &x);
      private:
        ThermalModelVDW *m_THM;
      };

      class BroydenSolutionCriteriumVDW : public Broyden::BroydenSolutionCriterium
      {
      public:
        BroydenSolutionCriteriumVDW(ThermalModelVDW *model, double relative_error = Broyden::TOL) : Broyden::BroydenSolutionCriterium(relative_error), m_THM(model) { }
        virtual bool IsSolved(const std::vector<double>& x, const std::vector<double>& f, const std::vector<double>& xdelta = std::vector<double>()) const;
      protected:
        ThermalModelVDW *m_THM;
      };
  };

  /// For backward compatibility
  typedef ThermalModelVDW ThermalModelVDWFull;

} // namespace thermalfist

#endif

