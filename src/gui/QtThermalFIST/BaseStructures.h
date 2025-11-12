/*
 * Thermal-FIST package
 *
 * Copyright (c) 2014-2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef BASESTRUCTURES_H
#define BASESTRUCTURES_H
#include <string>

#include "HRGBase/ThermalModelBase.h"
#include "configinteractions.h"

/**
*   A structure containing the current thermal model configuration.
*/
struct ThermalModelConfig {
  enum ThermalModelType { Ideal = 0, DiagonalEV = 1, CrosstermsEV = 2, QvdW = 3, CE = 4, SCE = 5, CCE = 6, EVSCE = 7, VDWSCE = 8, RealGas = 9 };

  enum ThermalInteraction { InteractionIdeal = 0, InteractionEVDiagonal = 1, InteractionEVCrossterms = 2, InteractionQVDW = 3, InteractionRealGas = 4 };
  enum ThermalEnsemble { EnsembleGCE = 0, EnsembleCE = 1, EnsembleSCE = 2, EnsembleCCE = 3 };

  int ModelType;        /**< 0 - Ideal, 1 - Diagonal EV, 2 - Crossterms EV, 3 - QvdW, 4 - CE, 5 - SCE, 6 - Real Gas  */
  int Ensemble;
  int InteractionModel;
  int QuantumStatistics; /**< 0 - Boltzmann, 1 - Quantum */
  int QuantumStatisticsType; /**< 0 - Cluster Expansion, 1 - Quadratures */
  int QuantumStatisticsInclude; /**< 0 - All, 1 - Only mesons, 2 - Only pions */


  int InteractionScaling; /**< 0 - constant EV, 1 - bag model EV, 2 - two-component EV, 3 - from file, 4 - BB, B\bar{B}, MB, MM */
  //double EVRadius;
  double vdWA; /// QvdW attraction (in MeV * fm3)
  double vdWB; /// QvdW repulsion (in fm3)
  std::string InteractionInput;
  int DisableMM;
  int DisableMB;
  int DisableBB;
  int DisableBantiB;

  double vdWaBB, vdWaBantiB, vdWaMB, vdWaMM;
  double vdWbBB, vdWbBantiB, vdWbMB, vdWbMM;

  QvdWParameters vdWparams;

  int RealGasExcludedVolumePrescription; // 0 - vdW, 1 - CS, 2 - virial, 3 - TVM

  /// Thermal parameters
  double T;
  double muB, muQ, muS, muC;
  double gq, gS, gC;
  double VolumeR, VolumeRSC;

  int B, S, Q, C;
  bool CanonicalB;
  bool CanonicalQ;
  bool CanonicalS;
  bool CanonicalC;

  /// Constraints on mu's
  double SoverB;
  double QoverB;
  double RhoB;
  double RhoQ;
  double RhoS;
  double RhoC;
  bool ConstrainMuB;
  bool ConstrainMuQ;
  bool ConstrainMuS;
  bool ConstrainMuC;
  int ConstrainMuBType; // 0 - Entropy per baryon, 1 - Baryon density

  /// Extra flags
  int FiniteWidth; /**< 0 - zero, 1 - BW-2Gamma, 2 - eBW */
  int WidthShape;  /**< 0 - Relativistic Breit-Wigner, 1 - Nonrelativistic Breit-Wigner */
  bool RenormalizeBR;
  bool ComputeFluctations;
  bool ResetMus;


  /// Partial chemical equilibrium
  bool UsePCE;
  double Tkin;
  bool PCEFreezeLongLived;
  double PCEWidthCut;
  bool PCESahaForNuclei;
  bool PCEAnnihilation;
  double PCEPionAnnihilationNumber;

  /// Whether to use rejection sampling instead of importance sampling for the EV multiplicity sampling
  bool fUseEVRejectionMultiplicity;

  /// Whether to use rejection sampling in the coordinate space to model EV effects
  bool fUseEVRejectionCoordinates;

  /// Whether to use the SPR (single-particle rejection) approximation for the EV effects in coordinate space
  bool fUseEVUseSPRApproximation;

  /// Effective mass model (pions)
  bool  UseEMMPions;
  double EMMPionFPi;

  /// Effective mass model (kaons)
  bool  UseEMMKaons;
  double EMMKaonFKa;

  double MagneticFieldB;
  int MagneticFieldLmax;

  static ThermalModelConfig fromThermalModel(thermalfist::ThermalModelBase *model);
};

struct ThermodynamicsCosmic {
  bool flag;
  double T, muB, muQ, mu_e, mu_mu, mu_tau;
  double pT4, eT4, IT4, sT3;
  double p, e, I, s;
  double rhoB, rhoQ, rhoE, rhoMu, rhoTau;
  std::vector<double> densities; // gamma, e+-, mu+-, tau+-, neutrinos
};

struct Thermodynamics {
  bool flag;
  double T, muB, muQ, muS, muC;
  double pT4, eT4, IT4, sT3;
  double p, e, I, s;
  double rhoB, rhoQ, rhoS, rhoC;
  double nh;
  double cs2, cVT3;
  std::vector< std::vector<double> > densities;
};

struct ChargesFluctuations {
  bool flag;
  double chi1B, chi2B, chi3B, chi4B;
  double chi1Q, chi2Q, chi3Q, chi4Q;
  double chi1S, chi2S, chi3S, chi4S;
  double chi1C, chi2C, chi3C, chi4C;
  double chi11BQ, chi11BS, chi11QS;
  double chi31BQ, chi31BS, chi31QS;
};

void SetThermalModelParameters(thermalfist::ThermalModelBase* model, const ThermalModelConfig& config);
void SetThermalModelConfiguration(thermalfist::ThermalModelBase *model, const ThermalModelConfig &config);
void SetThermalModelInteraction(thermalfist::ThermalModelBase *model, const ThermalModelConfig &config);

#endif