/*
 * Thermal-FIST package
 *
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef BASESTRUCTURES_H
#define BASESTRUCTURES_H
#include <string>

/**
*   A structure containing the current thermal model configuration.
*/
struct ThermalModelConfig {
	enum ThermalModelType { Ideal = 0, DiagonalEV = 1, CrosstermsEV = 2, QvdW = 3, CE = 4, SCE = 5, CCE = 6, EVSCE = 7, VDWSCE = 8 };

	enum ThermalInteraction { InteractionIdeal = 0, InteractionEVDiagonal = 1, InteractionEVCrossterms = 2, InteractionQVDW = 3 };
	enum ThermalEnsemble { EnsembleGCE = 0, EnsembleCE = 1, EnsembleSCE = 2, EnsembleCCE = 3 };

	int ModelType;				/**< 0 - Ideal, 1 - Diagonal EV, 2 - Crossterms EV, 3 - QvdW, 4 - CE, 5 - SCE  */
	int Ensemble;
	int InteractionModel;
	int QuantumStatistics; /**< 0 - Boltzmann, 1 - Quantum */
	int QuantumStatisticsType; /**< 0 - Cluster Expansion, 1 - Quadratures */
	int QuantumStatisticsInclude; /**< 0 - All, 1 - Only mesons, 2 - Only pions */
	int Interaction; /**< 0 - constant EV, 1 - bag model EV, 2 - two-component EV, 3 - from file */
	double EVRadius;
	std::string InteractionInput;

	/// Thermal parameters
	double T;
	double muB, muQ, muS, muC;
	double gq, gS, gC;
	double VolumeR, VolumeRSC;
	int B, S, Q, C;

	/// Constraints on mu's
	double QoverB;
	bool ConstrainMuQ;
	bool ConstrainMuS;
	bool ConstrainMuC;

	/// Extra flags
	int FiniteWidth; /**< 0 - zero, 1 - BW-2Gamma, 2 - eBW */
	bool RenormalizeBR;
	bool ComputeFluctations;
};

struct Thermodynamics {
  bool flag;
  double T, muB, muQ, muS, muC;
  double pT4, eT4, IT4, sT3;
  double nhT3;
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

#endif