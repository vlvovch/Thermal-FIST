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
	int FiniteWidth; /**< 0 - zero, 1 - BW */
	bool RenormalizeBR;
	bool ComputeFluctations;
};

#endif