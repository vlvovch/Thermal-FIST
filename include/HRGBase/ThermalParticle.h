#ifndef THERMALPARTICLE_H
#define THERMALPARTICLE_H

/** \file ThermalParticle.h
*   Contains classes and methods for basic information and calculations
*   involving a particle in ideal gas GCE.
*   Also contains information about decays.
*/

#include <string>
#include <vector>
#include <cmath>

#include "HRGBase/xMath.h"
#include "HRGBase/IdealGasFunctions.h"

// These are used for canonical ensemble calculations
// TODO: Determine dynamically, from the particle list
// Update: Obsolete
//const int BMAX = 1;								/**< Maximum baryon charge */
//const int QMAX = 2;								/**< Maximum electric charge */
//const int SMAX = 3;								/**< Maximum strangeness charge */
//const int CMAX = 1;								/**< Maximum charm charge */

/**
*   A structure containing all thermal parameters of the model.
*   Used in ThermalParticle class when calculating the ideal gas functions.
*/
struct ThermalModelParameters {
	double		T;				/**< Temperature [GeV] */
	double		muB;			/**< Baryon chemical potential [GeV] */
	double		muS;			/**< Strangeness chemical potential [GeV] */
	double		muQ;			/**< Electric charge chemical potential [GeV] */
	double		muC;			/**< Charm chemical potential [GeV] */
	double		gammaq;	/**< Chemical non-equilibirum fugacity of light quarks */
	double		gammaS;	/**< Chemical non-equilibirum fugacity of strange quarks */
	double		gammaC;	/**< Chemical non-equilibirum fugacity of charm quarks */
	double		V;				/**< Total system volume [fm^3] */
	double  SVc;			/**< Canonical correlation volume [fm^3] */
	int			B;				/**< Total baryon charge (CE) */
	int			Q;				/**< Total electric charge (CE) */
	int			S;				/**< Total strangeness charge (CE) */
	int			C;				/**< Total charm charge (CE) */

	ThermalModelParameters(double pT=0.1, double pmuB=0.5, double pmuS=0., double pmuQ=0., double pgammaS=1., double pV = 4000., double pSVc=1., int pB = 2, int pQ = 2, int pS = 0, int pC = 0) :
		T(pT), muB(pmuB), muS(pmuS), muQ(pmuQ), muC(0.), gammaq(1.), gammaS(pgammaS), gammaC(1.), V(pV), SVc(pSVc), B(pB), Q(pQ), S(pS), C(pC) {
	}

	ThermalModelParameters(double pT, double pgammaS, double pV, int pB, int pQ, int pS, int pC = 0):
		T(pT), gammaq(1.), gammaS(pgammaS), gammaC(1.), V(pV), SVc(pV), B(pB), Q(pQ), S(pS), C(pC) {
	}
};

/**
*   A structure containing information about a single decay channel of a particle.
*/
struct ParticleDecay
{
	double mBratio;								/**< Branching ratio */
	std::vector<int> mDaughters;		/**< PDGID numbers of daughter particles */
	double mM0;										/**< Sum of masses of decay products */

	double mPole;
	double mL;                      /**< Orbital angular momentum for decay, used in the eBW scheme */
	std::vector<double> mBratioVsM; /**< Mass-dependent branching ratios */
	double mBratioAverage;          /**< Average branching ratios after integrating with the Boltzmann factor*/

	/// Constructor.
	/**
	*   Takes branching ratio and vector of PDG IDs of daughter particles.
	*/
	ParticleDecay(double bratio=0., const std::vector<int> &daughters=std::vector<int>(0)) : 
			mBratio(bratio), mDaughters(daughters), mM0(0.), mPole(0.), mL(0.),
		  mBratioVsM(std::vector<double>(0)), mBratioAverage(bratio) {
	}

	double ModifiedWidth(double m) const;
};

/**
*   A class containing all information about a particle.
*   Also contains implementation of calculation of various quantities in ideal GCE gas.
*/
class ThermalParticle
{
	public:
		enum ResonanceWidthShape					{ RelativisticBreitWiger, NonRelativisticBreitWiger};
		enum ResonanceWidthIntegration		{ ZeroWidth, BWTwoGamma, FullInterval, FullIntervalWeighted, eBW };

		ThermalParticle(bool Stable_=true, std::string Name="hadron", int PDGID=0, double Deg=1., int Stat=0, double Mass=0.,
			int Strange=0, int Baryon=0, int Charge=0, double AbsS=0., double Width=0., double Threshold=0., int Charm = 0, double AbsC = 0., double radius = 0.5, int Quark = 0);
		~ThermalParticle(void);


		void FillCoefficients();

		// Fill coefficients for mass integration in the eBW scheme
		void FillCoefficientsDynamical();

		// Total width (eBW scheme) at a given mass
		double TotalWidtheBW(double M) const;

		// Energy-dependent branching ratios
		std::vector<double> BranchingRatiosM(double M, bool eBW = true) const;

		// Thermal mass distribution (not normalized!)
		double ThermalMassDistribution(double M, double T, double Mu, double width);
		double ThermalMassDistribution(double M, double T, double Mu);

		void NormalizeBranchingRatios();
		void RestoreBranchingRatios();// { m_Decays = m_DecaysOrig; }
		
		double Density(const ThermalModelParameters &params, IdealGasFunctions::Quantity type = IdealGasFunctions::ParticleDensity, bool useWidth = 0, double pMu = 0., double dMu = 0.) const;

		double DensityCluster(int n, const ThermalModelParameters &params, IdealGasFunctions::Quantity type = IdealGasFunctions::ParticleDensity, bool useWidth = 0, double pMu = 0., double dMu = 0.) const;

		double chi(int index, const ThermalModelParameters &params, bool useWidth=0, double pMu = 0., double dMu = 0.) const; 

		double ScaledVariance(const ThermalModelParameters &params, bool useWidth=0, double pMu = 0., double dMu = 0.) const;
		double Skewness(const ThermalModelParameters &params, bool useWidth=0, double pMu = 0., double dMu = 0.) const;
		double Kurtosis(const ThermalModelParameters &params, bool useWidth=0, double pMu = 0., double dMu = 0.) const;

		double FD(double k, double T, double mu, double m) const;

		double GetAbsQ() const;

		double GetCharge(int index) const;
		double GetAbsCharge(int index) const;

		bool IsNeutral() const;

		bool IsStable() const { return m_Stable; }
		void SetStable(bool stable = true) { m_Stable = stable; }

		bool IsAntiParticle() const { return m_AntiParticle; }
		void SetAntiParticle(bool antpar = true) { m_AntiParticle = antpar; }

		const std::string& Name() const { return m_Name; }
		void SetName(std::string &name) { m_Name = name; }

		int  PdgId() const { return m_PDGID; }
		void SetPdgId(int PdgId) { m_PDGID = PdgId; }

		int  Degeneracy() const { return m_Degeneracy; }
		void SetDegeneracy(double deg) { m_Degeneracy = deg; }

		int  Statistics() const { return m_Statistics; }
		void UseStatistics(bool enable);

		double Mass() const { return m_Mass; }
		void SetMass(double mass);// { m_Mass = mass; }

		int BaryonCharge() const { return m_Baryon; }
		void SetBaryonCharge(int chg) { m_Baryon = chg; }

		int ElectricCharge() const { return m_ElectricCharge; }
		void SetElectricCharge(int chg) { m_ElectricCharge = chg; }

		int Strangeness() const { return m_Strangeness; }
		void SetStrangenessCharge(int chg) { m_Strangeness = chg; }

		int Charm() const { return m_Charm; }
		void SetCharm(int chg) { m_Charm = chg; }

		double ArbitraryCharge() const { return m_ArbitraryCharge; }
		void SetArbitraryCharge(double arbchg) { m_ArbitraryCharge = arbchg; }

		double AbsoluteQuark() const { return m_AbsQuark; }
		void SetAbsoluteQuark(double abschg) { m_AbsQuark = abschg; }

		double AbsoluteStrangeness() const { return m_AbsS; }
		void SetAbsoluteStrangeness(double abschg) { m_AbsS = abschg; }

		double AbsoluteCharm() const { return m_AbsC; }
		void SetAbsoluteCharm(double abschg) { m_AbsC = abschg; }

		double ResonanceWidth() const { return m_Width; }
		void SetResonanceWidth(double width);// { m_Width = width; }

		double DecayThresholdMass() const { return m_Threshold; }
		void SetDecayThresholdMass(double threshold);// { m_Threshold = threshold; }

		// Threshold calculated from daughter masses
		double DecayThresholdMassDynamical() const { return m_ThresholdDynamical; }
		void CalculateAndSetDynamicalThreshold();

		ResonanceWidthShape GetResonanceWidthShape() const { return m_ResonanceWidthShape;  }
		void SetResonanceWidthShape(ResonanceWidthShape shape);// { m_ResonanceWidthShape = shape; }

		ResonanceWidthIntegration GetResonanceWidthIntegrationType() const { return m_ResonanceWidthIntegrationType; }
		void SetResonanceWidthIntegrationType(ResonanceWidthIntegration type);// { m_ResonanceWidthIntegrationType = type; }

		// Resonance Mass Distribution: Relativistic or Non-relativistic Breit-Wigner
		double MassDistribution(double m) const;

		// Resonance Mass Distribution with energy dependent width
		double MassDistribution(double m, double width) const;

		double Weight() const { return m_Weight; }
		void SetWeight(double weight) { m_Weight = weight; }

		int DecayType() const { return m_DecayType; }
		void SetDecayType(int type) { m_DecayType = type; }		

		const std::vector<ParticleDecay>& Decays() const { return m_Decays; }
				  std::vector<ParticleDecay>& Decays()       { return m_Decays; }
		void SetDecays(const std::vector<ParticleDecay> &Decays) { m_Decays = Decays; }
		void ClearDecays() { m_Decays.resize(0); }

		const std::vector<ParticleDecay>& DecaysOriginal() const		{ return m_DecaysOrig; }
				std::vector<ParticleDecay>& DecaysOriginal()					{ return m_DecaysOrig; }
		void SetDecaysOriginal(const std::vector<ParticleDecay> &DecaysOrig) { m_DecaysOrig = DecaysOrig; }

		void ReadDecays(std::string filename = "");

		const std::vector< std::pair<double, int> >&		DecayContributions()				const { return m_DecayContributions;					}
		std::vector< std::pair<double, int> >&					DecayContributions()							{ return m_DecayContributions; }

		const std::vector< std::pair<double, int> >& WeakDecayContributions()		const { return m_WeakDecayContributions;			}
		std::vector< std::pair<double, int> >& WeakDecayContributions()								{ return m_WeakDecayContributions; }

		const std::vector< std::pair<double, int> >& DecayContributionsSigmas() const { return m_DecayContributionsSigmas;		}
		std::vector< std::pair<double, int> >& DecayContributionsSigmas()							{ return m_DecayContributionsSigmas;		}

		const std::vector< std::pair< std::vector<double>, int> >& DecayCumulants()			const		{ return m_DecayCumulants; }
		std::vector< std::pair< std::vector<double>, int> >& DecayCumulants()										{ return m_DecayCumulants; }

		const std::vector< std::pair< std::vector<double>, int> >& DecayProbabilities() const		{ return m_DecayProbabilities; }
		std::vector< std::pair< std::vector<double>, int> >& DecayProbabilities()								{ return m_DecayProbabilities; }

		const std::vector< std::pair<double, std::vector<int> > >& DecayDistributions() const { return m_DecayDistributions; }
		std::vector< std::pair<double, std::vector<int> > >& DecayDistributions() { return m_DecayDistributions; }

		void CalculateThermalBranchingRatios(const ThermalModelParameters &params, bool useWidth = 0, double pMu = 0., double dMu = 0.);

		void SetCalculationType(IdealGasFunctions::QStatsCalculationType type)					{ m_QuantumStatisticsCalculationType = type; }
		IdealGasFunctions::QStatsCalculationType CalculationType()																const { return m_QuantumStatisticsCalculationType; }

		void SetClusterExpansionOrder(int order) { m_ClusterExpansionOrder = order; }
		int ClusterExpansionOrder() const { return m_ClusterExpansionOrder;  }

		std::vector<double> BranchingRatioWeights(const std::vector<double> & ms) const;


		const std::vector<double>& Nch() const { return m_Nch; }
		      std::vector<double>&  Nch()       { return m_Nch; }

		const std::vector<double>& DeltaNch() const { return m_DeltaNch; }
					std::vector<double>&  DeltaNch()       { return m_DeltaNch; }

	private:
		/**
		*	Auxiliary coefficients used for numerical integration using quadratures
		*/
		std::vector<double> m_xlag32, m_wlag32;
		std::vector<double> m_xleg,   m_wleg;
		std::vector<double> m_xleg32, m_wleg32;
		std::vector<double> m_brweight;


		/**
		*	For the eBW scheme
		*/
		std::vector<double> m_xlegdyn, m_wlegdyn, m_vallegdyn;
		std::vector<double> m_xlegpdyn, m_wlegpdyn, m_vallegpdyn;
		std::vector<double> m_xlagdyn, m_wlagdyn, m_vallagdyn;

		std::vector<double> m_xalldyn, m_walldyn, m_densalldyn;


		bool m_Stable;								/**< Flag whether particle is stable. */
		bool m_AntiParticle;					/**< Whether particle was created as an antiparticle to another one. */
		std::string m_Name;					/**< Particle name. */
		int m_PDGID;									/**< PDG (HEP) ID of a particle. */
		double m_Degeneracy;					/**< (Spin) Degeneracy factor. */
		int m_Statistics;						/**< Statistics used (Bolzmann or Quantum). */
		int m_StatisticsOrig;				/**< Particle's original Fermi/Bose statistics. */
		double m_Mass;								/**< Mass (GeV) */

		/** 
		*   0 - Cluster expansion (default), 1 - Numerical integration 
		*/
		IdealGasFunctions::QStatsCalculationType m_QuantumStatisticsCalculationType;	
		/** 
		*   Number of terms in cluster expansion.
		*   Default is 10 for m < 200 MeV (pions), 5 for m < 1000 MeV, and 3 otherwise
		*   
		*/
		int m_ClusterExpansionOrder;

		int m_Baryon;								/**< Baryon charge */
		int m_ElectricCharge;				/**< Electric charge */
		int m_Strangeness;						/**< Strangeness charge */
		int m_Charm;									/**< Charm charge */
		int m_Quark;									/**< Absolute quarks content */

		double m_ArbitraryCharge;		/**< Arbitrary charge associated with particle (external) */
		double m_AbsQuark;						/**< Absolute light quark content */
		double m_AbsS;								/**< Absolute strangeness content */
		double m_AbsC;								/**< Absolute charm content */

		double m_Width;							/**< Resonance width (GeV) */
		double m_Threshold;					/**< Lower decay threshold (GeV) */
		double m_ThresholdDynamical; /**< Lower decay threshold (GeV) calculated from daughter masses */
		ResonanceWidthShape m_ResonanceWidthShape;									/**< Either relativistic or non-relativitic Breit-Wigner */
		ResonanceWidthIntegration m_ResonanceWidthIntegrationType;	/**< Plus-minus TwoGamma or from m0 to infty */
		double m_Radius;							/**< Hard-core radius (fm) */
		double m_Vo;									/**< Eigenvolume parameter (fm^3). Obsolete. To be removed. */
		double m_Weight;							/**< Weight of a given particle. Default is 1 */

		int m_DecayType;													/**< 0 - stable, 1 - strong, 2 - weak	*/
		std::vector<ParticleDecay> m_Decays;			/**< All decay channels currently in use.	*/
																						
																						/**
																						*   All original decay channels. Needed to switch between renormalizing decay branching ratios and not.
																						*/
		std::vector<ParticleDecay> m_DecaysOrig;

		/**
		*   Contains information about decay chains of heavier particles resulting in production of a present particle.
		*   Contains indexes (0-based) and average yields resulting from corresponding decay chains.
		*/
		std::vector< std::pair<double, int> > m_DecayContributions;

		/**
		*   Contains information about decay chains of heavier particles resulting in production of a present particle.
		*   Contains indexes (0-based) and variance of yields resulting from corresponding decay chains.
		*/
		std::vector< std::pair<double, int> > m_DecayContributionsSigmas;

		/**
		*   Contains information about decay chains of heavier particles including weak decay feeddown resulting in production of a present particle.
		*   Contains indexes (0-based) and average yields resulting from corresponding decay chains.
		*/
		std::vector< std::pair<double, int> > m_WeakDecayContributions;

		/**
		*   Contains information about decay chains of heavier particles resulting in production of a present particle.
		*   Contains indexes (0-based) and first 4 moments (cumulants) of yields resulting from corresponding decay chains.
		*/
		std::vector< std::pair< std::vector<double>, int> > m_DecayCumulants;
		std::vector< std::pair< std::vector<double>, int> > m_DecayProbabilities;

		/**
		*   Contains all possible configurations which result from decays of a particle
		*/
		std::vector< std::pair<double, std::vector<int> > > m_DecayDistributions;

		/**
		*   For calculating final state charged particle multiplicities
		*   Arrays contain Nch, N+, N-
		*/
		std::vector<double> m_Nch;
		std::vector<double> m_DeltaNch;
};

#endif
