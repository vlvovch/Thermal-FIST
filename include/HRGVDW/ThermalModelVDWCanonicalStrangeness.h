#ifndef THERMALMODELVDWCANONICALSTRANGENESS_H
#define THERMALMODELVDWCANONICALSTRANGENESS_H

#include <map>
#include <iostream>

#include "HRGBase/ThermalModelCanonicalStrangeness.h"
#include "HRGVDW/ThermalModelVDWFull.h"

class ThermalModelVDWCanonicalStrangeness : public ThermalModelCanonicalStrangeness
{
	public:
		ThermalModelVDWCanonicalStrangeness(ThermalParticleSystem *TPS_, const ThermalModelParameters& params = ThermalModelParameters());

		virtual ~ThermalModelVDWCanonicalStrangeness(void);

		void PrepareModelVDW();  /**< Creates a ThermalModelVDW copy */
		void CleanModelVDW();		/**< Cleares the ThermalModelVDW copy */

		void FillVirial(const std::vector<double> & ri = std::vector<double>(0));
		void FillVirialEV(const std::vector< std::vector<double> > & bij = std::vector< std::vector<double> >(0));
		void FillAttraction(const std::vector< std::vector<double> > & aij = std::vector< std::vector<double> >(0));
		void SetVirial(int i, int j, double b) { if (i >= 0 && i < m_Virial.size() && j >= 0 && j < m_Virial[i].size()) m_Virial[i][j] = b; }
		void SetAttraction(int i, int j, double a) { if (i >= 0 && i < m_Attr.size() && j >= 0 && j < m_Attr[i].size())     m_Attr[i][j] = a; }

		double VirialCoefficient(int i, int j) const;
		double AttractionCoefficient(int i, int j) const;

		virtual void ReadInteractionParameters(const std::string &filename);
		virtual void WriteInteractionParameters(const std::string &filename);

		virtual void CalculateDensitiesGCE();

		virtual void CalculateEnergyDensitiesGCE();

		virtual void CalculatePressuresGCE();

		virtual void CalculateSums(const std::vector<double> &  Vcs);

		virtual void CalculateDensities();

		virtual double CalculateEnergyDensity();

		virtual double CalculateEntropyDensity();

		virtual double CalculatePressure();

	protected:
		// TODO: test
		virtual double MuShift(int id);

		ThermalModelVDWFull *m_modelVDW;
		std::vector< std::vector<double> > m_Virial;
		std::vector< std::vector<double> > m_Attr;
		double m_PNS;		/**< Pressure of all non-strange hadrons */
		std::vector<double> m_MuStar;
		std::vector<double> m_Suppression;		/**< Common suppression factor, from non-strange hadrons */
		//double m_EVNS;		/**< Total eigenvolume of all non-strange hadrons */
		//double m_EVS;		/**< Total eigenvolume of all strange hadrons */
};

#endif // THERMALMODELEVCANONICALSTRANGENESS_H
