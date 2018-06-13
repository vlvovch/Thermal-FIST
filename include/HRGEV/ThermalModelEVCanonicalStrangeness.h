#ifndef THERMALMODELEVCANONICALSTRANGENESS_H
#define THERMALMODELEVCANONICALSTRANGENESS_H

#include <map>
#include <iostream>

#include "HRGBase/ThermalModelCanonicalStrangeness.h"
#include "HRGEV/ThermalModelEVDiagonal.h"


class ThermalModelEVCanonicalStrangeness : public ThermalModelCanonicalStrangeness
{
	public:
		ThermalModelEVCanonicalStrangeness(ThermalParticleSystem *TPS_, const ThermalModelParameters& params = ThermalModelParameters());

		virtual ~ThermalModelEVCanonicalStrangeness(void);

		void PrepareModelEV();  /**< Creates the ThermalModelEV copy */

		void CleanModelEV();		/**< Cleares the ThermalModelEV copy */

		void SetRadius(double rad);
		void FillVirial(const std::vector<double> & ri = std::vector<double>(0));
		void FillVirialEV(const std::vector<double> & vi = std::vector<double>(0));

		virtual void ReadInteractionParameters(const std::string &filename);
		virtual void WriteInteractionParameters(const std::string &filename);
		double ExcludedVolume(int i) const;// { return m_v[i]; }
		virtual double CalculateEigenvolumeFraction();
		void SetRadius(int i, double rad);

		virtual void CalculateDensitiesGCE();

		virtual void CalculateEnergyDensitiesGCE();

		virtual void CalculatePressuresGCE();

		virtual void CalculateDensities();

		virtual double CalculateEnergyDensity();

		virtual double CalculateEntropyDensity();

		virtual double CalculatePressure();

	protected:
		ThermalModelEVDiagonal *m_modelEV;
		std::vector<double> m_v;                       /**< Vector of eigenvolumes of all hadrons */
		double m_PNS;		/**< Pressure of all non-strange hadrons */
		double m_Suppression;		/**< Common suppression factor, from non-strange hadrons */
		double m_EVNS;		/**< Total eigenvolume of all non-strange hadrons */
		double m_EVS;		/**< Total eigenvolume of all strange hadrons */
};

#endif // THERMALMODELEVCANONICALSTRANGENESS_H
