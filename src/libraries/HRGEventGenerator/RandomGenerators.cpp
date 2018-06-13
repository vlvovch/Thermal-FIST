#include "HRGEventGenerator/RandomGenerators.h"

#include "HRGBase/xMath.h"

namespace RandomGenerators {

	MTRand randgenMT;

	int RandomPoisson(double mean) {
	   int n;
	   if (mean <= 0) return 0;
	   if (mean < 25) {
		  double expmean = exp(-mean);
		  double pir = 1;
		  n = -1;
		  while(1) {
			 n++;
			 pir *= randgenMT.rand();
			 if (pir <= expmean) break;
		  }
		  return n;
	   }
	   // for large value we use inversion method
	   else {//if (mean < 1E9) {
		  double em, t, y;
		  double sq, alxm, g;
		  double pi = xMath::Pi();

		  sq = sqrt(2.0*mean);
		  alxm = log(mean);
		  g = mean*alxm - xMath::LogGamma(mean + 1.0);

		  do {
			 do {
				y = tan(pi*randgenMT.rand());
				em = sq*y + mean;
			 } while( em < 0.0 );

			 em = floor(em);
			 t = 0.9*(1.0 + y*y)* exp(em*alxm - xMath::LogGamma(em + 1.0) - g);
		  } while( randgenMT.rand() > t );

		  return static_cast<int> (em);

	   }
	   //else {
	   //   // use Gaussian approximation vor very large values
	   //   n = Int_t(Gaus(0,1)*TMath::Sqrt(mean) + mean + 0.5);
	   //   return n;
	   //}
	}


	int RandomPoisson(double mean,  MTRand &rangen) {
	   int n;
	   if (mean <= 0) return 0;
	   if (mean < 25) {
		  double expmean = exp(-mean);
		  double pir = 1;
		  n = -1;
		  while(1) {
			 n++;
			 pir *= rangen.rand();
			 if (pir <= expmean) break;
		  }
		  return n;
	   }
	   // for large value we use inversion method
	   else {//if (mean < 1E9) {
		  double em, t, y;
		  double sq, alxm, g;
		  double pi = xMath::Pi();

		  sq = sqrt(2.0*mean);
		  alxm = log(mean);
		  g = mean*alxm - xMath::LogGamma(mean + 1.0);

		  do {
			 do {
				y = tan(pi*rangen.rand());
				em = sq*y + mean;
			 } while( em < 0.0 );

			 em = floor(em);
			 t = 0.9*(1.0 + y*y)* exp(em*alxm - xMath::LogGamma(em + 1.0) - g);
		  } while( rangen.rand() > t );

		  return static_cast<int> (em);

	   }
	   //else {
	   //   // use Gaussian approximation vor very large values
	   //   n = Int_t(Gaus(0,1)*TMath::Sqrt(mean) + mean +0.5);
	   //   return n;
	   //}
	}


	double SiemensRasmussenGenerator::g(double x) const {
		double tp = -log(x);
		double talpha = alpha(tp);
		double en = w(tp);
		double sh = sinh(talpha);
		double shtalpha = 1.;
		if (talpha != 0.0)
			shtalpha = sh / talpha;
		double ch = sqrt(1. + sh*sh);
		return tp * tp * exp(-m_Gamma*en/m_T) * ( (1.+m_T/m_Gamma/en)*shtalpha - m_T/m_Gamma/en*ch ) / x;
	}

	double SiemensRasmussenGenerator::g2(double x, double tp) const {
		double talpha = alpha(tp);
		double en = w(tp);
		double sh = sinh(talpha);
		double shtalpha = 1.;
		if (talpha != 0.0)
			shtalpha = sh / talpha;
		double ch = sqrt(1. + sh*sh);
		return tp * tp * exp(-m_Gamma*en/m_T) * ( (1.+m_T/m_Gamma/en)*shtalpha - m_T/m_Gamma/en*ch ) / x;
	}

	void SiemensRasmussenGenerator::FixParameters() {
		double eps = 1e-8;
		double l = 0., r = 1.;
		double m1 = l + (r-l)/3.;
		double m2 = r - (r-l)/3.;
		int MAXITERS = 200;
		int iter = 0;
		while (fabs(m2-m1)>eps && iter<MAXITERS) {
			if (g(m1)<g(m2)) {
				l = m1;
			}
			else {
				r = m2;
			}
			m1 = l + (r-l)/3.;
			m2 = r - (r-l)/3.;
			iter++;
		}
		m_Max = g((m1+m2)/2.);
	}

	double SiemensRasmussenGenerator::GetRandom() const {
		while (1) {
			double x0 = randgenMT.randDblExc();
			double y0 = m_Max * randgenMT.randDblExc();
			if (y0<g(x0)) return -log(x0);
		}
		return 0.;
	}

	std::vector<double> SiemensRasmussenGenerator::GetMomentum() const {
		std::vector<double> ret(0);
		double tp   = GetRandom();
		double tphi = 2. * xMath::Pi() * randgenMT.rand();
		double cthe = 2. * randgenMT.rand() - 1.;
		double sthe = sqrt(1. - cthe*cthe);
		ret.push_back(tp*cos(tphi)*sthe); //px
		ret.push_back(tp*sin(tphi)*sthe); //py
		ret.push_back(tp*cthe);           //pz
		return ret;
	}



	void SSHGenerator::FixParameters() {
		{
			double eps = 1e-8;
			double l = 0., r = 1.;
			double m1 = l + (r-l)/3.;
			double m2 = r - (r-l)/3.;
			int MAXITERS = 200;
			int iter = 0;
			while (fabs(m2-m1)>eps && iter<MAXITERS) {
				if (g(m1)<g(m2)) {
					l = m1;
				}
				else {
					r = m2;
				}
				m1 = l + (r-l)/3.;
				m2 = r - (r-l)/3.;
				iter++;
			}
			m_MaxPt = g((m1+m2)/2.);

			m_dndpt.clear();
			double dx = m_dPt;
			for(double x = 0.5*dx; x<=1.; x += dx) {
				m_dndpt.add_val(x, g(x));
			}
		}

		{
			m_dndy.resize(0);
			m_MaxYs.resize(0);
			double dx = m_dPt;
			for(double x = 0.5*dx; x<=1.; x += dx) {
				
				double pt = -log(x);

				double eps = 1e-8;
				double l = -4. - m_EtaMax, r = 4. + m_EtaMax;
				double m1 = l + (r-l)/3.;
				double m2 = r - (r-l)/3.;
				int MAXITERS = 200;
				int iter = 0;
				while (fabs(m2-m1)>eps && iter<MAXITERS) {
					if (m_distr.dndy(m1,pt)<m_distr.dndy(m2,pt)) {
						l = m1;
					}
					else {
						r = m2;
					}
					m1 = l + (r-l)/3.;
					m2 = r - (r-l)/3.;
					iter++;
				}
				m_MaxYs.push_back(m_distr.dndy((m1+m2)/2.,pt));

				m_dndy.push_back(SplineFunction());
				double dy = m_dy;
				for(double ty = -4. - m_EtaMax + 0.5*dy; ty<=4. + m_EtaMax; ty += dy) {
					m_dndy[m_dndy.size()-1].add_val(ty, m_distr.dndy(ty,pt));
				}
			}
		}
	}

	void SSHGenerator::FixParameters2() {
		{
			double eps = 1e-8;
			double l = 0., r = 1.;
			double m1 = l + (r-l)/3.;
			double m2 = r - (r-l)/3.;
			int MAXITERS = 200;
			int iter = 0;
			while (fabs(m2-m1)>eps && iter<MAXITERS) {
				if (g(m1)<g(m2)) {
					l = m1;
				}
				else {
					r = m2;
				}
				m1 = l + (r-l)/3.;
				m2 = r - (r-l)/3.;
				iter++;
			}
			m_MaxPt = g((m1+m2)/2.);

			m_dndpt.clear();
			double dx = m_dPt;
			for(double x = 0.5*dx; x<=1.; x += dx) {
				m_dndpt.add_val(x, g(x));
			}
		}
		
		{
			m_dndy.resize(0);
			m_MaxYs.resize(0);
			double dx = m_dPt;
			for(double x = 0.5*dx; x <= 1.; x += dx) {
				
				double pt = -log(x);

				double eps = 1e-8;
				double l = -4. - m_EtaMax, r = 4. + m_EtaMax;
				double m1 = l + (r-l)/3.;
				double m2 = r - (r-l)/3.;
				int MAXITERS = 200;
				int iter = 0;
				while (fabs(m2-m1)>eps && iter<MAXITERS) {
					if (m_distr.dndysingle(m1,pt)<m_distr.dndysingle(m2,pt)) {
						l = m1;
					}
					else {
						r = m2;
					}
					m1 = l + (r-l)/3.;
					m2 = r - (r-l)/3.;
					iter++;
				}
				m_MaxYs.push_back(m_distr.dndysingle((m1+m2)/2.,pt));

				m_dndy.push_back(SplineFunction());
				double dy = m_dy;
				for(double ty = -4. + 0.5*dy; ty <= 4.; ty += dy) {
					m_dndy[m_dndy.size()-1].add_val(ty, m_distr.dndysingle(ty,pt));
				}
			}
		}

	}

	void SSHGenerator::FindMaximumPt() {
		
		double eps = 1e-8;
		double l = 0., r = 1.;
		double m1 = l + (r-l)/3.;
		double m2 = r - (r-l)/3.;
		int MAXITERS = 200;
		int iter = 0;
		while (fabs(m2-m1)>eps && iter<MAXITERS) {
			if (g(m1)<g(m2)) {
				l = m1;
			}
			else {
				r = m2;
			}
			m1 = l + (r-l)/3.;
			m2 = r - (r-l)/3.;
			iter++;
		}
		m_MaxPt = g((m1+m2)/2.);

		m_dndpt.clearall();
		double dx = 0.05;
		for(double x = 0.5*dx; x<=1.; x += dx) {
			m_dndpt.add_val(x, g(x));
		}
	}

	void SSHGenerator::FindMaximumY(double pt) {
		double eps = 1e-8;
		double l = -4. - m_EtaMax, r = 4. + m_EtaMax;
		double m1 = l + (r-l)/3.;
		double m2 = r - (r-l)/3.;
		int MAXITERS = 200;
		int iter = 0;
		while (fabs(m2-m1)>eps && iter<MAXITERS) {
			if (m_distr.dndy(m1,pt)<m_distr.dndy(m2,pt)) {
				l = m1;
			}
			else {
				r = m2;
			}
			m1 = l + (r-l)/3.;
			m2 = r - (r-l)/3.;
			iter++;
		}
		m_MaxY = m_distr.dndy((m1+m2)/2.,pt);
	}

	void SSHGenerator::FindMaximumY2(double pt) {
		double eps = 1e-8;
		double l = -4. - m_EtaMax, r = 4. + m_EtaMax;
		double m1 = l + (r-l)/3.;
		double m2 = r - (r-l)/3.;
		int MAXITERS = 200;
		int iter = 0;
		while (fabs(m2-m1)>eps && iter<MAXITERS) {
			if (m_distr.dndysingle(m1,pt)<m_distr.dndysingle(m2,pt)) {
				l = m1;
			}
			else {
				r = m2;
			}
			m1 = l + (r-l)/3.;
			m2 = r - (r-l)/3.;
			iter++;
		}
		m_MaxY = m_distr.dndysingle((m1+m2)/2.,pt);
	}

	std::pair<double, double> SSHGenerator::GetRandom() {
		double tpt = 0., ty = 0.;
		while (1) {
			double x0 = randgenMT.randDblExc();
			double y0 = m_MaxPt * randgenMT.randDblExc();
			if (y0<g2(x0)) {
				tpt = -log(x0);
				break;
			}
		}
		while (1) {
			int ind = (int)(exp(-tpt)/m_dPt);
			if (ind<0) ind = 0;
			if (ind>=m_dndy.size()) ind = m_dndy.size() - 1;
			double x0 = - 4. - m_EtaMax + (8. + 2. * m_EtaMax) * randgenMT.randDblExc();
			double y0 = m_MaxYs[ind] * randgenMT.randDblExc();
			if (y0<m_dndy[ind].f(x0)) {
				ty = x0;
				break;
			}
		}
		return std::make_pair(tpt, ty);
	}

	std::pair<double, double> SSHGenerator::GetRandom2() const {
		double tpt = 0., ty = 0., teta = 0.;
		while (1) {
			double x0 = randgenMT.randDblExc();
			double y0 = m_MaxPt * randgenMT.randDblExc();
			if (y0<g2(x0)) {
				tpt = -log(x0);
				break;
			}
		}
		while (1) {
			int ind = (int)(exp(-tpt)/m_dPt);
			if (ind<0) ind = 0;
			if (ind>=m_dndy.size()) ind = m_dndy.size() - 1;
			double x0 = - 4. + (8.) * randgenMT.randDblExc();
			double y0 = m_MaxYs[ind] * randgenMT.randDblExc();
			
			if (y0<m_dndy[ind].f(x0)) {
				ty = x0;
				teta = -m_EtaMax + 2. * m_EtaMax * randgenMT.randDblExc();
				break;
			}
		}
		return std::make_pair(tpt, ty-teta);
	}

	std::vector<double> SSHGenerator::GetMomentum() const {
		std::vector<double> ret(0);
		std::pair<double, double> pty = GetRandom2();
		double tpt  = pty.first;
		double ty   = pty.second;
		double tphi = 2. * xMath::Pi() * randgenMT.rand();
		ret.push_back(tpt*cos(tphi));                          //px
		ret.push_back(tpt*sin(tphi));                          //py
		ret.push_back(sqrt(tpt*tpt + m_Mass*m_Mass)*sinh(ty)); //pz
		return ret;
	}

	double BreitWignerGenerator::f(double x) const {
		return x / ((x*x-m_M*m_M)*(x*x-m_M*m_M) + m_M*m_M*m_Gamma*m_Gamma);
	}

	void BreitWignerGenerator::FixParameters() {
		double eps = 1e-8;
		double l = m_Mthr, r = m_M + 2.*m_Gamma;
		double m1 = l + (r-l)/3.;
		double m2 = r - (r-l)/3.;
		int MAXITERS = 200;
		int iter = 0;
		while (fabs(m2-m1)<eps && iter<MAXITERS) {
			if (f(m1)<f(m2)) {
				l = m1;
			}
			else {
				r = m2;
			}
			m1 = l + (r-l)/3.;
			m2 = r - (r-l)/3.;
			iter++;
		}
		m_Max = f((m1+m2)/2.);
	}

	double BreitWignerGenerator::GetRandom() const {
		//bool fl = true;
		if (m_Gamma<1e-7) return m_M;
		while (1) {
			double x0 = m_Mthr + (m_M + 2.*m_Gamma - m_Mthr) * randgenMT.rand();
			double y0 = m_Max * randgenMT.rand();
			if (y0<f(x0)) return x0;
		}
		return 0.;
	}

	void BreitWignerGenerator::SetParameters(double M, double gamma, double mthr) {
		m_M = M;
		m_Gamma = gamma;
		m_Mthr = mthr;
		FixParameters();
	}

}
