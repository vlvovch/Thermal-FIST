#ifndef SIMPLEPARTICLE_H
#define SIMPLEPARTICLE_H

#include <cmath>

struct SimpleParticle {
	double px,py,pz;
	double m;
	double p0;
	int PDGID;
	bool processed;
	SimpleParticle() { processed = false; }
	SimpleParticle(double px_, double py_, double pz_, double m_, int PDGID_):
		px(px_), 
		py(py_), 
		pz(pz_), 
		m(m_), 
		p0(sqrt(m*m+px*px+py*py+pz*pz)), 
		PDGID(PDGID_), processed(false) { }
	double GetP() const {
		return sqrt(p0*p0 - m*m);
	}
	double GetPt() const {
		return sqrt(px*px+py*py);
	}
	double GetMt() const {
		return sqrt(m*m+px*px+py*py);
	}
	double GetY() const {
		return 0.5*log((p0+pz)/(p0-pz));
	}
};

#endif
