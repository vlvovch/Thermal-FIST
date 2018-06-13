#ifndef SIMPLEEVENT_H
#define SIMPLEEVENT_H

#include <cmath>
#include <vector>
#include <fstream>

#include "HRGEventGenerator/SimpleParticle.h"

struct SimpleEvent {
	double weight, logweight;
	std::vector<SimpleParticle> Particles;
	SimpleEvent() { Particles.resize(0); weight = 1.; }
	void writeToFile(std::ofstream & fout, int eventnumber = 1);
};

#endif
