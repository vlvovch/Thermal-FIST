#include "HRGEventGenerator/HepMCEventWriter.h"

#include <string>
#include <iostream>

#include "ThermalFISTConfig.h"

namespace thermalfist {

  HepMCEventWriter::HepMCEventWriter(const std::string& filename)
  {
    OpenFile(filename);
  }

  HepMCEventWriter::~HepMCEventWriter()
  {
    CloseFile();
  }

  bool HepMCEventWriter::OpenFile(const std::string& filename)
  {
    if (!EventWriter::OpenFile(filename))
      return false;

    m_fout << std::scientific;
    m_fout.precision(16);

    if (!m_fout.is_open())
      return false;

    const std::string header = "HepMC::Version " + std::string("3.01.01") + "\nHepMC::Asciiv3-START_EVENT_LISTING\n";
    m_fout << header;

    m_fout << "W " << 1 << "\n";
    m_fout << "N " << "ImportanceSamplingWeight" << "\n";

    std::string version = std::to_string(ThermalFIST_VERSION_MAJOR) + "." + std::to_string(ThermalFIST_VERSION_MINOR);
    if (ThermalFIST_VERSION_DEVEL != 0) {
      version += "." + std::to_string(ThermalFIST_VERSION_DEVEL);
    }
    m_fout << "T " << "Thermal-FIST" << "\\|";
    m_fout << version << "\\|";
    m_fout << "ThermalFistEventGenerator" << "\n";

    return true;
  }

  void HepMCEventWriter::CloseFile()
  {
    if (m_fout.is_open()) {
      const std::string footer("HepMC::Asciiv3-END_EVENT_LISTING\n\n");
      m_fout << footer;
      m_fout.close();
    }
  }

  bool HepMCEventWriter::WriteEvent(const SimpleEvent& evt)
  {
    if (!m_fout.is_open())
      return false;

    m_fout << "E " << m_EventNumber++ << " " << 0 << " " << evt.Particles.size() << "\n";
    m_fout << "U " << "GEV" << " " << "MM" << "\n";
    m_fout << "W " << evt.weight << "\n";

    for (int ip = 0; ip < evt.Particles.size(); ++ip) {
      const SimpleParticle& part = evt.Particles[ip];
      m_fout << "P " << ip + 1 << " " << 0 << " ";

      m_fout << part.PDGID << " ";
      m_fout << part.px << " " << part.py << " " << part.pz << " " << part.p0 << " ";
      m_fout << part.m << " ";

      m_fout << 1 << "\n";
    }

    return true;
  }

}
