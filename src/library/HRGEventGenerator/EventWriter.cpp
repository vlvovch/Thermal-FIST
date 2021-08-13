#include "HRGEventGenerator/EventWriter.h"

#include <string> 
#include <iostream>
#include <iomanip>

#include "ThermalFISTConfig.h"

namespace thermalfist {

  EventWriter::EventWriter(const std::string& filename)
  {
    OpenFile(filename);
  }

  EventWriter::~EventWriter()
  {
    CloseFile();
  }

  bool EventWriter::OpenFile(const std::string& filename)
  {
    if (m_fout.is_open())
      CloseFile();

    m_fout.open(filename);

    if (!m_fout.is_open())
      return false;

    m_EventNumber = 0;

    return true;
  }

  void EventWriter::CloseFile()
  {
    if (m_fout.is_open()) {
      m_fout.close();
    }
  }

  bool EventWriter::WriteEvent(const SimpleEvent& evt)
  {
    if (!m_fout.is_open())
      return false;

    if (evt.weight != 1.0) {
      std::cout << "**WARNING** Writing a weighted event to a file. The information about the weight will be lost!" << std::endl;
    }

    m_fout << "Event " << ++m_EventNumber << std::endl;

    const int tabsize = 25;

    m_fout.precision(16);
    m_fout << std::scientific;

    m_fout << std::setw(tabsize) << "pdgid";

    m_fout << std::setw(tabsize) << "p0[GeV/c2]";

    m_fout << std::setw(tabsize) << "px[GeV/c]"
      << std::setw(tabsize) << "py[GeV/c]"
      << std::setw(tabsize) << "pz[GeV/c]";

    m_fout << std::endl;

    for (size_t i = 0; i < evt.Particles.size(); ++i) {
      const SimpleParticle& part = evt.Particles[i];

      m_fout << std::setw(tabsize) << part.PDGID;

      m_fout << std::setw(tabsize) << part.p0;

      m_fout << std::setw(tabsize) << part.px
        << std::setw(tabsize) << part.py
        << std::setw(tabsize) << part.pz;

      m_fout << std::endl;
    }

    return true;
  }

  bool EventWriterAsciiExtended::WriteEvent(const SimpleEvent& evt)
  {
    if (!m_fout.is_open())
      return false;
    ++m_EventNumber; 
    evt.writeToFile(m_fout, m_config, m_EventNumber);
    return true;
  }

  bool EventWriterForUrqmd::WriteEvent(const SimpleEvent& evt)
  {
    if (!m_fout.is_open())
      return false;
    ++m_EventNumber;
    evt.writeToFileForUrqmd(m_fout);
    return true;
  }

}
