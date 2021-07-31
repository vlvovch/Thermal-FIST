/*
 * Thermal-FIST package
 *
 * Copyright (c) 2019-2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef HEPMCEVENTWRITER_H
#define HEPMCEVENTWRITER_H

#include <fstream>

#include "SimpleEvent.h"


namespace thermalfist {

  /// \brief Writes the events in the HepMC::Asciiv3 format, see https://gitlab.cern.ch/hepmc/HepMC3
  /// 
  /// Assumes that all particles are final (status code = 1) and come from the root vertex.
  /// Has to be checked if this actually works
  class HepMCEventWriter
  {
  public:
    HepMCEventWriter(const std::string & filename = "");
    ~HepMCEventWriter();

    bool OpenFile(const std::string& filename);
    
    void CloseFile();

    bool WriteEvent(const SimpleEvent& evt);
  private:
    std::ofstream m_fout;
    int m_EventNumber;
  };

} // namespace thermalfist

#endif