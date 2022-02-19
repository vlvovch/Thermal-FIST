/*
 * Thermal-FIST package
 *
 * Copyright (c) 2019-2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef EVENTWRITER_H
#define EVENTWRITER_H

#include <fstream>

#include "SimpleEvent.h"


namespace thermalfist {

  /// \brief Base event writer class outputting the list of particle's four-momenta in a simple ascii format
  /// 
  class EventWriter
  {
  public:
    EventWriter(const std::string & filename = "");
    virtual ~EventWriter();

    virtual bool OpenFile(const std::string& filename);
    
    virtual void CloseFile();

    virtual bool WriteEvent(const SimpleEvent& evt);
  protected:
    std::ofstream m_fout;
    int m_EventNumber;
  };

  /// \brief Event writer class outputting the list of particles with options to include more information.
  ///
  ///  The possible extra information includes the event weight, particle coordinates, mother particle and decay epoch
  /// 
  class EventWriterAsciiExtended
    : public EventWriter
  {
  public:
    EventWriterAsciiExtended(const std::string& filename = "", const SimpleEvent::EventOutputConfig& config = SimpleEvent::EventOutputConfig()) : EventWriter(filename), m_config(config) { }
    void SetOutputConfig(const SimpleEvent::EventOutputConfig& config) { m_config = config; }
    bool WriteEvent(const SimpleEvent& evt);
  private:
    SimpleEvent::EventOutputConfig m_config;
  };

  /// \brief Event writer class outputting the list of particles suitable for UrQMD afterburner from https://github.com/jbernhard/urqmd-afterburner
  class EventWriterForUrqmd
    : public EventWriter
  {
  public:
    EventWriterForUrqmd(const std::string& filename = "") : EventWriter(filename) { }
    bool WriteEvent(const SimpleEvent& evt);
  };

} // namespace thermalfist



#endif