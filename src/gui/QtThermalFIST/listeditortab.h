/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef LISTEDITORTAB_H
#define LISTEDITORTAB_H

#include <QMainWindow>
#include <QLineEdit>
#include <QPushButton>
#include <QTableView>
#include <QSpinBox>
#include <QDoubleSpinBox>
#include <QCheckBox>
#include <QLabel>

#include "HRGBase/ThermalModelBase.h"
#include "BaseStructures.h"

class ListTableModel;

class ListEditorTab : public QWidget
{
  Q_OBJECT

  QTableView  *tableParticles;
  QPushButton *buttonNewPart;
  QPushButton *buttonRemovePart;
  QPushButton *buttonResetAll;

  QLineEdit      *lePDGID;
  QLineEdit      *leName;
  QDoubleSpinBox *spinMass;
  QDoubleSpinBox *spinDegeneracy;

  QSpinBox *spinBaryon, *spinCharge, *spinStrangeness, *spinCharm;
  QDoubleSpinBox *spinAbsoluteStrangeness, *spinAbsoluteCharm;
  QLabel *labelAbsoluteLightQuarkValue;
  QCheckBox *checkBoxStable;
  QPushButton *buttonDecayModes;
  QDoubleSpinBox *spinWidth, *spinThreshold;
  QPushButton *buttonSetThreshold;
  QPushButton *buttonReset;

  QPushButton *buttonMassesWidthsFromPDG;
  QPushButton *buttonSetAllThresholds;

  QPushButton *buttonApply;
  QPushButton *buttonSaveToFile;

  ListTableModel *myModel;

  thermalfist::ThermalModelBase *model;

  QString currentPath;

  int getCurrentRow();

public:
  ListEditorTab(QWidget *parent = 0, thermalfist::ThermalModelBase *model_in = NULL);
  ~ListEditorTab();
private slots:
  void changedRow();
  void editDecays();
  void PDGEdited();
  void SpinEdited();
  void NameEdited();
  void MassEdited();
  void ChargeEdited();
  void DecaysEdited();
  void StableEdited();
  void editDecaysDoubleClick(const QModelIndex &);
  void resetParticle();
  void resetAll();
  void addParticle();
  void removeParticle();
  void setParticleList(thermalfist::ThermalParticleSystem *TPS);
  void setParticleMassThreshold(int id, int type);
  void setThreshold();
  void setAllThresholds();
  void saveToFile();
  void loadMassesWidthsFromPdg();
public slots:
  void applyChanges();
  bool haveChangesToList();
public:
  void resetTPS();
  void setListPath(const QString& path) { currentPath = path; }
};

#endif