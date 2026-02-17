/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef CONFIGINTERACTIONS_H
#define CONFIGINTERACTIONS_H

#include <QWidget>
#include <QDialog>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QTableView>
#include <QSpinBox>
#include <QDoubleSpinBox>
#include <QComboBox>
#include <QCheckBox>
#include <QRadioButton>
#include <QTextEdit>
#include <QLayout>
#include <QGroupBox>

#include "HRGBase/ThermalModelBase.h"
#include "HRGFit/ThermalModelFit.h"
//#include "BaseStructures.h"

class TableModel;
class FitToExperimentTab;
class ModelConfigWidget;
class ThermalModelConfig;
class QTableWidget;

struct QvdWParameters {
  std::vector<std::vector<double>> m_bij;
  std::vector<std::vector<double>> m_aij;
  QvdWParameters(int N = 0) {
    m_aij = m_bij = std::vector<std::vector<double>>(N, std::vector<double>(N, 0.));
  }
  static QvdWParameters GetParameters(thermalfist::ThermalParticleSystem* TPS, const ThermalModelConfig* config);
};

class InteractionsDialog : public QDialog
{
  Q_OBJECT
public:
  ModelConfigWidget* m_parent;
  QRadioButton* radSet, *radMB, * radLoad;
  QDoubleSpinBox* spinB;
  QDoubleSpinBox* spinA;
  QLabel* labelRadiusValue;
  //QDoubleSpinBox *spinRadius;
  QComboBox* comboScaling;
  QComboBox* comboEVprescr;
  QCheckBox* CBMM, * CBMB, * CBBaB, * CBBB;

  QDoubleSpinBox* spinB_BB, *spinB_BaB, *spinB_MB, *spinB_MM;
  QDoubleSpinBox* spinA_BB, *spinA_BaB, *spinA_MB, *spinA_MM;

  QLineEdit* leFilePath;
  QPushButton* buttonChooseFile;

  QCheckBox* CBSearchMultipleSolutions;

  QGroupBox* groupMC;
  QCheckBox* CBEVMult, * CBEVCoord, * CBEVSPR;

public:
  explicit  InteractionsDialog(ModelConfigWidget* parent = 0);

private slots:
  void modeToggled();
  void chooseInputFile();
  void updateRadius();
  void updateSPR();
public slots:
  void OK();
  void Discard() { QDialog::reject(); };
};

class QvdWParametersTableDialog : public QDialog
{
  Q_OBJECT
  thermalfist::ThermalParticleSystem* TPS;
  ThermalModelConfig* config;
  QTableWidget* tablevdW;
  QComboBox* comboType;

public:
  explicit  QvdWParametersTableDialog(
    QWidget* parent = 0,
    ThermalModelConfig* iconfig = 0,
    thermalfist::ThermalParticleSystem* iTPS = NULL);

signals:

public slots:
  //void checkFixTableSize();
  void recalculate();

protected:
  bool eventFilter(QObject* obj, QEvent* ev);
};

#endif // CONFIGWIDGETS_H
