/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef CONFIGWIDGETS_H
#define CONFIGWIDGETS_H

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
#include "BaseStructures.h"
//#include "fittoexperimenttab.h"

class TableModel;
class FitToExperimentTab;

class ModelConfigWidget : public QWidget
{
  Q_OBJECT
public:
  QComboBox *comboModel;

  QComboBox *comboEnsemble;

  //QPushButton *buttonStatistics;
  QRadioButton *radioBoltz, *radioQuant;
  QComboBox *comboQuant;
  QCheckBox *CBQuadratures;

  QPushButton* buttonQvdWparameters;

  QPushButton *buttonConservationLaws;

  QPushButton *buttonInteractions;

  QComboBox *comboWidth;

  QPushButton *buttonOther;

  ThermalModelConfig currentConfig;

  bool m_eventGeneratorMode;
  bool m_thermalFitMode;

  thermalfist::ThermalModelBase *model;


  ModelConfigWidget(QWidget *parent = 0, thermalfist::ThermalModelBase *model = NULL, bool eventGeneratorMode = false, bool thermalFitMode = false);
  ~ModelConfigWidget();
  ThermalModelConfig updatedConfig();
private slots:
  void conservationLawsDialog();
  void interactionsDialog();
  void QvdWparametersDialog();
  void otherOptionsDialog();
  void modelTypeChanged();
  void ensembleChanged();
  void statChanged();
public:
  void setModel(thermalfist::ThermalModelBase *model);
  void setNewConfig(const ThermalModelConfig &config);
signals:
  void changed();
};


class ConservationLawsDialog : public QDialog
{
  Q_OBJECT
public:
  ModelConfigWidget *m_parent;
  QCheckBox *CBmuB, *CBmuQ, *CBmuS, *CBmuC;
  QCheckBox *CBmuBdens;
  QDoubleSpinBox *spinSBRatio, *spinQBRatio;

  QDoubleSpinBox *spinRhoB;

  QWidget *CBmuBfull, *CBmuBdensfull, *CBmuQfull;
  //QHBoxLayout *laymuB, *laymuQ;
  QLabel *labelNothing;
  
  QCheckBox *checkBConserve;
  QCheckBox *checkQConserve;
  QCheckBox *checkSConserve;
  QCheckBox *checkCConserve;
public:
  explicit  ConservationLawsDialog(ModelConfigWidget *parent = 0);
private slots:
  void updateControls();
  void toggleMuB();

public slots :
  void OK();
  void Discard() { QDialog::reject(); };
};


class OtherOptionsDialog : public QDialog
{
  Q_OBJECT
public:
  ModelConfigWidget *m_parent;
  QCheckBox *CBNormBratio, *CBFluctuations;

  QComboBox *comboBWShape;

  QCheckBox *CBUsePCE;
  QDoubleSpinBox* spinTkin;
  QCheckBox *CBFreezeLongLived;
  QDoubleSpinBox *spinWidthCut;
  QCheckBox *CBSahaNuclei;

  QCheckBox *CBAnnihilation;
  QDoubleSpinBox *spinPionsAnnihilation;

  QDoubleSpinBox *spinB;
  QSpinBox *spinLandauLevels;

  QCheckBox *CBPionInteractions;
  QDoubleSpinBox *spinfPi;
  QCheckBox *CBKaonInteractions;
  QDoubleSpinBox *spinfKa;

public:
  explicit  OtherOptionsDialog(ModelConfigWidget *parent = 0);
signals:

public slots :
  void OK();
  void Discard() { QDialog::reject(); };
  void UpdateControls();
};

#endif // CONFIGWIDGETS_H
