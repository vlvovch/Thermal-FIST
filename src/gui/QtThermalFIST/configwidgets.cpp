/*
 * Thermal-FIST package
 *
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "configwidgets.h"

#include <QLayout>
#include <QLabel>
#include <QHeaderView>
#include <QGroupBox>
#include <QApplication>
#include <QMessageBox>
#include <QElapsedTimer>
#include <QFileDialog>
#include <QDialogButtonBox>
#include <QDebug>
#include <cstdio>
#include <vector>

#include "ThermalFISTConfig.h"
#include "HRGBase/ThermalModelIdeal.h"
#include "HRGEV/ThermalModelEVDiagonal.h"
#include "HRGEV/ThermalModelEVCrossterms.h"
#include "HRGVDW/ThermalModelVDW.h"
#include "HRGBase/ThermalModelCanonical.h"
#include "HRGBase/ThermalModelCanonicalStrangeness.h"
#include "HRGEV/ThermalModelEVCanonicalStrangeness.h"
#include "HRGVDW/ThermalModelVDWCanonicalStrangeness.h"
#include "HRGBase/ThermalModelCanonicalCharm.h"
#include "HRGEV/ExcludedVolumeHelper.h"

#include "DebugText.h"
#include "tablemodel.h"
#include "particledialog.h"
#include "resultdialog.h"
#include "fittoexperimenttab.h"
#include "configinteractions.h"

using namespace thermalfist;

namespace {
  const QString ModelIdeal = "Ideal";
  const QString ModelDEV = "Excluded volume (Diagonal)";
  const QString ModelCRSEV = "Excluded volume (X-terms)";
  const QString ModelQvdW = "Quantum van der Waals";
  const QString ModelRealGas = "Quantum real gas";
}

ModelConfigWidget::ModelConfigWidget(QWidget* parent, ThermalModelBase* modelop, bool eventGeneratorMode, bool thermalFitMode)
  : QWidget(parent), m_eventGeneratorMode(eventGeneratorMode), m_thermalFitMode(thermalFitMode)
{
  model = modelop;

  currentConfig = ThermalModelConfig::fromThermalModel(model);

  QVBoxLayout* layout = new QVBoxLayout();
  layout->setContentsMargins(0, 0, 0, 0);
  layout->setAlignment(Qt::AlignTop);

  QHBoxLayout* layModelEnsemble = new QHBoxLayout();
  layModelEnsemble->setAlignment(Qt::AlignLeft);

  QLabel* labelModel = new QLabel(tr("Model:"));
  comboModel = new QComboBox();
  comboModel->addItem(ModelIdeal);
  comboModel->addItem(ModelDEV);
  comboModel->addItem(ModelCRSEV);
  comboModel->addItem(ModelQvdW);
  comboModel->addItem(ModelRealGas);
  comboModel->setCurrentIndex(0);
  comboModel->setToolTip(tr("Type of HRG model"));

  QLabel* labelEnsemble = new QLabel(tr("Ensemble:"));
  comboEnsemble = new QComboBox();
  comboEnsemble->addItem(tr("Grand-canonical"));
  comboEnsemble->addItem(tr("Canonical"));
  if (model->TPS()->hasStrange())
    comboEnsemble->addItem(tr("Strangeness-canonical"));
  if (model->TPS()->hasCharmed())
    comboEnsemble->addItem(tr("Charm-canonical"));
  comboEnsemble->setCurrentIndex(0);
  comboEnsemble->setToolTip(tr("Choice of ensemble"));

  layModelEnsemble->addWidget(labelModel);
  layModelEnsemble->addWidget(comboModel);
  layModelEnsemble->addSpacing(20);
  layModelEnsemble->addWidget(labelEnsemble);
  layModelEnsemble->addWidget(comboEnsemble);



  /// Quantum statistics
  QHBoxLayout* layStats = new QHBoxLayout();
  layStats->setAlignment(Qt::AlignLeft);
  QLabel* labelStats = new QLabel(tr("Statistics:"));
  radioBoltz = new QRadioButton(tr("Boltzmann"));
  radioQuant = new QRadioButton(tr("Quantum"));

  radioBoltz->setToolTip(tr("Maxwell-Boltzmann"));
  radioQuant->setToolTip(tr("Fermi-Dirac/Bose-Einstein"));

  QLabel* labelQuant = new QLabel(tr("for"));
  comboQuant = new QComboBox();
  comboQuant->addItem(tr("All particles"));
  comboQuant->addItem(tr("Mesons only"));
  comboQuant->addItem(tr("Pions only"));

  CBQuadratures = new QCheckBox(tr("Use quadratures"));

  CBQuadratures->setToolTip(tr("Use Gauss-Laguerre quadratures (slower but more reliable) or cluster expansion (faster but unreliable for large mu)"));

  if (!m_eventGeneratorMode)
    radioQuant->setChecked(true);
  else
    radioBoltz->setChecked(true);
  CBQuadratures->setChecked(true);



  if (1 || !m_eventGeneratorMode) {
    layStats->addWidget(labelStats);
    layStats->addWidget(radioBoltz);
    layStats->addWidget(radioQuant);
    layStats->addWidget(labelQuant);
    layStats->addWidget(comboQuant);
    layStats->addSpacing(10);
    layStats->addWidget(CBQuadratures);
  }


  QHBoxLayout* layOptions = new QHBoxLayout();
  layOptions->setAlignment(Qt::AlignLeft);

  //buttonStatistics = new QPushButton(tr("Quantum statistics..."));

  buttonConservationLaws = new QPushButton(tr("Conservation laws..."));
  connect(buttonConservationLaws, &QPushButton::clicked, this, &ModelConfigWidget::conservationLawsDialog);

  buttonInteractions = new QPushButton(tr("EV/vdW interactions..."));
  connect(buttonInteractions, &QPushButton::clicked, this, &ModelConfigWidget::interactionsDialog);

  buttonOther = new QPushButton(tr("PCE/Saha/Other..."));
  connect(buttonOther, &QPushButton::clicked, this, &ModelConfigWidget::otherOptionsDialog);


  layOptions->addWidget(buttonConservationLaws);
  if (1 || !m_eventGeneratorMode) {
    layOptions->addSpacing(20);
  }
  else {
    buttonConservationLaws->setVisible(false);
  }
  layOptions->addWidget(buttonInteractions);
  layOptions->addSpacing(20);
  layOptions->addWidget(buttonOther);




  QHBoxLayout* layOptions2 = new QHBoxLayout();
  layOptions2->setAlignment(Qt::AlignLeft);

  QLabel* labelWidth = new QLabel(tr("Resonance widths:"));
  comboWidth = new QComboBox();
  comboWidth->addItem(tr("Zero-width"));
  comboWidth->addItem(tr("Const Breit-Wigner"));
  comboWidth->addItem(tr("eBW"));
  comboWidth->addItem(tr("eBW (const BRs)"));
  comboWidth->setCurrentIndex(static_cast<int>(model->TPS()->ResonanceWidthIntegrationType()));
  comboWidth->setToolTip(tr("Prescription for the treatment of resonance widths"));

  buttonQvdWparameters = new QPushButton(tr("EV/vdW parameter list..."));
  connect(buttonQvdWparameters, &QPushButton::clicked, this, &ModelConfigWidget::QvdWparametersDialog);

  layOptions2->addWidget(labelWidth);
  layOptions2->addWidget(comboWidth);
  layOptions2->addWidget(buttonQvdWparameters);
  if (0 && m_eventGeneratorMode) {
    layOptions2->addSpacing(10);
    layOptions2->addWidget(labelStats);
    layOptions2->addWidget(radioBoltz);
    layOptions2->addWidget(radioQuant);
    layOptions2->addWidget(labelQuant);
    layOptions2->addWidget(comboQuant);
    layOptions2->addWidget(CBQuadratures);

    radioQuant->setChecked(false);
    radioBoltz->setChecked(true);
    radioQuant->setEnabled(false);
    labelQuant->setVisible(false);
    comboQuant->setVisible(false);
    CBQuadratures->setVisible(false);
  }

  layout->addLayout(layModelEnsemble);
  layout->addSpacing(10);
  layout->addLayout(layStats);
  layout->addSpacing(10);
  layout->addLayout(layOptions2);
  layout->addSpacing(10);
  layout->addLayout(layOptions);


  setLayout(layout);

  connect(comboModel, SIGNAL(currentIndexChanged(int)), this, SLOT(modelTypeChanged()));
  connect(comboEnsemble, SIGNAL(currentIndexChanged(int)), this, SLOT(ensembleChanged()));

  connect(radioBoltz, SIGNAL(toggled(bool)), this, SLOT(statChanged()));

  modelTypeChanged();
  ensembleChanged();
}

ModelConfigWidget::~ModelConfigWidget()
{
}



ThermalModelConfig ModelConfigWidget::updatedConfig()
{
  ThermalModelConfig ret = currentConfig;

  ret.ModelType = ThermalModelConfig::Ideal;
  if (comboEnsemble->currentText() == tr("Grand-canonical")) {
    if (comboModel->currentText() == ModelDEV)
      ret.ModelType = ThermalModelConfig::DiagonalEV;
    if (comboModel->currentText() == ModelCRSEV)
      ret.ModelType = ThermalModelConfig::CrosstermsEV;
    if (comboModel->currentText() == ModelQvdW)
      ret.ModelType = ThermalModelConfig::QvdW;
    if (comboModel->currentText() == ModelRealGas)
      ret.ModelType = ThermalModelConfig::RealGas;
  }

  if (comboEnsemble->currentText() == tr("Canonical")) {
    ret.ModelType = ThermalModelConfig::CE;
  }

  if (comboEnsemble->currentText() == tr("Strangeness-canonical")) {
    ret.ModelType = ThermalModelConfig::SCE;
    if (comboModel->currentText() == ModelDEV)
      ret.ModelType = ThermalModelConfig::EVSCE;
    if (comboModel->currentText() == ModelQvdW)
      ret.ModelType = ThermalModelConfig::VDWSCE;
  }

  if (comboEnsemble->currentText() == tr("Charm-canonical")) {
    ret.ModelType = ThermalModelConfig::CCE;
  }

  ret.InteractionModel = ThermalModelConfig::InteractionIdeal;
  if (ret.ModelType == ThermalModelConfig::DiagonalEV
    || ret.ModelType == ThermalModelConfig::EVSCE) {
    ret.InteractionModel = ThermalModelConfig::InteractionEVDiagonal;
  }
  if (ret.ModelType == ThermalModelConfig::CrosstermsEV) {
    ret.InteractionModel = ThermalModelConfig::InteractionEVCrossterms;
  }
  if (ret.ModelType == ThermalModelConfig::QvdW
    || ret.ModelType == ThermalModelConfig::VDWSCE) {
    ret.InteractionModel = ThermalModelConfig::InteractionQVDW;
  }
  if (ret.ModelType == ThermalModelConfig::RealGas) {
    ret.InteractionModel = ThermalModelConfig::InteractionRealGas;
  }

  // Event generator mode
  if (m_eventGeneratorMode) {
    ret.InteractionModel = ThermalModelConfig::InteractionIdeal;

    if (comboModel->currentText() == ModelDEV)
      ret.InteractionModel = ThermalModelConfig::InteractionEVDiagonal;
    if (comboModel->currentText() == ModelCRSEV)
      ret.InteractionModel = ThermalModelConfig::InteractionEVCrossterms;
    if (comboModel->currentText() == ModelQvdW)
      ret.InteractionModel = ThermalModelConfig::InteractionQVDW;
    if (comboModel->currentText() == ModelRealGas)
      ret.InteractionModel = ThermalModelConfig::InteractionRealGas;
  }

  ret.Ensemble = ThermalModelConfig::EnsembleGCE;
  if (ret.ModelType == ThermalModelConfig::CE)
    ret.Ensemble = ThermalModelConfig::EnsembleCE;
  if (ret.ModelType == ThermalModelConfig::CCE)
    ret.Ensemble = ThermalModelConfig::EnsembleCCE;
  if (ret.ModelType == ThermalModelConfig::SCE
    || ret.ModelType == ThermalModelConfig::EVSCE
    || ret.ModelType == ThermalModelConfig::VDWSCE)
    ret.Ensemble = ThermalModelConfig::EnsembleSCE;




  ret.QuantumStatistics = static_cast<int>(!radioBoltz->isChecked());
  ret.QuantumStatisticsType = static_cast<int>(CBQuadratures->isChecked());
  ret.QuantumStatisticsInclude = comboQuant->currentIndex();

  ret.FiniteWidth = comboWidth->currentIndex();

  if (!m_eventGeneratorMode) {
    ret.UsePCE &= (ret.Ensemble == ret.EnsembleGCE);
  }

  return (currentConfig = ret);
}

void ModelConfigWidget::conservationLawsDialog()
{
  currentConfig = updatedConfig();
  ConservationLawsDialog dialog(this);
  dialog.setWindowFlags(Qt::Window);
  dialog.exec();
  emit changed();
}

void ModelConfigWidget::interactionsDialog()
{
  currentConfig = updatedConfig();
  InteractionsDialog dialog(this);
  dialog.setWindowFlags(Qt::Window);
  dialog.exec();
  emit changed();
}

void ModelConfigWidget::QvdWparametersDialog()
{
  currentConfig = updatedConfig();
  QvdWParametersTableDialog dialog(this, &currentConfig, model->TPS());
  dialog.setWindowFlags(Qt::Window);
  dialog.showMaximized();
  dialog.exec();
  emit changed();
}

void ModelConfigWidget::otherOptionsDialog()
{
  currentConfig = updatedConfig();
  OtherOptionsDialog dialog(this);
  dialog.setWindowFlags(Qt::Window);
  dialog.exec();
  emit changed();
}

void ModelConfigWidget::modelTypeChanged()
{
  currentConfig = updatedConfig();

  QVector<QString> ensitems;
  for (int i = 0; i < comboEnsemble->count(); ++i)
    ensitems.push_back(comboEnsemble->itemText(i));
  QString currentText = comboEnsemble->currentText();

  QVector<QString> newitems;
  newitems.push_back(tr("Grand-canonical"));

  if (comboModel->currentText() == ModelIdeal || m_eventGeneratorMode)
    newitems.push_back(tr("Canonical"));

  if (comboModel->currentText() == ModelIdeal
    || comboModel->currentText() == ModelDEV
    || comboModel->currentText() == ModelQvdW || m_eventGeneratorMode) {
    if (model->TPS()->hasStrange())
      newitems.push_back(tr("Strangeness-canonical"));
  }

  if (comboModel->currentText() == ModelIdeal || m_eventGeneratorMode) {
    if (model->TPS()->hasCharmed())
      newitems.push_back(tr("Charm-canonical"));
  }

  if (newitems != ensitems) {
    comboEnsemble->blockSignals(true);
    comboEnsemble->clear();
    comboEnsemble->addItems(QStringList::fromVector(newitems));
    int index = comboEnsemble->findData(currentText, Qt::DisplayRole);
    if (index != -1)
      comboEnsemble->setCurrentIndex(index);
    comboEnsemble->blockSignals(false);
  }

  if (comboModel->currentText() == ModelIdeal) {
    buttonInteractions->setEnabled(false);
    buttonQvdWparameters->setEnabled(false);
  }
  else {
    buttonInteractions->setEnabled(true);
    buttonQvdWparameters->setEnabled(true);
  }

  emit changed();
}

void ModelConfigWidget::ensembleChanged()
{
  currentConfig = updatedConfig();

  QVector<QString> olditems;
  for (int i = 0; i < comboModel->count(); ++i)
    olditems.push_back(comboModel->itemText(i));
  QString currentText = comboModel->currentText();

  QVector<QString> newitems;
  newitems.push_back(ModelIdeal);

  if (comboEnsemble->currentText() == tr("Grand-canonical")
    || comboEnsemble->currentText() == tr("Strangeness-canonical") || m_eventGeneratorMode)
    newitems.push_back(ModelDEV);

  if (comboEnsemble->currentText() == tr("Grand-canonical") || m_eventGeneratorMode)
    newitems.push_back(ModelCRSEV);

  if (comboEnsemble->currentText() == tr("Grand-canonical")
    || comboEnsemble->currentText() == tr("Strangeness-canonical") || m_eventGeneratorMode)
    newitems.push_back(ModelQvdW);

  if (comboEnsemble->currentText() == tr("Grand-canonical") || m_eventGeneratorMode)
    newitems.push_back(ModelRealGas);

  if (newitems != olditems) {
    comboModel->blockSignals(true);
    comboModel->clear();
    comboModel->addItems(QStringList::fromVector(newitems));
    int index = comboModel->findData(currentText, Qt::DisplayRole);
    if (index != -1)
      comboModel->setCurrentIndex(index);
    comboModel->blockSignals(false);
  }

  //if (comboEnsemble->currentText() == tr("Canonical"))
  //  buttonConservationLaws->setEnabled(false);
  //else
  //  buttonConservationLaws->setEnabled(true);

  if (comboEnsemble->currentText() == tr("Canonical")) {
    CBQuadratures->setChecked(false);
    CBQuadratures->setEnabled(false);
  }
  else if (!radioBoltz->isChecked()) {
    CBQuadratures->setEnabled(true);
  }

  if (!(comboEnsemble->currentText() == tr("Canonical")) && m_eventGeneratorMode)
    buttonConservationLaws->setEnabled(false);
  else
    buttonConservationLaws->setEnabled(true);

  emit changed();
}

void ModelConfigWidget::statChanged()
{
  if (radioBoltz->isChecked()) {
    comboQuant->setEnabled(false);
    CBQuadratures->setEnabled(false);
  }
  else {
    comboQuant->setEnabled(true);
    if (comboEnsemble->currentText() == tr("Canonical")) {
      CBQuadratures->setEnabled(false);
    }
    else {
      CBQuadratures->setEnabled(true);
    }
  }
}

void ModelConfigWidget::setModel(ThermalModelBase* modelop) {
  model = modelop;
  modelTypeChanged();
  ensembleChanged();
}

void ModelConfigWidget::setNewConfig(const ThermalModelConfig& config)
{
  currentConfig = config;

  comboModel->blockSignals(true);
  comboModel->clear();
  comboModel->addItem(ModelIdeal);
  comboModel->addItem(ModelDEV);
  comboModel->addItem(ModelCRSEV);
  comboModel->addItem(ModelQvdW);
  comboModel->addItem(ModelRealGas);
  comboModel->setCurrentIndex(0);
  if (currentConfig.InteractionModel == ThermalModelConfig::DiagonalEV)
    comboModel->setCurrentIndex(1);
  if (currentConfig.InteractionModel == ThermalModelConfig::CrosstermsEV)
    comboModel->setCurrentIndex(2);
  if (currentConfig.InteractionModel == ThermalModelConfig::QvdW)
    comboModel->setCurrentIndex(3);
  if (currentConfig.InteractionModel == ThermalModelConfig::RealGas)
    comboModel->setCurrentIndex(4); 
  comboModel->blockSignals(false);


  comboEnsemble->blockSignals(true);
  comboEnsemble->clear();
  comboEnsemble->addItem(tr("Grand-canonical"));
  comboEnsemble->addItem(tr("Canonical"));
  if (model->TPS()->hasStrange())
    comboEnsemble->addItem(tr("Strangeness-canonical"));
  if (model->TPS()->hasCharmed())
    comboEnsemble->addItem(tr("Charm-canonical"));
  comboEnsemble->setCurrentIndex(0);
  if (currentConfig.Ensemble == ThermalModelConfig::EnsembleCE)
    comboEnsemble->setCurrentIndex(1);
  if (currentConfig.Ensemble == ThermalModelConfig::EnsembleSCE)
    comboEnsemble->setCurrentIndex(2);
  if (currentConfig.Ensemble == ThermalModelConfig::EnsembleCCE)
    comboEnsemble->setCurrentIndex(3);
  comboEnsemble->blockSignals(false);

  if (config.QuantumStatistics && !m_eventGeneratorMode)
    radioQuant->setChecked(true);
  else
    radioBoltz->setChecked(true);

  comboQuant->setCurrentIndex(config.QuantumStatisticsInclude);

  CBQuadratures->setChecked(config.QuantumStatisticsType);

  comboWidth->setCurrentIndex(config.FiniteWidth);

  modelTypeChanged();
  ensembleChanged();
  statChanged();
}

ConservationLawsDialog::ConservationLawsDialog(ModelConfigWidget* parent) : QDialog(parent), m_parent(parent)
{
  QVBoxLayout* layout = new QVBoxLayout();

  QGroupBox* grLaws = new QGroupBox(tr("Conservation laws"));

  QVBoxLayout* grLayout = new QVBoxLayout();

  QHBoxLayout* laymuB = new QHBoxLayout();
  laymuB->setAlignment(Qt::AlignLeft);
  CBmuB = new QCheckBox(tr("Constrain μB from entropy per baryon ratio, S/B:"));
  CBmuB->setChecked(m_parent->currentConfig.ConstrainMuB && m_parent->currentConfig.ConstrainMuBType == 0);
  connect(CBmuB, &QCheckBox::toggled, this, &ConservationLawsDialog::updateControls);
  spinSBRatio = new QDoubleSpinBox();
  spinSBRatio->setDecimals(3);
  spinSBRatio->setMinimum(0.001);
  spinSBRatio->setMaximum(100000.);
  spinSBRatio->setSingleStep(10.);
  spinSBRatio->setValue(m_parent->currentConfig.SoverB);

  laymuB->addWidget(CBmuB);
  laymuB->addWidget(spinSBRatio);
  laymuB->setContentsMargins(0, 0, 0, 0);


  CBmuBfull = new QWidget();
  CBmuBfull->setLayout(laymuB);


  QHBoxLayout* laymuBdens = new QHBoxLayout();
  laymuBdens->setAlignment(Qt::AlignLeft);
  CBmuBdens = new QCheckBox(tr("Constrain μB from baryon density, nB (fm^-3):"));
  CBmuBdens->setChecked(m_parent->currentConfig.ConstrainMuB && m_parent->currentConfig.ConstrainMuBType == 1);
  connect(CBmuBdens, &QCheckBox::toggled, this, &ConservationLawsDialog::toggleMuB);
  spinRhoB = new QDoubleSpinBox();
  spinRhoB->setDecimals(3);
  spinRhoB->setMinimum(-100.);
  spinRhoB->setMaximum(100.);
  spinRhoB->setSingleStep(0.01);
  spinRhoB->setValue(m_parent->currentConfig.RhoB);

  laymuBdens->addWidget(CBmuBdens);
  laymuBdens->addWidget(spinRhoB);
  laymuBdens->setContentsMargins(0, 0, 0, 0);


  CBmuBdensfull = new QWidget();
  CBmuBdensfull->setLayout(laymuBdens);

  QHBoxLayout* laymuQ = new QHBoxLayout();
  laymuQ->setAlignment(Qt::AlignLeft);
  CBmuQ = new QCheckBox(tr("Constrain μQ from electric-to-baryon charge ratio, Q/B:"));
  CBmuQ->setChecked(m_parent->currentConfig.ConstrainMuQ);
  connect(CBmuQ, &QCheckBox::toggled, this, &ConservationLawsDialog::updateControls);
  spinQBRatio = new QDoubleSpinBox();
  spinQBRatio->setDecimals(3);
  spinQBRatio->setMinimum(0.);
  spinQBRatio->setSingleStep(0.1);
  spinQBRatio->setValue(m_parent->currentConfig.QoverB);

  laymuQ->addWidget(CBmuQ);
  laymuQ->addWidget(spinQBRatio);
  laymuQ->setContentsMargins(0, 0, 0, 0);

  CBmuQfull = new QWidget();
  CBmuQfull->setLayout(laymuQ);

  CBmuS = new QCheckBox(tr("Constrain μS from the condition of zero net strangeness"));
  CBmuS->setChecked(m_parent->currentConfig.ConstrainMuS);

  CBmuC = new QCheckBox(tr("Constrain μC from the condition of zero net charm"));
  CBmuC->setChecked(m_parent->currentConfig.ConstrainMuC);

  //if (m_parent->model->TPS()->hasBaryons()
  //  && (m_parent->currentConfig.Ensemble != ThermalModelConfig::EnsembleCE ||
  //    !m_parent->currentConfig.CanonicalB))
  grLayout->addWidget(CBmuBfull);
  grLayout->addWidget(CBmuBdensfull);

  //if (m_parent->model->TPS()->hasCharged()
  //  && (m_parent->currentConfig.Ensemble != ThermalModelConfig::EnsembleCE ||
  //    !m_parent->currentConfig.CanonicalQ))
  grLayout->addWidget(CBmuQfull);

  //if (m_parent->model->TPS()->hasStrange() 
  //  && (m_parent->currentConfig.Ensemble != ThermalModelConfig::EnsembleCE ||
  //    !m_parent->currentConfig.CanonicalS)
  //  && m_parent->currentConfig.Ensemble != ThermalModelConfig::EnsembleSCE)
  grLayout->addWidget(CBmuS);

  //if (m_parent->model->TPS()->hasCharmed() 
  //  && (m_parent->currentConfig.Ensemble != ThermalModelConfig::EnsembleCE ||
  //    !m_parent->currentConfig.CanonicalC)
  //  && m_parent->currentConfig.Ensemble != ThermalModelConfig::EnsembleSCE
  //  && m_parent->currentConfig.Ensemble != ThermalModelConfig::EnsembleCCE)
  grLayout->addWidget(CBmuC);

  labelNothing = new QLabel(tr("It appears none of the chemical potentials can be constrained"));
  //grLayout->addWidget(labelNothing);


  grLaws->setLayout(grLayout);

  QDialogButtonBox* buttonBox = new QDialogButtonBox(QDialogButtonBox::Ok
    | QDialogButtonBox::Cancel);

  connect(buttonBox, &QDialogButtonBox::accepted, this, &ConservationLawsDialog::OK);
  connect(buttonBox, &QDialogButtonBox::rejected, this, &ConservationLawsDialog::Discard);

  if (!m_parent->m_eventGeneratorMode)
    layout->addWidget(grLaws);


  QGroupBox* grMixCan = new QGroupBox(tr("Mixed-canonical ensemble"));

  QVBoxLayout* layMixCan = new QVBoxLayout();

  checkBConserve = new QCheckBox(tr("Canonical treatment of baryon number"));
  checkQConserve = new QCheckBox(tr("Canonical treatment of electric charge"));
  checkSConserve = new QCheckBox(tr("Canonical treatment of strangeness"));
  checkCConserve = new QCheckBox(tr("Canonical treatment of charm"));
  checkBConserve->setChecked(m_parent->currentConfig.CanonicalB);
  checkQConserve->setChecked(m_parent->currentConfig.CanonicalQ);
  checkSConserve->setChecked(m_parent->currentConfig.CanonicalS);
  checkCConserve->setChecked(m_parent->currentConfig.CanonicalC);

  connect(checkBConserve, &QCheckBox::toggled, this, &ConservationLawsDialog::updateControls);
  connect(checkQConserve, &QCheckBox::toggled, this, &ConservationLawsDialog::updateControls);
  connect(checkSConserve, &QCheckBox::toggled, this, &ConservationLawsDialog::updateControls);
  connect(checkCConserve, &QCheckBox::toggled, this, &ConservationLawsDialog::updateControls);

  if (m_parent->model->TPS()->hasBaryons())
    layMixCan->addWidget(checkBConserve);

  if (m_parent->model->TPS()->hasCharged())
    layMixCan->addWidget(checkQConserve);

  if (m_parent->model->TPS()->hasStrange())
    layMixCan->addWidget(checkSConserve);

  if (m_parent->model->TPS()->hasCharmed())
    layMixCan->addWidget(checkCConserve);

  grMixCan->setLayout(layMixCan);

  if (m_parent->currentConfig.Ensemble == ThermalModelConfig::EnsembleCE && layMixCan->count() > 0)
    layout->addWidget(grMixCan);

  layout->addWidget(buttonBox);

  setLayout(layout);

  setWindowTitle(tr("Thermal model configuration"));

  updateControls();
}

void ConservationLawsDialog::updateControls()
{
  if (CBmuB->isChecked()) {
    spinSBRatio->setEnabled(true);
    if (CBmuBdens->isChecked())
      CBmuBdens->setChecked(false);
  }
  else {
    spinSBRatio->setEnabled(false);
  }

  if (CBmuBdens->isChecked()) {
    spinRhoB->setEnabled(true);
  }
  else {
    spinRhoB->setEnabled(false);
  }

  if (CBmuQ->isChecked()) {
    spinQBRatio->setEnabled(true);
  }
  else {
    spinQBRatio->setEnabled(false);
  }

  bool fl = false;

  if (m_parent->model->TPS()->hasBaryons()
    && (m_parent->currentConfig.Ensemble != ThermalModelConfig::EnsembleCE ||
      !checkBConserve->isChecked())) {
    CBmuBfull->setEnabled(true);
    CBmuBdens->setEnabled(true);
    fl = true;
  }
  else {
    CBmuBfull->setEnabled(false);
    CBmuBdens->setEnabled(false);
  }

  if (m_parent->model->TPS()->hasCharged()
    && (m_parent->currentConfig.Ensemble != ThermalModelConfig::EnsembleCE ||
      !checkQConserve->isChecked())) {
    CBmuQfull->setEnabled(true);
    fl = true;
  }
  else
    CBmuQfull->setEnabled(false);

  if (m_parent->model->TPS()->hasStrange()
    && (m_parent->currentConfig.Ensemble != ThermalModelConfig::EnsembleCE ||
      !checkSConserve->isChecked())
    && m_parent->currentConfig.Ensemble != ThermalModelConfig::EnsembleSCE) {
    CBmuS->setEnabled(true);
    fl = true;
  }
  else
    CBmuS->setEnabled(false);

  if (m_parent->model->TPS()->hasCharmed()
    && (m_parent->currentConfig.Ensemble != ThermalModelConfig::EnsembleCE ||
      !checkCConserve->isChecked())
    && m_parent->currentConfig.Ensemble != ThermalModelConfig::EnsembleSCE
    && m_parent->currentConfig.Ensemble != ThermalModelConfig::EnsembleCCE) {
    CBmuC->setEnabled(true);
    fl = true;
  }
  else
    CBmuC->setEnabled(false);
}

void ConservationLawsDialog::OK()
{
  m_parent->currentConfig.ConstrainMuB = CBmuB->isChecked() || CBmuBdens->isChecked();
  m_parent->currentConfig.ConstrainMuBType = static_cast<int>(CBmuBdens->isChecked());
  m_parent->currentConfig.SoverB = spinSBRatio->value();
  m_parent->currentConfig.RhoB = spinRhoB->value();
  m_parent->currentConfig.ConstrainMuQ = CBmuQ->isChecked();
  m_parent->currentConfig.QoverB = spinQBRatio->value();
  m_parent->currentConfig.ConstrainMuS = CBmuS->isChecked();
  m_parent->currentConfig.ConstrainMuC = CBmuC->isChecked();

  m_parent->currentConfig.CanonicalB = checkBConserve->isChecked();
  m_parent->currentConfig.CanonicalQ = checkQConserve->isChecked();
  m_parent->currentConfig.CanonicalS = checkSConserve->isChecked();
  m_parent->currentConfig.CanonicalC = checkCConserve->isChecked();
  QDialog::accept();
}

void ConservationLawsDialog::toggleMuB() {
  if (CBmuBdens->isChecked() && CBmuB->isChecked()) {
    CBmuB->setChecked(false);
  } else {
    updateControls();
  }
}

OtherOptionsDialog::OtherOptionsDialog(ModelConfigWidget* parent) : QDialog(parent), m_parent(parent)
{
  QVBoxLayout* layout = new QVBoxLayout(); 
  
  QGroupBox* grReso = new QGroupBox(tr("Resonances"));

  QVBoxLayout* grLayout = new QVBoxLayout();

  CBNormBratio = new QCheckBox(tr("Renormalize branching ratios to 100%"));
  CBNormBratio->setChecked(m_parent->currentConfig.RenormalizeBR);

  CBFluctuations = new QCheckBox(tr("Compute fluctuations and correlations"));
  CBFluctuations->setChecked(m_parent->currentConfig.ComputeFluctations);

  grLayout->addWidget(CBNormBratio);
  //layout->addWidget(CBFluctuations);

  QDialogButtonBox* buttonBox = new QDialogButtonBox(QDialogButtonBox::Ok
    | QDialogButtonBox::Cancel);

  connect(buttonBox, &QDialogButtonBox::accepted, this, &OtherOptionsDialog::OK);
  connect(buttonBox, &QDialogButtonBox::rejected, this, &OtherOptionsDialog::Discard);

  QHBoxLayout* layWidths = new QHBoxLayout();
  layWidths->setAlignment(Qt::AlignLeft);

  QLabel* labelWidths = new QLabel(tr("Breit-Wigner shape:"));
  comboBWShape = new QComboBox();
  comboBWShape->addItem(tr("Relativistic"));
  comboBWShape->addItem(tr("Nonrelativistic"));
  comboBWShape->setCurrentIndex(m_parent->currentConfig.WidthShape);

  layWidths->addWidget(labelWidths);
  layWidths->addWidget(comboBWShape);

  grLayout->addLayout(layWidths);


  //grLayout->addWidget(buttonBox);

  grReso->setLayout(grLayout);


  QGroupBox* grPCE = new QGroupBox(tr("Partial chemical equilibrium (PCE)"));
  QVBoxLayout* PCELayout = new QVBoxLayout();

  CBUsePCE = new QCheckBox(tr("Apply partial chemical equilibrium"));
  CBUsePCE->setChecked(m_parent->currentConfig.UsePCE);
  connect(CBUsePCE, SIGNAL(toggled(bool)), this, SLOT(UpdateControls()));

  QHBoxLayout* layTkin= new QHBoxLayout();
  layTkin->setAlignment(Qt::AlignLeft);
  QLabel* labelTkin = new QLabel(tr("T<sub>kin</sub> (MeV):"));
  labelTkin->setToolTip(tr("Make sure T<sub>kin</sub> >= T<sub>ch</sub>. The T<sub>kin</sub> > T<sub>ch</sub> case will still compute, but may not have much physical sense."));
  spinTkin = new QDoubleSpinBox();
  spinTkin->setRange(0., 10000.);
  spinTkin->setValue(m_parent->currentConfig.Tkin * 1.e3);

  layTkin->addWidget(labelTkin);
  layTkin->addWidget(spinTkin);

  CBFreezeLongLived = new QCheckBox(tr("Freeze long-lived resonances"));
  CBFreezeLongLived->setChecked(m_parent->currentConfig.PCEFreezeLongLived);
  CBFreezeLongLived->setToolTip(tr("Yields of resonances with Γ&lt;Γ<sub>lim</sub> will be frozen in the PCE"));
   
  QHBoxLayout* layGammaLim = new QHBoxLayout();
  layGammaLim->setAlignment(Qt::AlignLeft);
  QLabel* labelGammaLim = new QLabel(tr("Γ<sub>lim</sub> (MeV):"));
  spinWidthCut = new QDoubleSpinBox();
  spinWidthCut->setRange(0., 10000.);
  spinWidthCut->setValue(m_parent->currentConfig.PCEWidthCut * 1.e3);

  layGammaLim->addWidget(labelGammaLim);
  layGammaLim->addWidget(spinWidthCut);

  CBSahaNuclei = new QCheckBox(tr("Use Saha equation for light nuclei"));
  CBSahaNuclei->setChecked(m_parent->currentConfig.PCESahaForNuclei);
  CBSahaNuclei->setToolTip(tr("Check to calculate light nuclei yields at T<sub>kin</sub> using the Saha equation. Otherwise their yields are frozen at T<sub>ch</sub>."));

  CBAnnihilation = new QCheckBox(tr("Nucleon annihilation"));
  CBAnnihilation->setChecked(m_parent->currentConfig.PCEAnnihilation);
  CBAnnihilation->setToolTip(tr("Incorporate nucleon annihilation in PCE"));

  QHBoxLayout* layPionAnnihilation = new QHBoxLayout();
  layPionAnnihilation->setAlignment(Qt::AlignLeft);
  QLabel* labelPionAnnihilation = new QLabel(tr("Average number of pions in annihilation:"));
  spinPionsAnnihilation = new QDoubleSpinBox();
  spinPionsAnnihilation->setRange(0., 10000.);
  spinPionsAnnihilation->setValue(m_parent->currentConfig.PCEPionAnnihilationNumber);
  layPionAnnihilation->addWidget(labelPionAnnihilation);
  layPionAnnihilation->addWidget(spinPionsAnnihilation);

  PCELayout->addWidget(CBUsePCE);
  if (!parent->m_thermalFitMode)
    PCELayout->addLayout(layTkin);
  PCELayout->addSpacing(10);
  PCELayout->addWidget(CBFreezeLongLived);
  PCELayout->addLayout(layGammaLim);
  PCELayout->addSpacing(10);
  PCELayout->addWidget(CBSahaNuclei);
  PCELayout->addSpacing(10);
  PCELayout->addWidget(CBAnnihilation);
  PCELayout->addLayout(layPionAnnihilation);

  grPCE->setLayout(PCELayout);

  QGroupBox* grMagneticField = new QGroupBox(tr("Magnetic field"));
  QGridLayout* BLayout = new QGridLayout();
  QLabel* labelB = new QLabel(tr("eB [GeV<sup>2</sup>]:"));
  spinB = new QDoubleSpinBox();
  spinB->setMinimum(0.);
  spinB->setMaximum(100.);
  spinB->setSingleStep(0.01);
  spinB->setDecimals(4);
  spinB->setValue(m_parent->currentConfig.MagneticFieldB);
  QLabel* labelLandauLevels = new QLabel(tr("Landau levels l<sub>max</sub>:"));
  spinLandauLevels = new QSpinBox();
  spinLandauLevels->setMinimum(0);
  spinLandauLevels->setMaximum(100000);
  spinLandauLevels->setValue(m_parent->currentConfig.MagneticFieldLmax);
  BLayout->addWidget(labelB, 0, 0);
  BLayout->addWidget(spinB, 0, 1);
  BLayout->addWidget(labelLandauLevels, 1, 0);
  BLayout->addWidget(spinLandauLevels, 1, 1);
  grMagneticField->setLayout(BLayout);

  QGroupBox* grPions = new QGroupBox(tr("Pion/kaon interactions"));
  QVBoxLayout* PionsLayout = new QVBoxLayout();
  
  CBPionInteractions = new QCheckBox(tr("Pions: Effective mass model (a la ChPT)"));
  CBPionInteractions->setChecked(m_parent->currentConfig.UseEMMPions);
  connect(CBPionInteractions, SIGNAL(toggled(bool)), this, SLOT(UpdateControls()));
  QHBoxLayout* layPionInteractions = new QHBoxLayout();
  layPionInteractions->setAlignment(Qt::AlignLeft);
  QLabel* labelPionInteractions = new QLabel(tr("f<sub>π</sub> [MeV]:"));
  spinfPi = new QDoubleSpinBox();
  spinfPi->setRange(0., 10000.);
  spinfPi->setValue(m_parent->currentConfig.EMMPionFPi * 1.e3);
  layPionInteractions->addWidget(labelPionInteractions);
  layPionInteractions->addWidget(spinfPi);

  CBKaonInteractions = new QCheckBox(tr("Kaons: Effective mass model (a la ChPT)"));
  CBKaonInteractions->setChecked(m_parent->currentConfig.UseEMMKaons);
  connect(CBKaonInteractions, SIGNAL(toggled(bool)), this, SLOT(UpdateControls()));
  QHBoxLayout* layKaonInteractions = new QHBoxLayout();
  layKaonInteractions->setAlignment(Qt::AlignLeft);
  QLabel* labelKaonInteractions = new QLabel(tr("f<sub>K</sub> [MeV]:"));
  spinfKa = new QDoubleSpinBox();
  spinfKa->setRange(0., 10000.);
  spinfKa->setValue(m_parent->currentConfig.EMMKaonFKa * 1.e3);
  layKaonInteractions->addWidget(labelKaonInteractions);
  layKaonInteractions->addWidget(spinfKa);

  PionsLayout->addWidget(CBPionInteractions);
  PionsLayout->addLayout(layPionInteractions);
  PionsLayout->addWidget(CBKaonInteractions);
  PionsLayout->addLayout(layKaonInteractions);
  grPions->setLayout(PionsLayout);

  layout->addWidget(grReso);
  if (!parent->m_eventGeneratorMode && parent->currentConfig.Ensemble != ThermalModelConfig::EnsembleGCE) {
    grPCE->setEnabled(false);
    CBUsePCE->setChecked(false);
  }
  
  layout->addWidget(grPCE);
  layout->addWidget(grMagneticField);
  layout->addWidget(grPions);

  layout->addWidget(buttonBox);

  setLayout(layout);

  UpdateControls();

  setWindowTitle(tr("Thermal model configuration"));
}

void OtherOptionsDialog::UpdateControls()
{
  bool PCEControlsEnabled = CBUsePCE->isChecked();
  spinTkin->setEnabled(PCEControlsEnabled);
  CBFreezeLongLived->setEnabled(PCEControlsEnabled);
  spinWidthCut->setEnabled(PCEControlsEnabled);
  CBSahaNuclei->setEnabled(PCEControlsEnabled);
  CBAnnihilation->setEnabled(PCEControlsEnabled);
  spinPionsAnnihilation->setEnabled(PCEControlsEnabled);

  spinfPi->setEnabled(CBPionInteractions->isChecked());
  spinfKa->setEnabled(CBKaonInteractions->isChecked());
}


void OtherOptionsDialog::OK()
{
  ThermalModelConfig& config = m_parent->currentConfig;
  config.RenormalizeBR = CBNormBratio->isChecked();
  config.WidthShape = comboBWShape->currentIndex();
  config.UsePCE = CBUsePCE->isChecked();
  config.Tkin = spinTkin->value() * 1.e-3;
  config.PCEFreezeLongLived = CBFreezeLongLived->isChecked();
  config.PCEWidthCut = spinWidthCut->value() * 1.e-3;
  config.PCESahaForNuclei = CBSahaNuclei->isChecked();
  config.PCEAnnihilation = CBAnnihilation->isChecked();
  config.PCEPionAnnihilationNumber = spinPionsAnnihilation->value();
  config.UseEMMPions = CBPionInteractions->isChecked();
  config.EMMPionFPi  = spinfPi->value() * 1.e-3;
  config.UseEMMKaons = CBKaonInteractions->isChecked();
  config.EMMKaonFKa = spinfKa->value() * 1.e-3;
  config.MagneticFieldB = spinB->value();
  config.MagneticFieldLmax = spinLandauLevels->value();
  QDialog::accept();
}
