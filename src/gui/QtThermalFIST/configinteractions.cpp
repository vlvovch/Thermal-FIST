/*
 * Thermal-FIST package
 *
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "configinteractions.h"

#include <QLayout>
#include <QLabel>
#include <QTableWidget>
#include <QHeaderView>
#include <QTimer>
#include <QGroupBox>
#include <QApplication>
#include <QMessageBox>
#include <QElapsedTimer>
#include <QFileDialog>
#include <QDialogButtonBox>
#include <QClipboard>
#include <QMenu>
#include <QKeyEvent>
#include <QDebug>

#include <cstdio>
#include <vector>
#include <fstream>
#include <sstream>

#include "HelperRoutines.h"
#include "BaseStructures.h"


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
#include "configwidgets.h"

using namespace thermalfist;

namespace {
  const QString ModelIdeal = "Ideal";
  const QString ModelDEV   = "Excluded volume (Diagonal)";
  const QString ModelCRSEV = "Excluded volume (X-terms)";
  const QString ModelQvdW  = "Quantum van der Waals";

  const QString prescrEVvdW = "van der Waals";
  const QString prescrEVCS  = "Carnahan-Starling";
  const QString prescrEVVir = "Virial";
  const QString prescrEVTVM = "Trivirial";
}

InteractionsDialog::InteractionsDialog(ModelConfigWidget* parent) : QDialog(parent), m_parent(parent)
{
  QVBoxLayout* layout = new QVBoxLayout();

  if (parent->currentConfig.InteractionModel == ThermalModelConfig::InteractionRealGas) {
    QHBoxLayout* layPrescr = new QHBoxLayout();
    comboEVprescr = new QComboBox();
    comboEVprescr->addItem(prescrEVvdW);
    comboEVprescr->addItem(prescrEVCS);
    comboEVprescr->addItem(prescrEVVir);
    comboEVprescr->addItem(prescrEVTVM);
    comboEVprescr->setCurrentIndex(parent->currentConfig.RealGasExcludedVolumePrescription);
    
    layPrescr->setAlignment(Qt::AlignLeft);
    layPrescr->addWidget(new QLabel(tr("Excluded volume prescription:")));
    layPrescr->addWidget(comboEVprescr, 0, Qt::AlignLeft);

    layout->addLayout(layPrescr);
  }
  else {
    comboEVprescr = NULL;
  }

  radSet = new QRadioButton(tr("Uniform/scaled"));

  QHBoxLayout* layB = new QHBoxLayout();
  layB->setAlignment(Qt::AlignLeft);
  QLabel* labelvdWB = new QLabel(tr("b (fm<sup>3</sup>)"));
  spinB = new QDoubleSpinBox();
  spinB->setMinimum(0.);
  spinB->setMaximum(100.);
  spinB->setDecimals(4);
  spinB->setValue(m_parent->currentConfig.vdWB);
  QLabel* labelRadius = new QLabel(tr("Radius:"));
  labelRadiusValue = new QLabel(QString::number(CuteHRGHelper::rv(spinB->value()), 'g', 4) + " fm");

  connect(spinB, SIGNAL(valueChanged(double)), this, SLOT(updateRadius()));

  layB->addWidget(labelvdWB);
  layB->addWidget(spinB);
  layB->addWidget(labelRadius);
  layB->addWidget(labelRadiusValue);

  QHBoxLayout* layA = new QHBoxLayout();
  layA->setAlignment(Qt::AlignLeft);
  QLabel* labelvdWA = new QLabel(tr("a (MeV fm<sup>3</sup>)"));
  spinA = new QDoubleSpinBox();
  spinA->setMinimum(-10000.);
  spinA->setMaximum(10000.);
  spinA->setDecimals(4);
  spinA->setValue(m_parent->currentConfig.vdWA * 1.e3);

  layA->addWidget(labelvdWA);
  layA->addWidget(spinA);

  QHBoxLayout* layScaling = new QHBoxLayout();
  layScaling->setAlignment(Qt::AlignLeft);
  QLabel* labelScaling = new QLabel(tr("EV/vdW parameters scaling:"));
  comboScaling = new QComboBox();
  comboScaling->addItem(tr("Same for all particles"));
  comboScaling->addItem(tr("Mass-proportional (bag model)"));
  comboScaling->addItem(tr("Proportional to baryon content (point-like mesons)"));

  layScaling->addWidget(labelScaling);
  layScaling->addWidget(comboScaling);

  CBMM = new QCheckBox(tr("Switch off meson-meson EV/QvdW terms"));
  CBMB = new QCheckBox(tr("Switch off meson-(anti)baryon EV/QvdW terms"));
  CBBaB = new QCheckBox(tr("Switch off baryon-antibaryon EV/QvdW terms"));
  CBBB = new QCheckBox(tr("Switch off baryon-baryon EV/QvdW terms"));

  CBMM->setChecked(m_parent->currentConfig.DisableMM);
  CBMB->setChecked(m_parent->currentConfig.DisableMB);
  CBBB->setChecked(m_parent->currentConfig.DisableBB);
  CBBaB->setChecked(m_parent->currentConfig.DisableBantiB);

  radMB = new QRadioButton(tr("Baryon/meson"));
  QGridLayout* layMB = new QGridLayout();
  QLabel* labelBB = new QLabel(tr("baryon-baryon"));
  QLabel* labelBaB = new QLabel(tr("baryon-antibaryon"));
  QLabel* labelMB = new QLabel(tr("meson-baryon"));
  QLabel* labelMM = new QLabel(tr("meson-meson"));
  layMB->addWidget(labelBB, 0, 1);
  layMB->addWidget(labelBaB, 0, 2);
  layMB->addWidget(labelMB, 0, 3);
  layMB->addWidget(labelMM, 0, 4);

  layMB->addWidget(new QLabel(tr("b (fm<sup>3</sup>)")), 1, 0);
  spinB_BB = new QDoubleSpinBox();
  spinB_BB->setMinimum(0.);
  spinB_BB->setMaximum(100.);
  spinB_BB->setDecimals(4);
  spinB_BB->setValue(m_parent->currentConfig.vdWbBB);
  spinB_BaB = new QDoubleSpinBox();
  spinB_BaB->setMinimum(0.);
  spinB_BaB->setMaximum(100.);
  spinB_BaB->setDecimals(4);
  spinB_BaB->setValue(m_parent->currentConfig.vdWbBantiB);
  spinB_MB = new QDoubleSpinBox();
  spinB_MB->setMinimum(0.);
  spinB_MB->setMaximum(100.);
  spinB_MB->setDecimals(4);
  spinB_MB->setValue(m_parent->currentConfig.vdWbMB);
  spinB_MM = new QDoubleSpinBox();
  spinB_MM->setMinimum(0.);
  spinB_MM->setMaximum(100.);
  spinB_MM->setDecimals(4);
  spinB_MM->setValue(m_parent->currentConfig.vdWbMM);
  layMB->addWidget(spinB_BB, 1, 1);
  layMB->addWidget(spinB_BaB,1, 2);
  layMB->addWidget(spinB_MB, 1, 3);
  layMB->addWidget(spinB_MM, 1, 4);

  layMB->addWidget(new QLabel(tr("a (MeV fm<sup>3</sup>)")), 2, 0);
  spinA_BB = new QDoubleSpinBox();
  spinA_BB->setMinimum(-10000.);
  spinA_BB->setMaximum(10000.);
  spinA_BB->setDecimals(4);
  spinA_BB->setValue(m_parent->currentConfig.vdWaBB * 1.e3);
  spinA_BaB = new QDoubleSpinBox();
  spinA_BaB->setMinimum(-10000.);
  spinA_BaB->setMaximum(10000.);
  spinA_BaB->setDecimals(4);
  spinA_BaB->setValue(m_parent->currentConfig.vdWaBantiB * 1.e3);
  spinA_MB = new QDoubleSpinBox();
  spinA_MB->setMinimum(-10000.);
  spinA_MB->setMaximum(10000.);
  spinA_MB->setDecimals(4);
  spinA_MB->setValue(m_parent->currentConfig.vdWaMB * 1.e3);
  spinA_MM = new QDoubleSpinBox();
  spinA_MM->setMinimum(-10000.);
  spinA_MM->setMaximum(10000.);
  spinA_MM->setDecimals(4);
  spinA_MM->setValue(m_parent->currentConfig.vdWaMM * 1.e3);
  layMB->addWidget(spinA_BB, 2, 1);
  layMB->addWidget(spinA_BaB, 2, 2);
  layMB->addWidget(spinA_MB, 2, 3);
  layMB->addWidget(spinA_MM, 2, 4);

  radLoad = new QRadioButton(tr("Parameters read from external file"));


  QHBoxLayout* layFile = new QHBoxLayout();
  layFile->setAlignment(Qt::AlignLeft);
  leFilePath = new QLineEdit("");
  leFilePath->setReadOnly(true);
  leFilePath->setText(QString::fromStdString(m_parent->currentConfig.InteractionInput));
  leFilePath->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Minimum);
  buttonChooseFile = new QPushButton(tr("Choose file..."));
  connect(buttonChooseFile, SIGNAL(clicked()), this, SLOT(chooseInputFile()));
  layFile->addWidget(leFilePath);
  layFile->addWidget(buttonChooseFile);

  layout->addWidget(radSet, 0, Qt::AlignLeft);
  layout->addLayout(layB);
  layout->addLayout(layA);
  layout->addLayout(layScaling);
  layout->addWidget(CBMM, 0, Qt::AlignLeft);
  layout->addWidget(CBMB, 0, Qt::AlignLeft);
  layout->addWidget(CBBaB, 0, Qt::AlignLeft);
  layout->addWidget(CBBB, 0, Qt::AlignLeft);
  layout->addSpacing(30);
  layout->addWidget(radMB, 0, Qt::AlignLeft);
  layout->addLayout(layMB);
  layout->addSpacing(30);
  layout->addWidget(radLoad, 0, Qt::AlignLeft);
  layout->addLayout(layFile);


  // QGroupBox* groupMultipleSolutions = new QGroupBox(tr("Multiple solutions"));
  QVBoxLayout* layoutMultipleSolutions = new QVBoxLayout();
  CBSearchMultipleSolutions = new QCheckBox(tr("Search for multiple solutions (in case of phase transition)"));
  CBSearchMultipleSolutions->setChecked(m_parent->currentConfig.SearchMultipleSolutions);
  layoutMultipleSolutions->addWidget(CBSearchMultipleSolutions);
  // groupMultipleSolutions->setLayout(layoutMultipleSolutions);
  layout->addLayout(layoutMultipleSolutions);
  // layout->addWidget(groupMultipleSolutions);


  groupMC = new QGroupBox(tr("Event generator options"));
  QVBoxLayout* layoutMC = new QVBoxLayout();
  CBEVMult = new QCheckBox(tr("Use rejection sampling for excluded volume multiplicities"));
  CBEVMult->setChecked(parent->currentConfig.fUseEVRejectionMultiplicity);
  CBEVCoord = new QCheckBox(tr("Use rejection sampling for excluded volume in coordinate space"));
  CBEVCoord->setChecked(parent->currentConfig.fUseEVRejectionCoordinates);
  CBEVSPR = new QCheckBox(tr("Apply the SPR approximation"));
  CBEVSPR->setChecked(parent->currentConfig.fUseEVUseSPRApproximation);
  connect(CBEVCoord, SIGNAL(toggled(bool)), this, SLOT(updateSPR()));
  layoutMC->addWidget(CBEVMult);
  layoutMC->addWidget(CBEVCoord);
  layoutMC->addWidget(CBEVSPR);
  groupMC->setLayout(layoutMC);
  layout->addWidget(groupMC);

  QDialogButtonBox* buttonBox = new QDialogButtonBox(QDialogButtonBox::Ok
    | QDialogButtonBox::Cancel);

  connect(buttonBox, &QDialogButtonBox::accepted, this, &InteractionsDialog::OK);
  connect(buttonBox, &QDialogButtonBox::rejected, this, &InteractionsDialog::Discard);

  layout->addWidget(buttonBox);

  setLayout(layout);

  setWindowTitle(tr("Excluded volume/van der Waals parameters"));


  if (m_parent->currentConfig.InteractionScaling == 4) {
    radMB->setChecked(true);
  }
  else if (m_parent->currentConfig.InteractionScaling < 3) {
    radSet->setChecked(true);
    comboScaling->setCurrentIndex(m_parent->currentConfig.InteractionScaling);
  }
  else
    radLoad->setChecked(true);

  connect(radSet, &QRadioButton::toggled, this, &InteractionsDialog::modeToggled);
  connect(radMB, &QRadioButton::toggled, this, &InteractionsDialog::modeToggled);
  modeToggled();
}

void InteractionsDialog::modeToggled()
{
  CBMM->setEnabled(false);
  CBMB->setEnabled(false);
  CBBaB->setEnabled(false);
  CBBB->setEnabled(false);

  if (radSet->isChecked()) {
    leFilePath->setEnabled(false);
    buttonChooseFile->setEnabled(false);

    spinB_BB->setEnabled(false);
    spinB_BaB->setEnabled(false);
    spinB_MB->setEnabled(false);
    spinB_MM->setEnabled(false);
    spinA_BB->setEnabled(false);
    spinA_BaB->setEnabled(false);
    spinA_MB->setEnabled(false);
    spinA_MM->setEnabled(false);

    spinB->setEnabled(true);
    spinA->setEnabled(true);
    comboScaling->setEnabled(true);
    if (m_parent->currentConfig.InteractionModel == ThermalModelConfig::InteractionQVDW || 
      m_parent->currentConfig.InteractionModel == ThermalModelConfig::InteractionRealGas) {
      spinA->setEnabled(true);
    }
    else {
      spinA->setEnabled(false);
    }
    if (m_parent->currentConfig.InteractionModel == ThermalModelConfig::InteractionEVCrossterms
      || m_parent->currentConfig.InteractionModel == ThermalModelConfig::InteractionQVDW ||
      m_parent->currentConfig.InteractionModel == ThermalModelConfig::InteractionRealGas)
    {
      CBMM->setEnabled(true);
      CBMB->setEnabled(true);
      CBBaB->setEnabled(true);
      CBBB->setEnabled(true);
    }
  }
  else if (radMB->isChecked()) {
    leFilePath->setEnabled(false);
    buttonChooseFile->setEnabled(false);

    spinB_BB->setEnabled(true);
    spinB_BaB->setEnabled(true);
    spinB_MB->setEnabled(true);
    spinB_MM->setEnabled(true);
    spinA_BB->setEnabled(true);
    spinA_BaB->setEnabled(true);
    spinA_MB->setEnabled(true);
    spinA_MM->setEnabled(true);

    if (m_parent->currentConfig.InteractionModel == ThermalModelConfig::InteractionEVCrossterms) {
      spinA_BB->setEnabled(false);
      spinA_BaB->setEnabled(false);
      spinA_MB->setEnabled(false);
      spinA_MM->setEnabled(false);
    }

    spinB->setEnabled(false);
    spinA->setEnabled(false);
    comboScaling->setEnabled(false);

    CBMM->setEnabled(false);
    CBMB->setEnabled(false);
    CBBaB->setEnabled(false);
    CBBB->setEnabled(false);
  }
  else {
    leFilePath->setEnabled(true);
    buttonChooseFile->setEnabled(true);

    spinB->setEnabled(false);
    spinA->setEnabled(false);
    comboScaling->setEnabled(false);

    CBMM->setEnabled(false);
    CBMB->setEnabled(false);
    CBBaB->setEnabled(false);
    CBBB->setEnabled(false);


    spinB_BB->setEnabled(false);
    spinB_BaB->setEnabled(false);
    spinB_MB->setEnabled(false);
    spinB_MM->setEnabled(false);
    spinA_BB->setEnabled(false);
    spinA_BaB->setEnabled(false);
    spinA_MB->setEnabled(false);
    spinA_MM->setEnabled(false);

    if (leFilePath->text() == "") {
      //  chooseInputFile();
    }
  }

  if (!m_parent->m_eventGeneratorMode || m_parent->currentConfig.InteractionModel == ThermalModelConfig::InteractionIdeal)
  {
    groupMC->setVisible(false);
  }
  else {
    groupMC->setVisible(true);
    CBEVMult->setVisible(true);
    CBEVCoord->setVisible(true);
    CBEVSPR->setVisible(true);
  }

  if (m_parent->currentConfig.InteractionModel == ThermalModelConfig::InteractionQVDW ||
    m_parent->currentConfig.InteractionModel == ThermalModelConfig::InteractionRealGas) {
    CBEVMult->setVisible(false);
  }

  updateSPR();
}


void InteractionsDialog::chooseInputFile()
{
  QString listpathprefix = QString(ThermalFIST_INPUT_FOLDER) + "/interaction";
  if (leFilePath->text().size() != 0)
    listpathprefix = QString(leFilePath->text());
  QString path = QFileDialog::getOpenFileName(this, tr("Open file with EV/vdW parameters"), listpathprefix);
  if (path.length() > 0)
  {
    leFilePath->setText(path);
  }
}

void InteractionsDialog::updateRadius()
{
  labelRadiusValue->setText(QString::number(CuteHRGHelper::rv(spinB->value()), 'g', 4) + " fm");
}

void InteractionsDialog::updateSPR()
{
  CBEVSPR->setEnabled(CBEVCoord->isChecked());
}


void InteractionsDialog::OK()
{
  if (radSet->isChecked())
    m_parent->currentConfig.InteractionScaling = comboScaling->currentIndex();
  else if (radMB->isChecked()) {
    m_parent->currentConfig.InteractionScaling = 4;
  }
  else {
    m_parent->currentConfig.InteractionScaling = 3;
  }

  m_parent->currentConfig.vdWA = spinA->value() * 1.e-3;
  m_parent->currentConfig.vdWB = spinB->value();
  m_parent->currentConfig.InteractionInput = leFilePath->text().toStdString();

  m_parent->currentConfig.DisableMM = CBMM->isChecked();
  m_parent->currentConfig.DisableMB = CBMB->isChecked();
  m_parent->currentConfig.DisableBB = CBBB->isChecked();
  m_parent->currentConfig.DisableBantiB = CBBaB->isChecked();

  m_parent->currentConfig.vdWaBB = spinA_BB->value() * 1.e-3;
  m_parent->currentConfig.vdWaBantiB = spinA_BaB->value() * 1.e-3;
  m_parent->currentConfig.vdWaMB = spinA_MB->value() * 1.e-3;
  m_parent->currentConfig.vdWaMM = spinA_MM->value() * 1.e-3;
  m_parent->currentConfig.vdWbBB = spinB_BB->value();
  m_parent->currentConfig.vdWbBantiB = spinB_BaB->value();
  m_parent->currentConfig.vdWbMB = spinB_MB->value();
  m_parent->currentConfig.vdWbMM = spinB_MM->value();

  m_parent->currentConfig.SearchMultipleSolutions = CBSearchMultipleSolutions->isChecked();

  m_parent->currentConfig.fUseEVRejectionMultiplicity = CBEVMult->isChecked();
  m_parent->currentConfig.fUseEVRejectionCoordinates = CBEVCoord->isChecked();
  m_parent->currentConfig.fUseEVUseSPRApproximation = CBEVSPR->isChecked();

  m_parent->currentConfig.vdWparams = QvdWParameters::GetParameters(m_parent->model->TPS(), &m_parent->currentConfig);

  if (comboEVprescr != NULL)
    m_parent->currentConfig.RealGasExcludedVolumePrescription = comboEVprescr->currentIndex();

  QDialog::accept();
}

QvdWParameters QvdWParameters::GetParameters(ThermalParticleSystem* TPS, const ThermalModelConfig* config)
{
  QvdWParameters ret;
  int N = TPS->ComponentsNumber();
  ret.m_aij = ret.m_bij = 
    std::vector<std::vector<double>>(N, std::vector<double>(N, 0.));

  if (config->InteractionScaling == 4) {
    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < N; ++j) {
        int B1 = TPS->Particle(i).BaryonCharge();
        int B2 = TPS->Particle(j).BaryonCharge();

        if (B1 * B2 > 0) {
          ret.m_bij[i][j] = config->vdWbBB;
          ret.m_aij[i][j] = config->vdWaBB;
        } 
        else if (B1 * B2 < 0) {
          ret.m_bij[i][j] = config->vdWbBantiB;
          ret.m_aij[i][j] = config->vdWaBantiB;
        }
        else if (B1 == 0 && B2 == 0) {
          ret.m_bij[i][j] = config->vdWbMM;
          ret.m_aij[i][j] = config->vdWaMM;
        }
        else {
          ret.m_bij[i][j] = config->vdWbMB;
          ret.m_aij[i][j] = config->vdWaMB;
        }
      }
    }
    return ret;
  }

  if (config->InteractionScaling != 3) {
    // First repulsion
    double radius = CuteHRGHelper::rv(config->vdWB);
    std::vector<double> radii(N, 0.);
    // Uniform EV
    if (config->InteractionScaling == 0) {
      std::fill(radii.begin(), radii.end(), radius);
    }

    // Bag model EV, nucleon mass here equals 0.938 GeV
    if (config->InteractionScaling == 1) {
      for (int i = 0; i < TPS->Particles().size(); ++i) {
        ThermalParticle& part = TPS->Particle(i);
        radii[i] = radius * pow(part.Mass() / xMath::mnucleon(), 1. / 3.);
      }
    }

    // Linear in baryon content
    if (config->InteractionScaling == 2) {
      for (int i = 0; i < TPS->Particles().size(); ++i) {
        ThermalParticle& part = TPS->Particle(i);
        if (part.BaryonCharge() != 0)
          radii[i] = radius * pow(abs(part.BaryonCharge()), 1. / 3.);
        else
          radii[i] = 0.;
      }
    }

    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < N; ++j) {
        int B1 = TPS->Particle(i).BaryonCharge();
        int B2 = TPS->Particle(j).BaryonCharge();

        bool fl = true;
        if (B1 * B2 > 0)
          fl &= !config->DisableBB;
        else if (B1 * B2 < 0)
          fl &= !config->DisableBantiB;
        else if (B1 == 0 && B2 == 0)
          fl &= !config->DisableMM;
        else
          fl &= !config->DisableMB;

        if (config->InteractionModel == ThermalModelConfig::InteractionEVDiagonal)
          fl = true;

        if (!fl) {
          ret.m_bij[i][j] = 0.;
        }
        else {
          ret.m_bij[i][j] = CuteHRGHelper::btilrr(radii[i], radii[j]);
        }

      }
    }


    // Now the attraction
    //if (config.InteractionModel == ThermalModelConfig::InteractionQVDW) {
      std::vector<double> as(TPS->Particles().size(), config->vdWA);

      // Mass-proportional, nucleon mass here equals 0.938 GeV
      if (config->InteractionScaling == 1) {
        for (int i = 0; i < TPS->Particles().size(); ++i) {
          ThermalParticle& part = TPS->Particle(i);
          as[i] = config->vdWA * part.Mass() / xMath::mnucleon();
        }
      }

      // Linear in baryon content
      if (config->InteractionScaling == 2) {
        for (int i = 0; i < TPS->Particles().size(); ++i) {
          ThermalParticle& part = TPS->Particle(i);
          if (part.BaryonCharge() != 0)
            as[i] = config->vdWA * abs(part.BaryonCharge());
          else
            as[i] = 0.;
        }
      }

      for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
          int B1 = TPS->Particle(i).BaryonCharge();
          int B2 = TPS->Particle(j).BaryonCharge();

          bool fl = true;
          if (B1 * B2 > 0)
            fl &= !config->DisableBB;
          else if (B1 * B2 < 0)
            fl &= !config->DisableBantiB;
          else if (B1 == 0 && B2 == 0)
            fl &= !config->DisableMM;
          else
            fl &= !config->DisableMB;

          if (!fl) {
            ret.m_aij[i][j] = 0.;
          }
          else {
            // Fill aij assuming "chemistry rule" aij = \sqrt{ai * aj}
            ret.m_aij[i][j] = sqrt(as[i] * as[j]);
          }
        }
      }
    //}
      return ret;
  }

  // Read from file
  if (config->InteractionScaling == 3) {
    std::ifstream fin(config->InteractionInput.c_str());
    char cc[2000];
    while (!fin.eof()) {
      fin.getline(cc, 2000);
      std::string tmp = std::string(cc);
      std::vector<std::string> elems = CuteHRGHelper::split(tmp, '#');
      if (elems.size() < 1)
        continue;

      std::istringstream iss(elems[0]);
      
      if (config->ModelType == ThermalModelConfig::DiagonalEV) {
        long long pdgid;
        double b;
        if (iss >> pdgid >> b) {
          int ind = TPS->PdgToId(pdgid);
          if (ind != -1)
            std::fill(ret.m_bij[ind].begin(), ret.m_bij[ind].end(), b);
        }
      }
      else if (config->ModelType == ThermalModelConfig::CrosstermsEV || config->ModelType == ThermalModelConfig::QvdW ||
        config->ModelType == ThermalModelConfig::InteractionRealGas) {
        long long pdgid1, pdgid2;
        double b, a;
        if (iss >> pdgid1 >> pdgid2 >> b) {
          if (!(iss >> a))
            a = 0.;
          int ind1 = TPS->PdgToId(pdgid1);
          int ind2 = TPS->PdgToId(pdgid2);
          if (ind1 != -1 && ind2 != -1) {
            ret.m_bij[ind1][ind2] = b;
            ret.m_aij[ind1][ind2] = a;
          }
        }
      }
    }
    fin.close();

    return ret;
  }

  return ret;
}

QvdWParametersTableDialog::QvdWParametersTableDialog(
  QWidget* parent,
  ThermalModelConfig* iconfig,
  thermalfist::ThermalParticleSystem* iTPS) :
  QDialog(parent), TPS(iTPS), config(iconfig)
{
  QVBoxLayout* layout = new QVBoxLayout();

  QHBoxLayout* laySel = new QHBoxLayout();
  laySel->setAlignment(Qt::AlignLeft);

  QLabel* labelType = new QLabel(tr("Parameters:"));
  comboType = new QComboBox();
  //QString bname = tr("b<sup>ij</sup> [fm<sup>3</sup>]");
  QString bname = tr("b_ij [fm^3]");
  if (config->ModelType == ThermalModelConfig::DiagonalEV) {
    bname = tr("v_i [fm^3]");
  }
  comboType->addItem(bname);
  if (config->ModelType == ThermalModelConfig::QvdW || config->ModelType == ThermalModelConfig::RealGas) {
    comboType->addItem(tr("a_ij [MeV fm^3]"));
  }
  comboType->setCurrentIndex(0);
  connect(comboType, SIGNAL(currentIndexChanged(int)), this, SLOT(recalculate()));

  laySel->addWidget(labelType);
  laySel->addWidget(comboType);

  tablevdW = new QTableWidget();
  tablevdW->installEventFilter(this);
  tablevdW->setEditTriggers(QAbstractItemView::NoEditTriggers);
  tablevdW->setHorizontalScrollBarPolicy(Qt::ScrollBarAsNeeded);

  layout->addLayout(laySel);
  layout->addWidget(tablevdW);

  setLayout(layout);

  this->setWindowTitle(tr("van der Waals parameters"));

  recalculate();
  //QTimer::singleShot(0, this, SLOT(checkFixTableSize()));
}

void QvdWParametersTableDialog::recalculate()
{
  tablevdW->clear();

  QVector<long long> idtopdg(0);
  for (int i = 0; i < TPS->ComponentsNumber(); ++i) {
    const ThermalParticle& part = TPS->Particle(i);
    idtopdg.push_back(part.PdgId());
  }

  tablevdW->setRowCount(idtopdg.size());
  if (config->ModelType == ThermalModelConfig::DiagonalEV) {
    tablevdW->setColumnCount(1);
  }
  else {
    tablevdW->setColumnCount(idtopdg.size());
  }

  for (int i = 0; i < idtopdg.size(); ++i) {
    QString name = TPS->ParticleByPDG(idtopdg[i]).Name().c_str();
    tablevdW->setHorizontalHeaderItem(i, new QTableWidgetItem(name));
    if (!(config->ModelType == ThermalModelConfig::DiagonalEV)) {
      tablevdW->setVerticalHeaderItem(i, new QTableWidgetItem(name));
    }
  }

  for (int i = 0; i < TPS->ComponentsNumber(); ++i) {
    for (int j = 0; j < TPS->ComponentsNumber(); ++j) {
      double val = config->vdWparams.m_bij[i][j];
      if (comboType->currentIndex() == 1)
        val = config->vdWparams.m_aij[i][j] * 1.e3;

      tablevdW->setItem(i, j, new QTableWidgetItem(QString::number(val)));

      if ((config->ModelType == ThermalModelConfig::DiagonalEV))
        break;
    }
  }
}

bool QvdWParametersTableDialog::eventFilter(QObject* obj, QEvent* event)
{
  if (obj == tablevdW) {
    if (event->type() == QEvent::KeyPress && static_cast<QKeyEvent*>(event)->matches(QKeySequence::Copy)) {
      copyTableViewSelectionToClipBoard(tablevdW);
      return true;
    }
    else {
      //return false;
      return QDialog::eventFilter(obj, event);
    }
  }
  else {
    // pass the event on to the parent class
    return QDialog::eventFilter(obj, event);
  }
}
