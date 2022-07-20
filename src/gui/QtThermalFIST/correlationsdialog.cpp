/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "correlationsdialog.h"

#include <QLayout>
#include <QHeaderView>
#include <QTimer>
#include <QLabel>
#include <QApplication>
#include <QClipboard>
#include <QMenu>
#include <QKeyEvent>
#include <QDebug>

#include "HelperRoutines.h"

#include "HRGBase/ThermalModelBase.h"

using namespace thermalfist;

CorrelationsDialog::CorrelationsDialog(QWidget* parent, ThermalModelBase* mod) :
  QDialog(parent), model(mod)
{
  QVBoxLayout* layout = new QVBoxLayout();

  QHBoxLayout* laySel = new QHBoxLayout();
  laySel->setAlignment(Qt::AlignLeft);

  QLabel* labelPart = new QLabel(tr("Feeddown:"));
  comboFeeddown = new QComboBox();
  comboFeeddown->addItem(tr("Primordial"));
  comboFeeddown->addItem(tr("Final"));
  comboFeeddown->setCurrentIndex(1);
  connect(comboFeeddown, SIGNAL(currentIndexChanged(int)), this, SLOT(recalculate()));

  QLabel* labelType = new QLabel(tr("Type:"));
  comboType = new QComboBox();
  comboType->addItem(tr("Particle"));
  comboType->addItem(tr("Net-particle"));
  comboType->addItem(tr("Conserved charges"));
  comboType->addItem(tr("Particle/Conserved-charge"));
  comboType->addItem(tr("Net-particle/Conserved-charge"));
  comboType->setCurrentIndex(0);
  connect(comboType, SIGNAL(currentIndexChanged(int)), this, SLOT(recalculate()));

  QLabel* labelQuantity = new QLabel(tr("Quantity:"));
  comboQuantity = new QComboBox();
  comboQuantity->addItem(tr("Susceptibility"));
  comboQuantity->addItem(tr("Moment"));
  comboQuantity->addItem(tr("Scaled moment"));
  comboQuantity->addItem(tr("Pearson correlation"));
  comboQuantity->addItem(tr("Delta[N1,N2]"));
  comboQuantity->addItem(tr("Sigma[N1,N2]"));
  comboQuantity->setCurrentIndex(0);
  connect(comboQuantity, SIGNAL(currentIndexChanged(int)), this, SLOT(recalculate()));

  laySel->addWidget(labelPart);
  laySel->addWidget(comboFeeddown);
  laySel->addSpacing(20);
  laySel->addWidget(labelType);
  laySel->addWidget(comboType);
  laySel->addSpacing(20);
  laySel->addWidget(labelQuantity);
  laySel->addWidget(comboQuantity);

  tableCorr = new QTableWidget();
  tableCorr->installEventFilter(this);
  tableCorr->setEditTriggers(QAbstractItemView::NoEditTriggers);
  tableCorr->setHorizontalScrollBarPolicy(Qt::ScrollBarAsNeeded);

  layout->addLayout(laySel);
  layout->addWidget(tableCorr);


  setLayout(layout);

  this->setWindowTitle(tr("Correlations"));

  recalculate();
  //QTimer::singleShot(0, this, SLOT(checkFixTableSize()));
}

void CorrelationsDialog::checkFixTableSize() {
  int tle = 0;
  for (int i = 0; i < tableCorr->horizontalHeader()->count(); ++i)
    tle += tableCorr->horizontalHeader()->sectionSize(i);
  tableCorr->setMinimumWidth(tle + 23);
  tle = 0;
  for (int i = 0; i < tableCorr->horizontalHeader()->count(); ++i)
    tle += tableCorr->horizontalHeader()->sectionSize(i);
  tableCorr->setMinimumWidth(tle + 23);
}

void CorrelationsDialog::recalculate()
{
  tableCorr->clear();

  if (comboType->currentIndex() <= 1) {
    QVector<long long> idtopdg(0);
    for (int i = 0; i < model->TPS()->ComponentsNumber(); ++i) {
      const ThermalParticle& part = model->TPS()->Particle(i);
      if (comboFeeddown->currentIndex() == 0 || (comboFeeddown->currentIndex() == 1 && part.IsStable())) {
        if (comboType->currentIndex() == 0 || (comboType->currentIndex() == 1 && !part.IsNeutral() && part.PdgId() > 0)) {
          idtopdg.push_back(part.PdgId());
        }
      }
    }

    tableCorr->setColumnCount(idtopdg.size());
    tableCorr->setRowCount(idtopdg.size());

    for (int i = 0; i < idtopdg.size(); ++i) {
      QString name = model->TPS()->ParticleByPDG(idtopdg[i]).Name().c_str();
      if (comboType->currentIndex() == 1)
        name = "net-" + name;
      tableCorr->setVerticalHeaderItem(i, new QTableWidgetItem(name));
      tableCorr->setHorizontalHeaderItem(i, new QTableWidgetItem(name));
    }

    for (int i = 0; i < idtopdg.size(); ++i) {
      for (int j = 0; j < idtopdg.size(); ++j) {
        long long pdg1 = idtopdg[i], pdg2 = idtopdg[j];

        double corr = 1.;

        if (comboType->currentIndex() == 0) {
          if (comboFeeddown->currentIndex() == 0)
            corr = model->TwoParticleSusceptibilityPrimordialByPdg(pdg1, pdg2);
          else
            corr = model->TwoParticleSusceptibilityFinalByPdg(pdg1, pdg2);
        }
        else {
          if (comboFeeddown->currentIndex() == 0)
            corr = model->TwoParticleSusceptibilityPrimordialByPdg(pdg1, pdg2)
            - model->TwoParticleSusceptibilityPrimordialByPdg(pdg1, -pdg2)
            - model->TwoParticleSusceptibilityPrimordialByPdg(-pdg1, pdg2)
            + model->TwoParticleSusceptibilityPrimordialByPdg(-pdg1, -pdg2);
          else
            corr = model->TwoParticleSusceptibilityFinalByPdg(pdg1, pdg2)
            - model->TwoParticleSusceptibilityFinalByPdg(pdg1, -pdg2)
            - model->TwoParticleSusceptibilityFinalByPdg(-pdg1, pdg2)
            + model->TwoParticleSusceptibilityFinalByPdg(-pdg1, -pdg2);
        }

        if (comboQuantity->currentIndex() == 0)
          tableCorr->setItem(i, j, new QTableWidgetItem(QString::number(corr)));
        else if (comboQuantity->currentIndex() == 1)
          tableCorr->setItem(i, j, new QTableWidgetItem(QString::number(corr * model->Volume() * pow(model->Parameters().T, 3) * pow(xMath::GeVtoifm(), 3))));
        //else if (comboQuantity->currentIndex() == 2) {
        //  double val = corr * model->Volume() * pow(model->Parameters().T, 3) * pow(xMath::GeVtoifm(), 3);
        //  double yld1 = model->GetYield(pdg1, Feeddown::Primordial);
        //  double yld2 = model->GetYield(pdg2, Feeddown::Primordial);
        //  if (comboFeeddown->currentIndex() == 1) {
        //    yld1 = model->GetYield(pdg1, Feeddown::StabilityFlag);
        //    yld2 = model->GetYield(pdg2, Feeddown::StabilityFlag);
        //  }
        //  tableCorr->setItem(i, j, new QTableWidgetItem(QString::number(val / sqrt(yld1 * yld2))));
        //}
        //else if (comboQuantity->currentIndex() == 3) {
        //  double val = corr;
        //  double var1 = model->TwoParticleSusceptibilityPrimordialByPdg(pdg1, pdg1);
        //  double var2 = model->TwoParticleSusceptibilityPrimordialByPdg(pdg2, pdg2);
        //  if (comboFeeddown->currentIndex() == 1) {
        //    var1 = model->TwoParticleSusceptibilityFinalByPdg(pdg1, pdg1);
        //    var2 = model->TwoParticleSusceptibilityFinalByPdg(pdg2, pdg2);
        //  }
        //  tableCorr->setItem(i, j, new QTableWidgetItem(QString::number(val / sqrt(var1 * var2))));
        //}
        else {

          double N1 = 0., N2 = 0.;
          double wn1 = 0., wn2 = 0.;

          if (comboType->currentIndex() == 0) {
            if (comboFeeddown->currentIndex() == 0) {
              N1 = model->GetDensity(pdg1, Feeddown::Primordial) * model->Volume();
              N2 = model->GetDensity(pdg2, Feeddown::Primordial) * model->Volume();
              wn1 = model->ScaledVariancePrimordial(model->TPS()->PdgToId(pdg1));
              wn2 = model->ScaledVariancePrimordial(model->TPS()->PdgToId(pdg2));
            }
            else {
              N1 = model->GetDensity(pdg1, Feeddown::StabilityFlag) * model->Volume();
              N2 = model->GetDensity(pdg2, Feeddown::StabilityFlag) * model->Volume();
              wn1 = model->ScaledVarianceTotal(model->TPS()->PdgToId(pdg1));
              wn2 = model->ScaledVarianceTotal(model->TPS()->PdgToId(pdg2));
            }
          }
          else {
            if (comboFeeddown->currentIndex() == 0) {
              N1 = model->GetDensity(pdg1, Feeddown::Primordial) * model->Volume()
                - model->GetDensity(-pdg1, Feeddown::Primordial) * model->Volume();
              N2 = model->GetDensity(pdg2, Feeddown::Primordial) * model->Volume()
                - model->GetDensity(-pdg2, Feeddown::Primordial) * model->Volume();
              wn1 = (model->TwoParticleSusceptibilityPrimordialByPdg(pdg1, pdg1)
                - model->TwoParticleSusceptibilityPrimordialByPdg(pdg1, -pdg1)
                - model->TwoParticleSusceptibilityPrimordialByPdg(-pdg1, pdg1)
                + model->TwoParticleSusceptibilityPrimordialByPdg(-pdg1, -pdg1)) / N1;
              wn1 *= model->Volume() * pow(model->Parameters().T, 3) * pow(xMath::GeVtoifm(), 3);
              wn2 = (model->TwoParticleSusceptibilityPrimordialByPdg(pdg2, pdg2)
                - model->TwoParticleSusceptibilityPrimordialByPdg(pdg2, -pdg2)
                - model->TwoParticleSusceptibilityPrimordialByPdg(-pdg2, pdg2)
                + model->TwoParticleSusceptibilityPrimordialByPdg(-pdg2, -pdg2)) / N2;
              wn2 *= model->Volume() * pow(model->Parameters().T, 3) * pow(xMath::GeVtoifm(), 3);
            }
            else {
              N1 = model->GetDensity(pdg1, Feeddown::StabilityFlag) * model->Volume()
                - model->GetDensity(-pdg1, Feeddown::StabilityFlag) * model->Volume();
              N2 = model->GetDensity(pdg2, Feeddown::StabilityFlag) * model->Volume()
                - model->GetDensity(-pdg2, Feeddown::StabilityFlag) * model->Volume();
              wn1 = (model->TwoParticleSusceptibilityFinalByPdg(pdg1, pdg1)
                - model->TwoParticleSusceptibilityFinalByPdg(pdg1, -pdg1)
                - model->TwoParticleSusceptibilityFinalByPdg(-pdg1, pdg1)
                + model->TwoParticleSusceptibilityFinalByPdg(-pdg1, -pdg1)) / N1;
              wn1 *= model->Volume() * pow(model->Parameters().T, 3) * pow(xMath::GeVtoifm(), 3);
              wn2 = (model->TwoParticleSusceptibilityFinalByPdg(pdg2, pdg2)
                - model->TwoParticleSusceptibilityFinalByPdg(pdg2, -pdg2)
                - model->TwoParticleSusceptibilityFinalByPdg(-pdg2, pdg2)
                + model->TwoParticleSusceptibilityFinalByPdg(-pdg2, -pdg2)) / N2;
              wn2 *= model->Volume() * pow(model->Parameters().T, 3) * pow(xMath::GeVtoifm(), 3);
            }
          }

          corr *= model->Volume() * pow(model->Parameters().T, 3) * pow(xMath::GeVtoifm(), 3);

          // Scaled moment
          if (comboQuantity->currentIndex() == 2) {
            if (pdg1 == pdg2) {
              qDebug() << pdg1 << " " << corr << " " << N1 << " " << N2 << endl;
            }
            if (i != j) {
              tableCorr->setItem(i, j, new QTableWidgetItem("N/A"));
            }
            else {
              tableCorr->setItem(i, j, new QTableWidgetItem(QString::number(corr / N1)));
            }
          }
          // Pearson
          else if (comboQuantity->currentIndex() == 3) {
            double var1 = wn1 * N1, var2 = wn2 * N2;
            tableCorr->setItem(i, j, new QTableWidgetItem(QString::number(corr / sqrt(var1 * var2))));
          }
          else {
            if (pdg1 == pdg2) {
              tableCorr->setItem(i, j, new QTableWidgetItem(QString::number(0.)));
              continue;
            }

            // Delta
            if (comboQuantity->currentIndex() == 4) {
              double DeltaN1N2 = -(N1 * wn2 - N2 * wn1) / (N2 - N1);
              tableCorr->setItem(i, j, new QTableWidgetItem(QString::number(DeltaN1N2)));
              continue;
            }
            // Sigma
            else {
              double SigmaN1N2 = (N1 * wn2 + N2 * wn1 - 2. * corr) / (N2 + N1);
              tableCorr->setItem(i, j, new QTableWidgetItem(QString::number(SigmaN1N2)));
              continue;
            }
          }
        }
      }
    }
  }

  // Conserved charges
  if (comboType->currentIndex() == 2) {
    std::vector<int> flags(4, 1);
    flags[0] = model->TPS()->hasBaryons();
    flags[1] = model->TPS()->hasCharged();
    flags[2] = model->TPS()->hasStrange();
    flags[3] = model->TPS()->hasCharmed();
    {
      int totnum = 0;
      for (int i = 0; i < flags.size(); ++i)
        totnum += flags[i];
      tableCorr->setColumnCount(totnum);
      tableCorr->setRowCount(totnum);
    }
    std::vector<QString> names = { "B", "Q", "S", "C" };
    int i1 = 0, i2 = 0;
    for (int i = 0; i < flags.size(); ++i) {
      if (flags[i]) {
        tableCorr->setHorizontalHeaderItem(i1, new QTableWidgetItem(names[i1]));
        i2 = 0;
        for (int j = 0; j < flags.size(); ++j) {
          if (flags[j]) {
            if (i1 == 0)
              tableCorr->setVerticalHeaderItem(i2, new QTableWidgetItem(names[i2]));

            // Susceptibility
            if (comboQuantity->currentIndex() == 0) {
              tableCorr->setItem(i1, i2, new QTableWidgetItem(QString::number(model->Susc((ConservedCharge::Name)i, (ConservedCharge::Name)j))));
            }
            else if (comboQuantity->currentIndex() == 1) {
              tableCorr->setItem(i1, i2, new QTableWidgetItem(QString::number(model->Susc((ConservedCharge::Name)i, (ConservedCharge::Name)j) * model->Volume() * pow(model->Parameters().T, 3) * pow(xMath::GeVtoifm(), 3))));
            }
            else {
              double N1 = model->ConservedChargeDensity((ConservedCharge::Name)i) * model->Volume();
              double N2 = model->ConservedChargeDensity((ConservedCharge::Name)j) * model->Volume();
              double wn1 = model->Susc((ConservedCharge::Name)i, (ConservedCharge::Name)i) * model->Volume() * pow(model->Parameters().T, 3) * pow(xMath::GeVtoifm(), 3) / N1;
              double wn2 = model->Susc((ConservedCharge::Name)j, (ConservedCharge::Name)j) * model->Volume() * pow(model->Parameters().T, 3) * pow(xMath::GeVtoifm(), 3) / N2;
              double corr = model->Susc((ConservedCharge::Name)i, (ConservedCharge::Name)j) * model->Volume() * pow(model->Parameters().T, 3) * pow(xMath::GeVtoifm(), 3);
              // Scaled moment
              if (comboQuantity->currentIndex() == 2) {
                if (i != j) {
                  tableCorr->setItem(i, j, new QTableWidgetItem("N/A"));
                }
                else {
                  tableCorr->setItem(i, j, new QTableWidgetItem(QString::number(corr / N1)));
                }
              }
              // Pearson
              else if (comboQuantity->currentIndex() == 3) {
                double var1 = wn1 * N1, var2 = wn2 * N2;
                tableCorr->setItem(i, j, new QTableWidgetItem(QString::number(corr / sqrt(var1 * var2))));
              }
              // Delta
              else if (comboQuantity->currentIndex() == 4) {
                double DeltaN1N2 = -(N1 * wn2 - N2 * wn1) / (N2 - N1);
                tableCorr->setItem(i1, i2, new QTableWidgetItem(QString::number(DeltaN1N2)));
              }
              // Sigma
              else {
                double SigmaN1N2 = (N1 * wn2 + N2 * wn1 - 2. * corr) / (N2 + N1);
                tableCorr->setItem(i1, i2, new QTableWidgetItem(QString::number(SigmaN1N2)));
              }
            }
            i2++;
          }
        }
        i1++;
      }
    }
  }

  // (net)particle-charge correlations
  if (comboType->currentIndex() >= 3 && comboType->currentIndex() <= 4) {
    std::vector<int> flags(4, 1);
    flags[0] = model->TPS()->hasBaryons();
    flags[1] = model->TPS()->hasCharged();
    flags[2] = model->TPS()->hasStrange();
    flags[3] = model->TPS()->hasCharmed();
    {
      int totnum = 0;
      for (int i = 0; i < flags.size(); ++i)
        totnum += flags[i];
      tableCorr->setColumnCount(totnum);
    }
    std::vector<QString> names = { "B", "Q", "S", "C" };

    QVector<long long> idtopdg(0);
    for (int i = 0; i < model->TPS()->ComponentsNumber(); ++i) {
      const ThermalParticle& part = model->TPS()->Particle(i);
      if (comboFeeddown->currentIndex() == 0 || (comboFeeddown->currentIndex() == 1 && part.IsStable())) {
        if (comboType->currentIndex() == 3 || (comboType->currentIndex() == 4 && !part.IsNeutral() && part.PdgId() > 0)) {
          idtopdg.push_back(part.PdgId());
        }
      }
    }

    tableCorr->setRowCount(idtopdg.size());

    for (int i = 0; i < idtopdg.size(); ++i) {
      int i2 = 0;
      QString name = model->TPS()->ParticleByPDG(idtopdg[i]).Name().c_str();
      long long pdg1 = idtopdg[i];
      if (comboType->currentIndex() == 4)
        name = "net-" + name;
      tableCorr->setVerticalHeaderItem(i, new QTableWidgetItem(name));
      for (int cc = 0; cc < flags.size(); ++cc) {
        if (flags[cc]) {
          if (i == 0) {
            tableCorr->setHorizontalHeaderItem(i2, new QTableWidgetItem(names[cc]));
          }

          double corr = 1.;
          if (comboType->currentIndex() == 3) {
            if (comboFeeddown->currentIndex() == 0)
              corr = model->PrimordialParticleChargeSusceptibilityByPdg(pdg1, (ConservedCharge::Name)cc);
            else
              corr = model->FinalParticleChargeSusceptibilityByPdg(pdg1, (ConservedCharge::Name)cc);
          }
          else {
            if (comboFeeddown->currentIndex() == 0)
              corr = model->PrimordialNetParticleChargeSusceptibilityByPdg(pdg1, (ConservedCharge::Name)cc);
            else
              corr = model->FinalNetParticleChargeSusceptibilityByPdg(pdg1, (ConservedCharge::Name)cc);
          }

          if (comboQuantity->currentIndex() == 0)
            tableCorr->setItem(i, i2, new QTableWidgetItem(QString::number(corr)));
          else if (comboQuantity->currentIndex() == 1)
            tableCorr->setItem(i, i2, new QTableWidgetItem(QString::number(corr * model->Volume() * pow(model->Parameters().T, 3) * pow(xMath::GeVtoifm(), 3))));
          else {
            double N1 = 0., N2 = 0.;
            double wn1 = 0., wn2 = 0.;

            if (comboType->currentIndex() == 3) {
              if (comboFeeddown->currentIndex() == 0) {
                N1 = model->GetDensity(pdg1, Feeddown::Primordial) * model->Volume();
                N2 = model->ConservedChargeDensity((ConservedCharge::Name)cc) * model->Volume();
                wn1 = model->ScaledVariancePrimordial(model->TPS()->PdgToId(pdg1));
                wn2 = model->Susc((ConservedCharge::Name)cc, (ConservedCharge::Name)cc) * model->Volume() * pow(model->Parameters().T, 3) * pow(xMath::GeVtoifm(), 3) / N2;
              }
              else {
                N1 = model->GetDensity(pdg1, Feeddown::StabilityFlag) * model->Volume();
                N2 = model->ConservedChargeDensity((ConservedCharge::Name)cc) * model->Volume();
                wn1 = model->ScaledVarianceTotal(model->TPS()->PdgToId(pdg1));
                wn2 = model->Susc((ConservedCharge::Name)cc, (ConservedCharge::Name)cc) * model->Volume() * pow(model->Parameters().T, 3) * pow(xMath::GeVtoifm(), 3) / N2;
              }
            }
            else {
              if (comboFeeddown->currentIndex() == 0) {
                N1 = model->GetDensity(pdg1, Feeddown::Primordial) * model->Volume()
                  - model->GetDensity(-pdg1, Feeddown::Primordial) * model->Volume();
                N2 = model->ConservedChargeDensity((ConservedCharge::Name)cc) * model->Volume();
                wn1 = (model->TwoParticleSusceptibilityPrimordialByPdg(pdg1, pdg1)
                  - model->TwoParticleSusceptibilityPrimordialByPdg(pdg1, -pdg1)
                  - model->TwoParticleSusceptibilityPrimordialByPdg(-pdg1, pdg1)
                  + model->TwoParticleSusceptibilityPrimordialByPdg(-pdg1, -pdg1)) / N1;
                wn1 *= model->Volume() * pow(model->Parameters().T, 3) * pow(xMath::GeVtoifm(), 3);
                wn2 = model->Susc((ConservedCharge::Name)cc, (ConservedCharge::Name)cc) * model->Volume() * pow(model->Parameters().T, 3) * pow(xMath::GeVtoifm(), 3) / N2;
              }
              else {
                N1 = model->GetDensity(pdg1, Feeddown::StabilityFlag) * model->Volume()
                  - model->GetDensity(-pdg1, Feeddown::StabilityFlag) * model->Volume();
                N2 = model->ConservedChargeDensity((ConservedCharge::Name)cc) * model->Volume();
                wn1 = (model->TwoParticleSusceptibilityFinalByPdg(pdg1, pdg1)
                  - model->TwoParticleSusceptibilityFinalByPdg(pdg1, -pdg1)
                  - model->TwoParticleSusceptibilityFinalByPdg(-pdg1, pdg1)
                  + model->TwoParticleSusceptibilityFinalByPdg(-pdg1, -pdg1)) / N1;
                wn1 *= model->Volume() * pow(model->Parameters().T, 3) * pow(xMath::GeVtoifm(), 3);
                wn2 = model->Susc((ConservedCharge::Name)cc, (ConservedCharge::Name)cc) * model->Volume() * pow(model->Parameters().T, 3) * pow(xMath::GeVtoifm(), 3) / N2;
              }
            }

            corr *= model->Volume() * pow(model->Parameters().T, 3) * pow(xMath::GeVtoifm(), 3);

            // Scaled moment
            if (comboQuantity->currentIndex() == 2) {
              tableCorr->setItem(i, i2, new QTableWidgetItem("N/A"));
            }
            // Pearson
            else if (comboQuantity->currentIndex() == 3) {
              double var1 = wn1 * N1, var2 = wn2 * N2;
              tableCorr->setItem(i, i2, new QTableWidgetItem(QString::number(corr / sqrt(var1 * var2))));
            }
            // Delta
            else if (comboQuantity->currentIndex() == 4) {
              double DeltaN1N2 = -(N1 * wn2 - N2 * wn1) / (N2 - N1);
              tableCorr->setItem(i, i2, new QTableWidgetItem(QString::number(DeltaN1N2)));
            }
            // Sigma
            else {
              double SigmaN1N2 = (N1 * wn2 + N2 * wn1 - 2. * corr) / (N2 + N1);
              tableCorr->setItem(i, i2, new QTableWidgetItem(QString::number(SigmaN1N2)));
            }
          }

          i2++;
        }
      }
    }
  }
}

bool CorrelationsDialog::eventFilter(QObject* obj, QEvent* event)
{
  if (obj == tableCorr) {
    if (event->type() == QEvent::KeyPress && static_cast<QKeyEvent*>(event)->matches(QKeySequence::Copy)) {
      copyTableViewSelectionToClipBoard(tableCorr);
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

