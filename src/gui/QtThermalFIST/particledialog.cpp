/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "particledialog.h"

#include <QLayout>
#include <QHeaderView>
#include <QTimer>
#include <QLabel>
#include <QApplication>
#include <QClipboard>
#include <QMenu>
#include <QKeyEvent>
#include <QDebug>

#include "SpectralFunctionDialog.h"
#include "HelperRoutines.h"

#include "HRGBase/ThermalModelBase.h"

using namespace thermalfist;



ParticleDialog::ParticleDialog(QWidget* parent, ThermalModelBase* mod, int ParticleID) :
  QDialog(parent), model(mod), pid(ParticleID)
{
  QVBoxLayout* layout = new QVBoxLayout();

  QFont font = QApplication::font();
  font.setPointSize(font.pointSize() + 2);
  QLabel* labInfo = new QLabel(tr("Information:"));
  labInfo->setFont(font);
  data = new QTextEdit();
  QFont font2("Monospace");
  font2.setStyleHint(QFont::TypeWriter);
  data->setFont(font2);
  data->setReadOnly(true);
  data->setPlainText(GetParticleInfo());
  data->setMinimumHeight(350);
  data->setMinimumWidth(300);

  QLabel* labDecays = new QLabel(QString(model->TPS()->Particles()[pid].Name().c_str()) + tr(" decays:"));
  labDecays->setFont(font);
  myModel = new DecayTableModel(0, &model->TPS()->Particle(pid), model->TPS());
  tableDecays = new QTableView();
  tableDecays->setModel(myModel);
  tableDecays->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);

  QHBoxLayout* layoutButtons = new QHBoxLayout();
  layoutButtons->setAlignment(Qt::AlignLeft);

  layout->addWidget(labInfo);
  layout->addWidget(data);
  if (model->TPS()->Particles()[pid].Decays().size() > 0) {
    layout->addWidget(labDecays);
    layout->addWidget(tableDecays);
  }

  if (model->TPS()->Particles()[pid].ResonanceWidth() > 0.) {
    buttonSpectralFunction = new QPushButton(tr("Spectral function..."));
    connect(buttonSpectralFunction, SIGNAL(clicked()), this, SLOT(showSpectralFunction()));
    layout->addWidget(buttonSpectralFunction, 0, Qt::AlignRight);
  }

  QLabel* labProd = new QLabel(tr("Production of ") + QString::fromStdString(model->TPS()->Particles()[pid].Name()));
  labProd->setFont(font);

  QLabel* labelPrimDensity = new QLabel(tr("Primordial density = "));
  if (model->IsCalculated()) labelPrimDensity->setText(labelPrimDensity->text() +
    QString::number(model->Densities()[pid]) +
    " fm<sup>-3</sup>");

  QLabel* labelPrimMultiplicity = new QLabel(tr("Primordial yield = "));
  if (model->IsCalculated()) labelPrimMultiplicity->setText(labelPrimMultiplicity->text() +
    QString::number(model->Densities()[pid] * model->Volume()));

  QLabel* labelTotalMultiplicity = new QLabel(tr("Total yield = "));
  if (model->IsCalculated()) labelTotalMultiplicity->setText(labelTotalMultiplicity->text() +
    QString::number(model->TotalDensities()[pid] * model->Volume()));

  QLabel* labelStrongMultiplicity = new QLabel(tr("Primordial + strong decays = "));
  if (model->IsCalculated()) labelStrongMultiplicity->setText(labelStrongMultiplicity->text() +
    QString::number(model->GetDensity(model->TPS()->Particle(pid).PdgId(), Feeddown::Strong) * model->Volume()));

  QLabel* labelEMMultiplicity = new QLabel(tr("Primordial + strong + EM decays = "));
  if (model->IsCalculated()) labelEMMultiplicity->setText(labelEMMultiplicity->text() +
    QString::number(model->GetDensity(model->TPS()->Particle(pid).PdgId(), Feeddown::Electromagnetic) * model->Volume()));

  QLabel* labelWeakMultiplicity = new QLabel(tr("Primordial + strong + EM + weak decays = "));
  if (model->IsCalculated()) labelWeakMultiplicity->setText(labelWeakMultiplicity->text() +
    QString::number(model->GetDensity(model->TPS()->Particle(pid).PdgId(), Feeddown::Weak) * model->Volume()));

  tableSources = new QTableWidget();
  tableSources->setColumnCount(3);
  tableSources->setRowCount(model->TPS()->DecayContributionsByFeeddown()[Feeddown::StabilityFlag][pid].size() + 1);
  tableSources->verticalHeader()->hide();
  tableSources->setHorizontalHeaderItem(0, new QTableWidgetItem(tr("Source")));
  tableSources->setHorizontalHeaderItem(1, new QTableWidgetItem(tr("Multiplicity")));
  tableSources->setHorizontalHeaderItem(2, new QTableWidgetItem(tr("Fraction (%)")));

  if (model->IsCalculated()) {

    std::vector< std::pair<double, int> > sources = model->TPS()->DecayContributionsByFeeddown()[Feeddown::StabilityFlag][pid];
    for (int i = 0; i < sources.size(); ++i)
      sources[i].first *= model->Densities()[sources[i].second];

    qSort(sources.begin(), sources.end());

    tableSources->setItem(0, 0, new QTableWidgetItem(tr("Primordial")));
    tableSources->setItem(0, 1, new QTableWidgetItem(QString::number(model->Densities()[pid] * model->Volume())));
    tableSources->setItem(0, 2, new QTableWidgetItem(QString::number(model->Densities()[pid] / model->TotalDensities()[pid] * 100.)));


    for (int j = 0; j < 3; ++j) {
      QTableWidgetItem* item = tableSources->item(0, j);
      item->setFlags(item->flags() & ~Qt::ItemIsEditable);
    }

    for (int i = 0; i < sources.size(); ++i)
    {
      int tindex = sources.size() - i - 1;
      int decaypartid = sources[tindex].second;
      tableSources->setItem(i + 1, 0, new QTableWidgetItem(tr("Decays from primordial ") + QString(model->TPS()->Particles()[decaypartid].Name().c_str())));
      tableSources->setItem(i + 1, 1, new QTableWidgetItem(QString::number(sources[tindex].first * model->Volume())));
      tableSources->setItem(i + 1, 2, new QTableWidgetItem(QString::number(sources[tindex].first / model->TotalDensities()[pid] * 100.)));

      for (int j = 0; j < 3; ++j) {
        QTableWidgetItem* item = tableSources->item(i + 1, j);
        item->setFlags(item->flags() & ~Qt::ItemIsEditable);
      }
    }
  }

  tableSources->resizeColumnsToContents();
  tableSources->setContextMenuPolicy(Qt::CustomContextMenu);
  connect(tableSources, SIGNAL(customContextMenuRequested(QPoint)), this, SLOT(contextMenuRequest(QPoint)));
  tableSources->installEventFilter(this);

  QVBoxLayout* layoutProd = new QVBoxLayout();
  layoutProd->setAlignment(Qt::AlignTop);
  layoutProd->addWidget(labProd);
  layoutProd->addWidget(labelPrimDensity);
  layoutProd->addWidget(labelPrimMultiplicity);
  layoutProd->addWidget(labelTotalMultiplicity);
  layoutProd->addWidget(labelStrongMultiplicity);
  layoutProd->addWidget(labelEMMultiplicity);
  layoutProd->addWidget(labelWeakMultiplicity);
  layoutProd->addWidget(tableSources);

  QHBoxLayout* mainLayout = new QHBoxLayout();

  mainLayout->addLayout(layout);
  mainLayout->addLayout(layoutProd);

  setLayout(mainLayout);

  this->setWindowTitle(QString(model->TPS()->Particles()[pid].Name().c_str()));

  QTimer::singleShot(0, this, SLOT(checkFixTableSize()));
}

QString ParticleDialog::GetParticleInfo() {
  QString ret = "";

  ret += tr("Name").leftJustified(20) + " = ";
  ret += QString(model->TPS()->Particles()[pid].Name().c_str()) + "\r\n";

  ret += tr("PDG ID").leftJustified(20) + " = ";
  ret += QString::number(model->TPS()->Particles()[pid].PdgId()) + "\r\n";

  ret += tr("Mass").leftJustified(20) + " = ";
  ret += QString::number(model->TPS()->Particles()[pid].Mass() * 1.e3) + " " + tr("MeV") + "\r\n";

  ret += tr("Type").leftJustified(20) + " = ";
  if (model->TPS()->Particles()[pid].BaryonCharge() == 0) ret += tr("Meson") + "\r\n";
  else if (model->TPS()->Particles()[pid].BaryonCharge() == 1) ret += tr("Baryon") + "\r\n";
  else if (model->TPS()->Particles()[pid].BaryonCharge() == -1) ret += "Anti-Baryon\r\n";
  else if (model->TPS()->Particles()[pid].BaryonCharge() > 1) ret += tr("Light nucleus") + "\r\n";
  else if (model->TPS()->Particles()[pid].BaryonCharge() < -1) ret += tr("Light antinucleus") + "\r\n";

  ret += tr("Stability flag").leftJustified(20) + " = ";
  if (model->TPS()->Particles()[pid].IsStable()) ret += "1\r\n";
  else ret += "0\r\n";
  //if (model->TPS()->Particles()[pid].IsStable()) ret += tr("Yes") + "\r\n";
  //else ret += "No\r\n";

  ret += tr("Decay type").leftJustified(20) + " = ";
  ParticleDecayType::DecayType dtype = model->TPS()->Particles()[pid].DecayType();
  if (dtype == ParticleDecayType::Strong) ret += "Strong\r\n";
  else if (dtype == ParticleDecayType::Electromagnetic) ret += "Electromagnetic\r\n";
  else if (dtype == ParticleDecayType::Weak) ret += "Weak\r\n";
  else ret += "Stable\r\n";

  ret += tr("Neutral?").leftJustified(20) + " = ";
  if (model->TPS()->Particles()[pid].IsNeutral()) ret += tr("Yes") + "\r\n";
  else ret += "No\r\n";

  ret += tr("Statistics").leftJustified(20) + " = ";
  double tstat = model->TPS()->Particles()[pid].Statistics();
  if (tstat > 0)
    ret += "+";
  ret += QString::number(tstat);
  if (tstat == 1.)
    ret += " (Fermi-Dirac)";
  else if (tstat == -1.)
    ret += " (Bose-Einstein)";
  else if (tstat == 0.)
    ret += " (Boltzmann)";
  ret += "\r\n";

  ret += tr("Degeneracy").leftJustified(20) + " = ";
  ret += QString::number(model->TPS()->Particles()[pid].Degeneracy()) + "\r\n";

  ret += tr("Electric charge").leftJustified(20) + " = ";
  ret += QString::number(model->TPS()->Particles()[pid].ElectricCharge()) + "\r\n";

  ret += tr("Strangeness").leftJustified(20) + " = ";
  ret += QString::number(model->TPS()->Particles()[pid].Strangeness()) + "\r\n";

  ret += tr("Charm").leftJustified(20) + " = ";
  ret += QString::number(model->TPS()->Particles()[pid].Charm()) + "\r\n";

  ret += tr("|u,d|").leftJustified(20) + " = ";
  ret += QString::number(model->TPS()->Particles()[pid].AbsoluteQuark()) + "\r\n";

  ret += tr("|s|").leftJustified(20) + " = ";
  ret += QString::number(model->TPS()->Particles()[pid].AbsoluteStrangeness()) + "\r\n";

  ret += tr("|c|").leftJustified(20) + " = ";
  ret += QString::number(model->TPS()->Particles()[pid].AbsoluteCharm()) + "\r\n";

  ret += tr("Isospin").leftJustified(20) + " = ";
  ret += QString::number(2 * model->TPS()->Particles()[pid].ElectricCharge() - model->TPS()->Particles()[pid].BaryonCharge() - model->TPS()->Particles()[pid].Strangeness() - model->TPS()->Particles()[pid].Charm()) + "\r\n";

  ret += tr("I3").leftJustified(20) + " = ";
  ret += QString::number(model->TPS()->Particles()[pid].ElectricCharge() - (model->TPS()->Particles()[pid].BaryonCharge() + model->TPS()->Particles()[pid].Strangeness() + model->TPS()->Particles()[pid].Charm()) / 2.) + "\r\n";

  if (model->TPS()->Particles()[pid].ResonanceWidth() > 1e-5) {
    ret += tr("Decay width").leftJustified(20) + " = ";
    ret += QString::number(model->TPS()->Particles()[pid].ResonanceWidth() * 1.e3) + " " + tr("MeV") + "\r\n";

    if (model->TPS()->Particles()[pid].DecayThresholdMass() > 1e-5) {
      ret += tr("Threshold mass").leftJustified(20) + " = ";
      ret += QString::number(model->TPS()->Particles()[pid].DecayThresholdMass() * 1.e3) + " " + tr("MeV");
      ret += "\r\n";
    }
  }

  return ret;
}

void ParticleDialog::checkFixTableSize() {
  int tle = 0;
  for (int i = 0; i < tableDecays->horizontalHeader()->count(); ++i)
    tle += tableDecays->horizontalHeader()->sectionSize(i);
  tableDecays->setMinimumWidth(tle + 23);
  tle = 0;
  for (int i = 0; i < tableSources->horizontalHeader()->count(); ++i)
    tle += tableSources->horizontalHeader()->sectionSize(i);
  tableSources->setMinimumWidth(tle + 23);
}

void ParticleDialog::addDecay() {
  myModel->addDecay();
}

void ParticleDialog::removeDecay() {
  QModelIndexList selectedList = tableDecays->selectionModel()->selectedRows();
  //qDebug() << selectedList.count() << "\n";
  for (unsigned int i = 0; i < selectedList.count(); ++i)
    myModel->removeDecay(selectedList.at(i).row());
}

void ParticleDialog::addColumn() {
  myModel->addColumn();
  checkFixTableSize();
}

void ParticleDialog::removeColumn() {
  myModel->removeColumn();
}

void ParticleDialog::showSpectralFunction()
{
  SpectralFunctionDialog dialog(this, &model->TPS()->Particle(pid), model->Parameters().T, model->ChemicalPotential(pid), static_cast<int>(model->TPS()->ResonanceWidthIntegrationType() == ThermalParticle::eBW));
  dialog.setWindowFlags(Qt::Window);
  dialog.setMinimumSize(QSize(800, 400));
  dialog.exec();
}


void ParticleDialog::contextMenuRequest(QPoint pos)
{
  QMenu* menu = new QMenu(this);
  menu->setAttribute(Qt::WA_DeleteOnClose);

  menu->addAction("Copy table to clipboard", this, SLOT(copyFeeddownTable()));

  menu->popup(tableSources->mapToGlobal(pos));
}

void ParticleDialog::copyFeeddownTable()
{
  QString table_text = "";

  for (int i = 0; i < tableSources->rowCount(); ++i) {
    for (int j = 0; j < tableSources->columnCount(); ++j) {
      if (tableSources->item(i, j) != NULL)
        table_text.append(tableSources->item(i, j)->text());
      if (j != tableSources->columnCount() - 1)
        table_text.append('\t');
      else
        table_text.append('\n');
    }
  }

  QApplication::clipboard()->setText(table_text);
}

bool ParticleDialog::eventFilter(QObject* obj, QEvent* event)
{
  if (obj == tableSources) {
    if (event->type() == QEvent::KeyPress && static_cast<QKeyEvent*>(event)->matches(QKeySequence::Copy)) {
      copyTableViewSelectionToClipBoard(tableSources);
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

