/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "listeditortab.h"

#include <algorithm>
#include <utility>
#include <sstream>

#include <QLayout>
#include <QFileDialog>
#include <QLabel>
#include <QHeaderView>
#include <QGroupBox>
#include <QApplication>
#include <QMessageBox>
#include <QInputDialog>
#include <QDebug>
#include "listtablemodel.h"
#include "decayseditor.h"

#include "ThermalFISTConfig.h"
#include "HRGBase/ThermalParticleSystem.h"

using namespace thermalfist;


ListEditorTab::ListEditorTab(QWidget *parent, thermalfist::ThermalModelBase *model_in)
  : QWidget(parent)
{
  QHBoxLayout *DataEditLay = new QHBoxLayout();

  QVBoxLayout *dataLayv = new QVBoxLayout();

  myModel = new ListTableModel(0);
  tableParticles = new QTableView();
  tableParticles->setModel(myModel);
  tableParticles->setSelectionBehavior(QAbstractItemView::SelectRows);
  tableParticles->setSelectionMode(QAbstractItemView::SingleSelection);
  tableParticles->resizeColumnsToContents();
  connect(tableParticles->selectionModel(), SIGNAL(selectionChanged(QItemSelection, QItemSelection)), this, SLOT(changedRow()));
  connect(tableParticles, SIGNAL(doubleClicked(QModelIndex)), this, SLOT(editDecaysDoubleClick(QModelIndex)));
  //tableParticles->setSortingEnabled(true);
  //tableParticles->show();
  //tableParticles->setColumnCount(7);
  //tableParticles->setRowCount(3);

  QHBoxLayout *dataLay2 = new QHBoxLayout();
  dataLay2->setAlignment(Qt::AlignLeft);
  buttonNewPart = new QPushButton(tr("Add particle"));
  buttonRemovePart = new QPushButton(tr("Remove particle"));

  connect(buttonNewPart, SIGNAL(clicked()), this, SLOT(addParticle()));
  connect(buttonRemovePart, SIGNAL(clicked()), this, SLOT(removeParticle()));

  dataLay2->addWidget(buttonNewPart);
  dataLay2->addWidget(buttonRemovePart);

  dataLayv->addWidget(tableParticles);
  dataLayv->addLayout(dataLay2);

  // QGroupBox *grEditor = new QGroupBox(tr("Editor"));

  QVBoxLayout *editorLay = new QVBoxLayout();
  editorLay->setContentsMargins(15, 0, 0, 0);
  //editorLay->setMargin(10);
  editorLay->setAlignment(Qt::AlignTop);

  QGroupBox *grParticle = new QGroupBox(tr("Particle properties:"));

  QVBoxLayout *particleLay = new QVBoxLayout();

  QHBoxLayout *edBaseLay = new QHBoxLayout();
  edBaseLay->setAlignment(Qt::AlignLeft);
  QLabel *labelPDGID = new QLabel(tr("PDG ID:"));

  lePDGID = new QLineEdit();
  QRegExp rx("-?\\d{1,15}");
  QValidator *validator = new QRegExpValidator(rx, this);
  lePDGID->setValidator(validator);
  lePDGID->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
  connect(lePDGID, SIGNAL(editingFinished()), this, SLOT(PDGEdited()));

  QLabel *labelName = new QLabel(tr("Name:"));
  leName = new QLineEdit();
  leName->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
  leName->setMaximumWidth(150);
  connect(leName, SIGNAL(editingFinished()), this, SLOT(NameEdited()));
  
  edBaseLay->addWidget(labelPDGID);
  //edBaseLay->addWidget(spinPDGID);
  edBaseLay->addWidget(lePDGID, 0);
  edBaseLay->addSpacing(15);
  edBaseLay->addWidget(labelName);
  edBaseLay->addWidget(leName);
  

  QHBoxLayout *degeneracyLay = new QHBoxLayout();
  degeneracyLay->setAlignment(Qt::AlignLeft);

  QLabel *labelMass = new QLabel(tr("Mass [GeV]:"));
  spinMass = new QDoubleSpinBox();
  spinMass->setButtonSymbols(QAbstractSpinBox::NoButtons);
  spinMass->setDecimals(6);
  spinMass->setMinimum(0.);
  spinMass->setMaximum(127.);
  connect(spinMass, SIGNAL(editingFinished()), this, SLOT(MassEdited()));

  QLabel *labelDegeneracy = new QLabel(tr("Degeneracy:"));
  spinDegeneracy = new QDoubleSpinBox();
  spinDegeneracy->setMinimum(0.);
  connect(spinDegeneracy, SIGNAL(editingFinished()), this, SLOT(SpinEdited()));

  //QLabel *labelRadius = new QLabel(tr("Hard core radius: "));
  //spinRadius = new QDoubleSpinBox();
  //spinRadius->setMinimum(0.);
  //connect(spinRadius, SIGNAL(editingFinished()), this, SLOT(RadiusEdited()));
  degeneracyLay->addWidget(labelMass);
  degeneracyLay->addWidget(spinMass);
  degeneracyLay->addSpacing(15);
  degeneracyLay->addWidget(labelDegeneracy);
  degeneracyLay->addWidget(spinDegeneracy);
  //degeneracyLay->addSpacing(15);
  //degeneracyLay->addWidget(labelRadius);
  //degeneracyLay->addWidget(spinRadius);


  QGroupBox *grQuant = new QGroupBox(tr("Quantum numbers:"));

  QVBoxLayout *layAllQuant = new QVBoxLayout();
  QHBoxLayout *layQuant = new QHBoxLayout();
  //QGridLayout *layQuant = new QGridLayout();
  layQuant->setAlignment(Qt::AlignLeft);
  //layQuant->setSpacing(15);
  //layQuant->setHorizontalSpacing(15);
  QLabel *labelBaryon = new QLabel(tr("B:"));
  spinBaryon = new QSpinBox();
  spinBaryon->setMinimum(0);
  spinBaryon->setMaximum(10);
  QLabel *labelCharge = new QLabel(tr("Q:"));
  spinCharge = new QSpinBox();
  spinCharge->setMinimum(-10);
  spinCharge->setMaximum(10);

  QLabel *labelStrangeness = new QLabel(tr("S:"));
  spinStrangeness = new QSpinBox();
  spinStrangeness->setMinimum(-10);
  spinStrangeness->setMaximum(10);
  QLabel *labelCharm = new QLabel(tr("C:"));
  spinCharm = new QSpinBox();
  spinCharm->setMinimum(-10);
  spinCharm->setMaximum(10);

  connect(spinBaryon, SIGNAL(editingFinished()), this, SLOT(ChargeEdited()));
  connect(spinCharge, SIGNAL(editingFinished()), this, SLOT(ChargeEdited()));
  connect(spinStrangeness, SIGNAL(editingFinished()), this, SLOT(ChargeEdited()));
  connect(spinCharm, SIGNAL(editingFinished()), this, SLOT(ChargeEdited()));


  

  layQuant->addWidget(labelBaryon);
  layQuant->addWidget(spinBaryon);
  layQuant->addSpacing(15);
  //layQuant->addSpacing(15);
  layQuant->addWidget(labelCharge);
  layQuant->addWidget(spinCharge);
  layQuant->addSpacing(15);
  //layQuant->addSpacing(15);
  layQuant->addWidget(labelStrangeness);
  layQuant->addWidget(spinStrangeness);
  layQuant->addSpacing(15);
  layQuant->addWidget(labelCharm);
  layQuant->addWidget(spinCharm);

  
  QHBoxLayout *layAbsQuant = new QHBoxLayout();
  layAbsQuant->setAlignment(Qt::AlignLeft);
  //layAbsQuant->setSpacing(15);
  QLabel *labelAbsoluteLightQuark = new QLabel(tr("|q<sub>l</sub>|:"));
  labelAbsoluteLightQuarkValue = new QLabel(tr("0"));
  //spinAbsoluteStrangeness = new QDoubleSpinBox();
  //spinAbsoluteStrangeness->setDecimals(3);
  //spinAbsoluteStrangeness->setMinimum(0.);
  //spinAbsoluteStrangeness->setMaximum(10);
  //connect(spinAbsoluteStrangeness, SIGNAL(editingFinished()), this, SLOT(ChargeEdited()));

  QLabel *labelAbsoluteStrangeness = new QLabel(tr("|s|:"));
  spinAbsoluteStrangeness = new QDoubleSpinBox();
  spinAbsoluteStrangeness->setDecimals(5);
  spinAbsoluteStrangeness->setMinimum(0.);
  //spinAbsoluteStrangeness->setMaximum(10);
  connect(spinAbsoluteStrangeness, SIGNAL(editingFinished()), this, SLOT(ChargeEdited()));

  QLabel *labelAbsoluteCharm = new QLabel(tr("|c|:"));
  spinAbsoluteCharm = new QDoubleSpinBox();
  spinAbsoluteCharm->setDecimals(5);
  spinAbsoluteCharm->setMinimum(0.);
  //spinAbsoluteStrangeness->setMaximum(10);
  connect(spinAbsoluteCharm, SIGNAL(editingFinished()), this, SLOT(ChargeEdited()));

  layAbsQuant->addWidget(labelAbsoluteLightQuark);
  layAbsQuant->addWidget(labelAbsoluteLightQuarkValue);
  layAbsQuant->addSpacing(15);
  layAbsQuant->addWidget(labelAbsoluteStrangeness);
  layAbsQuant->addWidget(spinAbsoluteStrangeness);
  layAbsQuant->addSpacing(15);
  layAbsQuant->addWidget(labelAbsoluteCharm);
  layAbsQuant->addWidget(spinAbsoluteCharm);

  layAllQuant->addLayout(layQuant);
  layAllQuant->addLayout(layAbsQuant);
  grQuant->setLayout(layAllQuant);

  checkBoxStable = new QCheckBox(tr("Stable"));
  connect(checkBoxStable, SIGNAL(toggled(bool)), this, SLOT(StableEdited()));

  QGroupBox *grDecays = new QGroupBox(tr("Decays:"));

  QVBoxLayout *layDecays = new QVBoxLayout();
  QHBoxLayout *layDecays1 = new QHBoxLayout();
  layDecays1->setAlignment(Qt::AlignLeft);
  
  QLabel *labelWidth = new QLabel(tr("Width [MeV]:"));
  spinWidth = new QDoubleSpinBox();
  spinWidth->setButtonSymbols(QAbstractSpinBox::NoButtons);
  spinWidth->setMinimum(0.);
  spinWidth->setMaximum(100000.);
  spinWidth->setDecimals(6);
  QLabel *labelThreshold = new QLabel(tr("Threshold [GeV]:"));
  spinThreshold = new QDoubleSpinBox();
  spinThreshold->setButtonSymbols(QAbstractSpinBox::NoButtons);
  spinThreshold->setMinimum(0.);
  spinThreshold->setDecimals(6);

  connect(spinWidth, SIGNAL(editingFinished()), this, SLOT(DecaysEdited()));
  connect(spinThreshold, SIGNAL(editingFinished()), this, SLOT(DecaysEdited()));

  QHBoxLayout *layDecays2 = new QHBoxLayout();
  layDecays2->setAlignment(Qt::AlignLeft);

  buttonDecayModes = new QPushButton(tr("Decay modes..."));
  connect(buttonDecayModes, SIGNAL(clicked()), this, SLOT(editDecays()));

  buttonSetThreshold = new QPushButton(tr("Set threshold from decay channels..."));
  connect(buttonSetThreshold, SIGNAL(clicked()), this, SLOT(setThreshold()));

  
  layDecays1->addWidget(labelWidth);
  layDecays1->addWidget(spinWidth);
  layDecays1->addSpacing(15);
  layDecays1->addWidget(labelThreshold);
  layDecays1->addWidget(spinThreshold);

  layDecays2->addWidget(buttonDecayModes);
  layDecays2->addSpacing(15);
  layDecays2->addWidget(buttonSetThreshold);

  layDecays->addLayout(layDecays1);
  layDecays->addLayout(layDecays2);

  grDecays->setLayout(layDecays);

  buttonReset = new QPushButton(tr("Discard changes"));
  connect(buttonReset, SIGNAL(clicked()), this, SLOT(resetParticle()));

  particleLay->addLayout(edBaseLay);
  particleLay->addLayout(degeneracyLay);
  particleLay->addWidget(grQuant);
  particleLay->addWidget(grDecays);
  particleLay->addWidget(checkBoxStable);
  particleLay->addWidget(buttonReset, 0, Qt::AlignLeft);

  grParticle->setLayout(particleLay);

  QGroupBox *grGlobal = new QGroupBox(tr("Global modifications:"));
  QVBoxLayout *layoutGlobal = new QVBoxLayout();
  layoutGlobal->setAlignment(Qt::AlignLeft);

  
  buttonMassesWidthsFromPDG = new QPushButton(tr("Set masses/widths/names from PDG file..."));
  connect(buttonMassesWidthsFromPDG, SIGNAL(clicked()), this, SLOT(loadMassesWidthsFromPdg()));

  buttonSetAllThresholds = new QPushButton(tr("Set thresholds from decay channels..."));
  connect(buttonSetAllThresholds, SIGNAL(clicked()), this, SLOT(setAllThresholds()));

  buttonResetAll = new QPushButton(tr("Discard changes"));
  connect(buttonResetAll, SIGNAL(clicked()), this, SLOT(resetAll()));

  //layoutGlobal->addLayout(dataLay2);
  layoutGlobal->addWidget(buttonMassesWidthsFromPDG, 0, Qt::AlignLeft);
  layoutGlobal->addWidget(buttonSetAllThresholds, 0, Qt::AlignLeft);
  layoutGlobal->addWidget(buttonResetAll, 0, Qt::AlignLeft);

  grGlobal->setLayout(layoutGlobal);

  QHBoxLayout *layApplyOrSave = new QHBoxLayout();
  layApplyOrSave->setAlignment(Qt::AlignLeft);
  buttonApply = new QPushButton(tr("Apply changes"));
  connect(buttonApply, SIGNAL(clicked()), this, SLOT(applyChanges()));
  buttonSaveToFile = new QPushButton(tr("Save list to file..."));
  connect(buttonSaveToFile, SIGNAL(clicked()), this, SLOT(saveToFile()));

  layApplyOrSave->addWidget(buttonApply);
  layApplyOrSave->addWidget(buttonSaveToFile);

  editorLay->addWidget(grParticle);
  editorLay->addWidget(grGlobal);
  editorLay->addLayout(layApplyOrSave);

  DataEditLay->addLayout(dataLayv);
  DataEditLay->addLayout(editorLay);

  //tableParticles->verticalHeader()->hide();

  QVBoxLayout *mainLayout = new QVBoxLayout;
  mainLayout->addLayout(DataEditLay);
  setLayout(mainLayout);

  model = model_in;

  setParticleList(model->TPS());

  currentPath = "";
}

ListEditorTab::~ListEditorTab()
{
  delete myModel;
}

int ListEditorTab::getCurrentRow()
{
  QModelIndexList selectedList = tableParticles->selectionModel()->selectedRows();
  if (selectedList.count() == 0) return -1;
  return selectedList.at(0).row();
}


void ListEditorTab::applyChanges()
{
  if (!myModel->haveChanges()) {
    QMessageBox::information(this, tr("Apply changes"), tr("No new changes detected!"));
    return;
  }
  
  QMessageBox msgBox;
  QPushButton *applyButton = msgBox.addButton(tr("Apply"), QMessageBox::ApplyRole);
  QPushButton *abortButton = msgBox.addButton(tr("No"), QMessageBox::NoRole);
  msgBox.setWindowTitle(tr("Apply changes"));
  msgBox.setText(tr("Do you want to apply changes to the particle list to use in the current session?"));

  msgBox.exec();
  
  if (msgBox.clickedButton() == abortButton || !myModel->haveChanges())
    return;

  std::vector<ThermalParticle> particles = myModel->GetTPS()->Particles();
  myModel->GetTPS()->SetTableFromVector(particles, true);
  
  myModel->applyChanges();

  *model->TPS() = *myModel->GetTPS();
  model->ChangeTPS(model->TPS());

  resetTPS();

}

void ListEditorTab::changedRow()
{
  QModelIndexList selectedList = tableParticles->selectionModel()->selectedRows();
  for (int i = 0; i<selectedList.count() && i<1; i++)
  {
    //QMessageBox::information(this,"", QString::number(selectedList.at(i).row()));
    int row = selectedList.at(i).row();
    ThermalParticleSystem *TPS = myModel->GetTPS();
    int index = myModel->GetRowToParticle()[row];
    //spinPDGID->setValue(TPS->Particles()[index].PdgId());
    lePDGID->setText(QString::number(TPS->Particles()[index].PdgId()));
    leName->setText(TPS->Particles()[index].Name().c_str());
    spinMass->setValue(TPS->Particles()[index].Mass());
    spinDegeneracy->setValue(TPS->Particles()[index].Degeneracy());
    //spinRadius->setValue(TPS->Particles()[index].fRadius);
    spinBaryon->setValue(TPS->Particles()[index].BaryonCharge());
    spinCharge->setValue(TPS->Particles()[index].ElectricCharge());
    spinStrangeness->setValue(TPS->Particles()[index].Strangeness());
    spinAbsoluteStrangeness->setValue(TPS->Particles()[index].AbsoluteStrangeness());
    checkBoxStable->setChecked(TPS->Particles()[index].IsStable());
    spinWidth->setValue(TPS->Particles()[index].ResonanceWidth() * 1.e3);
    spinThreshold->setValue(TPS->Particles()[index].DecayThresholdMass());
    spinCharm->setValue(TPS->Particles()[index].Charm());
    spinAbsoluteCharm->setValue(TPS->Particles()[index].AbsoluteCharm());
    labelAbsoluteLightQuarkValue->setText(QString::number(TPS->Particles()[index].AbsoluteQuark()));
  }
}

void ListEditorTab::editDecays()
{
  int row = getCurrentRow();
  if (row >= 0) {
    decayseditor dialog(this, &myModel->GetTPS()->Particle(myModel->GetRowToParticle()[row]), myModel->GetTPS());
    //dialog.setWordCount(document().wordCount());
    dialog.exec();
  }
}

void ListEditorTab::PDGEdited()
{
  int row = getCurrentRow();
  if (row >= 0) {
    long long pdgedit = lePDGID->text().toLongLong();
    if (myModel->GetTPS()->PdgToId(pdgedit) == -1) {
      int index = myModel->GetRowToParticle()[row];

      ThermalParticle &part = myModel->GetTPS()->Particle(index);

      int oldPDGID = part.PdgId();
      part.SetPdgId(pdgedit);

      if (!part.IsNeutral()) {
        int indexa = myModel->GetTPS()->PdgToId(-oldPDGID);
        if (indexa != -1) {
          ThermalParticle &partanti = myModel->GetTPS()->Particle(indexa);
          partanti.SetPdgId(-pdgedit);
        }
      }

      myModel->GetTPS()->FillPdgMap();

      myModel->updateParticle(row);
    }
  }
}

void ListEditorTab::SpinEdited()
{
  int row = getCurrentRow();
  if (row >= 0) {
    int index = myModel->GetRowToParticle()[row];
    myModel->GetTPS()->Particle(index).SetDegeneracy(spinDegeneracy->value());
    myModel->updateParticle(row);
  }
}

void ListEditorTab::NameEdited()
{
  int row = getCurrentRow();
  if (row >= 0) {
    int index = myModel->GetRowToParticle()[row];
    ThermalParticle &part = myModel->GetTPS()->Particle(index);

    part.SetName(leName->text().toStdString());
    if (!part.IsNeutral()) {
      int indexa = myModel->GetTPS()->PdgToId(-myModel->GetTPS()->Particles()[index].PdgId());
      if (indexa != -1) {
        std::string tname = leName->text().toStdString();
        int bary = part.BaryonCharge();
        if (bary == 0 && tname[tname.size() - 1] == '+') tname[tname.size() - 1] = '-';
        else if (bary == 0 && tname[tname.size() - 1] == '-') tname[tname.size() - 1] = '+';
        else tname += "anti-" + tname;
        myModel->GetTPS()->Particle(indexa).SetName(tname);
      }
    }
    myModel->updateParticle(row);
  }
}

void ListEditorTab::MassEdited()
{
  int row = getCurrentRow();
  if (row >= 0) {
    int index = myModel->GetRowToParticle()[row];
    myModel->GetTPS()->Particle(index).SetMass(spinMass->value());
    myModel->updateParticle(row);
  }
}

void ListEditorTab::DecaysEdited()
{
  int row = getCurrentRow();
  if (row >= 0) {
    int index = myModel->GetRowToParticle()[row];
    ThermalParticle &part = myModel->GetTPS()->Particle(index);
    part.SetResonanceWidth(spinWidth->value() * 1.e-3);
    part.SetDecayThresholdMass(spinThreshold->value());
    myModel->updateParticle(row);
  }
}

void ListEditorTab::ChargeEdited()
{
  int row = getCurrentRow();
  if (row >= 0) {
    int index = myModel->GetRowToParticle()[row];
    ThermalParticle &part = myModel->GetTPS()->Particle(index);
    part.SetBaryonCharge(spinBaryon->value());
    part.SetElectricCharge(spinCharge->value());
    part.SetStrangenessCharge(spinStrangeness->value());
    part.SetCharm(spinCharm->value());
    part.SetAbsoluteStrangeness(spinAbsoluteStrangeness->value());
    part.SetAbsoluteCharm(spinAbsoluteCharm->value());
    if (part.BaryonCharge() & 1) part.SetStatistics(1);
    else part.SetStatistics(-1);

    if (!myModel->GetTPS()->Particles()[index].IsNeutral()) {
      int pdgid = myModel->GetTPS()->Particles()[index].PdgId();
      if (myModel->GetTPS()->PdgToId(-pdgid) != -1) {
        int indexa = myModel->GetTPS()->PdgToId(-pdgid);
        ThermalParticle &partanti = myModel->GetTPS()->Particle(indexa);
        partanti.SetBaryonCharge(-spinBaryon->value());
        partanti.SetElectricCharge(-spinCharge->value());
        partanti.SetStrangenessCharge(-spinStrangeness->value());
        partanti.SetCharm(-spinCharm->value());
        partanti.SetAbsoluteStrangeness(spinAbsoluteStrangeness->value());
        partanti.SetAbsoluteCharm(spinAbsoluteCharm->value());
        if (partanti.BaryonCharge() & 1) partanti.SetStatistics(1);
        else partanti.SetStatistics(-1);
      }
      else {
        std::string tname = myModel->GetTPS()->Particles()[index].Name();
        int bary = myModel->GetTPS()->Particles()[index].BaryonCharge();
        if (bary == 0 && tname[tname.size() - 1] == '+') tname[tname.size() - 1] = '-';
        else if (bary == 0 && tname[tname.size() - 1] == '-') tname[tname.size() - 1] = '+';
        else tname += std::string("bar");

        myModel->GetTPS()->AddParticle(ThermalParticle(
          myModel->GetTPS()->Particles()[index].IsStable(),
          tname,
          -pdgid,
          part.Degeneracy(),
          part.Statistics(),
          part.Mass(),
          -part.Strangeness(),
          -part.BaryonCharge(),
          -part.ElectricCharge(),
          part.AbsoluteStrangeness(),
          part.ResonanceWidth(),
          part.DecayThresholdMass(),
          -part.Charm(),
          part.AbsoluteCharm()));

        myModel->GetTPS()->FillPdgMap();
      }
    }

    myModel->updateParticle(row);

    labelAbsoluteLightQuarkValue->setText(QString::number(part.AbsoluteQuark()));
  }
}

void ListEditorTab::StableEdited()
{
  int row = getCurrentRow();
  if (row >= 0) {
    int index = myModel->GetRowToParticle()[row];
    myModel->GetTPS()->Particle(index).SetStable(checkBoxStable->isChecked());

    myModel->updateParticle(row);
  }
}

void ListEditorTab::editDecaysDoubleClick(const QModelIndex & index) {
  int row = index.row();
  if (row >= 0) {
    decayseditor dialog(this, &myModel->GetTPS()->Particle(myModel->GetRowToParticle()[row]), myModel->GetTPS());
    dialog.exec();
  }
}

void ListEditorTab::resetParticle()
{
  int row = getCurrentRow();
  int index = myModel->GetRowToParticle()[row];
  if (row >= 0) {
    int pdgid = myModel->GetTPS()->Particles()[index].PdgId();
    if (myModel->GetTPSorig()->PdgToId(pdgid) == -1) {
      //QMessageBox::show(QString("No particle with PDG ID %1 in original database!").arg(pdgid));
      QMessageBox::warning(this, "", QString("No particle with PDG ID %1 in the original list!").arg(pdgid));
    }
    else {
      myModel->GetTPS()->Particle(index) = myModel->GetTPSorig()->Particle(myModel->GetTPSorig()->PdgToId(pdgid));
      if (!myModel->GetTPS()->Particle(index).IsNeutral()) {
        int tind     = myModel->GetTPS()->PdgToId(-pdgid);
        int tindOrig = myModel->GetTPSorig()->PdgToId(-pdgid);
        if (tind != -1 && tindOrig != -1)
          myModel->GetTPS()->Particle(tind) = myModel->GetTPSorig()->Particle(tindOrig);
      }
    }

    myModel->updateParticle(row);

    changedRow();
  }
}

void ListEditorTab::resetAll() {
  myModel->reset();
  //leDatabase->setText(path);
  tableParticles->resizeColumnsToContents();
}

void ListEditorTab::addParticle() {
  if (myModel->GetTPS()->PdgToId(0) == -1) {
    myModel->addParticle(ThermalParticle());
    tableParticles->selectRow(myModel->GetRowToParticle().size() - 1);
  }
  else QMessageBox::warning(this, "", "New particle already added!");
}

void ListEditorTab::removeParticle() {
  int row = getCurrentRow();
  if (row >= 0) {
    myModel->removeParticle(row);
  }
}

void ListEditorTab::setParticleList(thermalfist::ThermalParticleSystem * TPS)
{
  myModel->setTPS(TPS);
  tableParticles->resizeColumnsToContents();
}

namespace {
  bool cmpBranchingRatios(const std::pair<double, double> &a, const std::pair<double, double> &b) {
    return (a.second < b.second);
  }
}

void ListEditorTab::setParticleMassThreshold(int id, int type)
{
  ThermalParticle &part = myModel->GetTPS()->Particle(id);
  
  std::vector< std::pair<double, double> > thresholds(0);
  for (int i = 0; i < part.Decays().size(); ++i) {
    double thre = 0.;
    for (int j = 0; j < part.Decays()[i].mDaughters.size(); ++j) {
      int tpdgid = part.Decays()[i].mDaughters[j];
      if (myModel->GetTPS()->PdgToId(tpdgid) != -1)
        thre += myModel->GetTPS()->Particle(myModel->GetTPS()->PdgToId(tpdgid)).Mass();
    }
    thresholds.push_back(std::make_pair(thre, part.Decays()[i].mBratio));
  }

  if (thresholds.size() == 0)
    return;

  double thr = 0.;

  if (type == 0) // Lowest channel
  {
    thr = (*std::min_element(thresholds.begin(), thresholds.end())).first;
  }
  else if (type == 1) // Dominant channel
  {
    thr = (*std::max_element(thresholds.begin(), thresholds.end(), cmpBranchingRatios)).first;
  }
  else // Weighted sum of channels
  {
    thr = 0.;
    double brsum = 0.;
    for (int j = 0; j < thresholds.size(); ++j) {
      thr += thresholds[j].second * thresholds[j].first;
      brsum += thresholds[j].second;
    }
    thr *= 1. / brsum;
  }

  part.SetDecayThresholdMass(thr);
}

void ListEditorTab::setThreshold()
{
  QStringList items;
  QString lowest = tr("Lowest decay channel");
  QString dominant = tr("Dominant decay channel");
  QString weighted = tr("Weighted sum of channels");
  items << tr("Lowest decay channel") << tr("Dominant decay channel") << tr("Weighted sum of channels");

  bool ok;
  QString item = QInputDialog::getItem(this, tr("Set threshold mass"),
    tr("Set threshold mass from:"), items, 0, false, &ok);
  if (ok && !item.isEmpty()) {

    int row = getCurrentRow();
    if (row >= 0) {
      int index = myModel->GetRowToParticle()[row];
      if (item == lowest)
        setParticleMassThreshold(index, 0);
      else if (item == dominant)
        setParticleMassThreshold(index, 1);
      else
        setParticleMassThreshold(index, 2);

      // Check if it has antiparticle
      if (!myModel->GetTPS()->Particles()[index].IsNeutral()) {
        int pdgid = myModel->GetTPS()->Particles()[index].PdgId();
        if (myModel->GetTPS()->PdgToId(-pdgid) != -1) {
          int indexa = myModel->GetTPS()->PdgToId(-pdgid);
          if (item == lowest)
            setParticleMassThreshold(indexa, 0);
          else if (item == dominant)
            setParticleMassThreshold(indexa, 1);
          else
            setParticleMassThreshold(indexa, 2);
        }
      }

      spinThreshold->setValue(myModel->GetTPS()->Particles()[index].DecayThresholdMass());
    }
    
  }

}

void ListEditorTab::setAllThresholds()
{
  QStringList items;
  QString lowest = tr("Lowest decay channel");
  QString dominant = tr("Dominant decay channel");
  QString weighted = tr("Weighted sum of channels");
  items << tr("Lowest decay channel") << tr("Dominant decay channel") << tr("Weighted sum of channels");

  bool ok;
  QString item = QInputDialog::getItem(this, tr("Set threshold mass"),
    tr("Set threshold mass from:"), items, 0, false, &ok);

  if (ok && !item.isEmpty()) {
    for (int ind = 0; ind < myModel->GetTPS()->Particles().size(); ++ind) {
      if (item == lowest)
        setParticleMassThreshold(ind, 0);
      else if (item == dominant)
        setParticleMassThreshold(ind, 1);
      else
        setParticleMassThreshold(ind, 2);
    }

    int row = getCurrentRow();
    if (row >= 0) {
      int index = myModel->GetRowToParticle()[row];
      spinThreshold->setValue(myModel->GetTPS()->Particles()[index].DecayThresholdMass());
    }
  }
}

void ListEditorTab::saveToFile()
{
  QString listpathprefix = QString(ThermalFIST_INPUT_FOLDER) + "/list/list.dat";
  if (currentPath.size() != 0) {
    listpathprefix = currentPath;
  }
  QString path = QFileDialog::getSaveFileName(this, tr("Save particle list as"), listpathprefix);
  if (path.length()>0)
  {
    myModel->GetTPS()->WriteTableToFile(path.toStdString());

    listpathprefix = QFileInfo(path).absolutePath() + "/decays.dat";

    path = QFileDialog::getSaveFileName(this, tr("Save particle decays list as"), listpathprefix);

    if (path.length() > 0)
    {
      myModel->GetTPS()->WriteDecaysToFile(path.toStdString());
      currentPath = path;
    }
  }

  if (myModel->haveChanges())
    applyChanges();

}

void ListEditorTab::loadMassesWidthsFromPdg()
{
  QString listpathprefix = QString(ThermalFIST_INPUT_FOLDER) + "/list/mass_width_2020.mcd";
  if (currentPath.size() != 0) {
    listpathprefix = QFileInfo(currentPath).absolutePath() + "/pdglisting/";
  }
  QString path = QFileDialog::getOpenFileName(this, tr("Open file with particle list from PDG website"), 
    listpathprefix, 
    tr("PDG file (*.mcd)"));
  if (path.length() > 0)
  {
    std::ifstream fin(path.toStdString());
    if (fin.is_open()) {
      char cc[2000];
      while (!fin.eof()) {
        fin.getline(cc, 2000);
        std::string tmp = std::string(cc);
        if (tmp.size() == 0 || tmp[0] == '*')
          continue;
        std::string strpdg = tmp.substr(0, 34);
        std::string strvals = tmp.substr(34);
        std::string strnam = tmp.substr(107);
        if (strvals.size() == 0)
          continue;

        std::istringstream isspdg(strpdg);
        std::vector<long long> pdgs;
        long long tpdg;
        while (isspdg >> tpdg) {
          pdgs.push_back(tpdg);
        }
        if (pdgs.size() == 0)
          continue;

        std::istringstream issvals(strvals);
        double mass = 0., width = 0., tmpval = 0.;
        if (!(issvals >> mass))
          continue;

        if (!(issvals >> tmpval >> tmpval >> width))
          width = -1;

        std::istringstream issnam(strnam);
        std::string pname = "";
        if (!(issnam >> pname))
          pname = "";

        for (int i = 0; i < pdgs.size(); ++i) {
          tpdg = pdgs[i];
          int tid = myModel->GetTPS()->PdgToId(tpdg);
          if (tid != -1) {
            ThermalParticle& part = myModel->GetTPS()->Particle(tid);
            part.SetMass(mass);
            if (width > 0.)
              part.SetResonanceWidth(width);

            if (pname != "") {
              std::string suff = "";
              std::string origname = part.Name();
              for (int ic = origname.size() - 1; ic >= 0; ic--) {
                if (origname[ic] != '+' && origname[ic] != '0' && origname[ic] != '-')
                  break;
                else
                  suff += origname[ic];
              }
              reverse(suff.begin(), suff.end());
              part.SetName(pname + suff);
            }
          }
          tid = myModel->GetTPS()->PdgToId(-tpdg);
          if (tid != -1) {
            ThermalParticle& part = myModel->GetTPS()->Particle(tid);
            part.SetMass(mass);
            if (width > 0.)
              part.SetResonanceWidth(width);
            if (pname != "") {
              part.SetName("anti-" + myModel->GetTPS()->ParticleByPDG(tpdg).Name());
            }
          }
        }
      }
      fin.close();

      //myModel->reset();
    }
  }
}

bool ListEditorTab::haveChangesToList()
{
  return myModel->haveChanges();
}

void ListEditorTab::resetTPS()
{
  model->ChangeTPS(model->TPS());
  myModel->setTPS(model->TPS());
  tableParticles->resizeColumnsToContents();
}
