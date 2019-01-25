/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "quantitydialog.h"

#include <QLayout>

using namespace thermalfist;

QuantityDialog::QuantityDialog(QWidget *parent, ThermalModelBase *mod, FittedQuantity *quantity) :
    QDialog(parent), model(mod), quant(quantity)
{
    QVBoxLayout *layout = new QVBoxLayout();

    QGridLayout *layGrid = new QGridLayout();

    checkType = new QCheckBox(tr("Ratio?"));
    checkType->setChecked(quant->type);
    connect(checkType, SIGNAL(toggled(bool)), this, SLOT(changeType(bool)));
    QLabel *labPDG1 = new QLabel(tr("PDGID 1"));
    QLabel *labPDG2 = new QLabel(tr("PDGID 2"));

    lePDG1 = new QLineEdit();
    QRegExp rx("-?\\d{1,15}");
    QValidator *validator1 = new QRegExpValidator(rx, this);
    lePDG1->setValidator(validator1);
    lePDG1->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    if (!quant->type) lePDG1->setText(QString::number(quant->mult.fPDGID));
    else lePDG1->setText(QString::number(quant->ratio.fPDGID1));

    lePDG2 = new QLineEdit();
    QValidator *validator2 = new QRegExpValidator(rx, this);
    lePDG2->setValidator(validator2);
    lePDG2->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    if (!quant->type) lePDG2->setEnabled(false);
    else lePDG2->setText(QString::number(quant->ratio.fPDGID2));

    name1 = new QLabel();
    name2 = new QLabel();
    name1->setText(QString(model->TPS()->GetNameFromPDG(lePDG1->text().toLongLong()).c_str()));
    name2->setText(QString(model->TPS()->GetNameFromPDG(lePDG2->text().toLongLong()).c_str()));
    layGrid->addWidget(labPDG1, 0, 0);
    layGrid->addWidget(lePDG1, 0, 1);
    layGrid->addWidget(name1, 0, 2);
    layGrid->addWidget(labPDG2, 1, 0);
    layGrid->addWidget(lePDG2, 1, 1);
    layGrid->addWidget(name2, 1, 2);

    connect(lePDG1, SIGNAL(textChanged(QString)), this, SLOT(PDGChanged1()));
    connect(lePDG2, SIGNAL(textChanged(QString)), this, SLOT(PDGChanged2()));

    QLabel *labFeeddown1 = new QLabel(tr("Feeddown 1"));
    comboFeeddown1 = new QComboBox();
    comboFeeddown1->addItem(tr("Primordial"));
    comboFeeddown1->addItem(tr("Stability flags"));
    comboFeeddown1->addItem(tr("Strong+EM+weak decays"));
    comboFeeddown1->addItem(tr("Strong+EM decays"));
    comboFeeddown1->addItem(tr("Strong decays"));
    QLabel *labFeeddown2 = new QLabel(tr("Feeddown 2"));
    comboFeeddown2 = new QComboBox();
    comboFeeddown2->addItem(tr("Primordial"));
    comboFeeddown2->addItem(tr("Stability flags"));
    comboFeeddown2->addItem(tr("Strong+EM+weak decays"));
    comboFeeddown2->addItem(tr("Strong+EM decays"));
    comboFeeddown2->addItem(tr("Strong decays"));

    if (!quant->type) {
      comboFeeddown1->setCurrentIndex(quant->mult.fFeedDown);
      comboFeeddown2->setCurrentIndex(1);
    }
    else {
      comboFeeddown1->setCurrentIndex(quant->ratio.fFeedDown1);
      comboFeeddown2->setCurrentIndex(quant->ratio.fFeedDown2);
    }

    if (!quant->type) comboFeeddown2->setEnabled(false);

    layGrid->addWidget(labFeeddown1, 2, 0);
    layGrid->addWidget(comboFeeddown1, 2, 1);
    layGrid->addWidget(labFeeddown2, 3, 0);
    layGrid->addWidget(comboFeeddown2, 3, 1);
    

    QLabel *labVal = new QLabel(tr("Value"));
    spinValue = new QDoubleSpinBox();
    spinValue->setMaximum(100000);
    spinValue->setDecimals(8);
    if (!quant->type) spinValue->setValue(quant->mult.fValue);
    else spinValue->setValue(quant->ratio.fValue);

    QLabel *labErr = new QLabel(tr("Error"));
    spinError = new QDoubleSpinBox();
    spinError->setDecimals(8);
    if (!quant->type) spinError->setValue(quant->mult.fError);
    else spinError->setValue(quant->ratio.fError);

    

    layGrid->addWidget(labVal, 4, 0);
    layGrid->addWidget(spinValue, 4, 1);
    layGrid->addWidget(labErr, 5, 0);
    layGrid->addWidget(spinError, 5, 1);

    buttonOK = new QPushButton(tr("OK"));
    buttonOK->setDefault(true);
    connect(buttonOK, SIGNAL(clicked()), this, SLOT(returnOK()));

    buttonCancel = new QPushButton(tr("Cancel"));
    connect(buttonCancel, SIGNAL(clicked()), this, SLOT(returnCancel()));

    buttonBox = new QDialogButtonBox(Qt::Horizontal);
    buttonBox->addButton(buttonOK, QDialogButtonBox::ActionRole);
    buttonBox->addButton(buttonCancel, QDialogButtonBox::ActionRole);


    layout->addWidget(checkType);
    layout->addLayout(layGrid);
    layout->addWidget(buttonBox,0, Qt::AlignHCenter);

    setLayout(layout);

}

void QuantityDialog::changeType(bool flag) {
  if (!flag) {
    lePDG2->setEnabled(false);
    comboFeeddown2->setEnabled(false);
  }
  else {
    lePDG2->setEnabled(true);
    comboFeeddown2->setEnabled(true);
  }
}

void QuantityDialog::returnOK() {
  if (!checkType->isChecked()) {
    quant->type  = FittedQuantity::Multiplicity;
    quant->mult  = ExperimentMultiplicity(lePDG1->text().toLongLong(), spinValue->value(), spinError->value());
    quant->mult.fFeedDown = static_cast<Feeddown::Type>(comboFeeddown1->currentIndex());
  }
  else {
    quant->type  = FittedQuantity::Ratio;
    quant->ratio = ExperimentRatio(lePDG1->text().toLongLong(), lePDG2->text().toLongLong(), spinValue->value(), spinError->value());
    quant->ratio.fFeedDown1 = static_cast<Feeddown::Type>(comboFeeddown1->currentIndex());
    quant->ratio.fFeedDown2 = static_cast<Feeddown::Type>(comboFeeddown2->currentIndex());
  }
  this->accept();
}

void QuantityDialog::returnCancel() {
    this->reject();
}

void QuantityDialog::PDGChanged1() {
  long long pdgid = lePDG1->text().toLongLong();
  name1->setText(QString(model->TPS()->GetNameFromPDG(pdgid).c_str()));
}

void QuantityDialog::PDGChanged2() {
  long long pdgid = lePDG2->text().toLongLong();
  name2->setText(QString(model->TPS()->GetNameFromPDG(pdgid).c_str()));
}
