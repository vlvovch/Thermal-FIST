/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
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
    spinPDG1 = new QSpinBox();
    spinPDG1->setMinimum(-2000000000);
    spinPDG1->setMaximum(2000000000);
    if (!quant->type) spinPDG1->setValue(quant->mult.fPDGID);
    else spinPDG1->setValue(quant->ratio.fPDGID1);
    spinPDG2 = new QSpinBox();
    spinPDG2->setMinimum(-2000000000);
    spinPDG2->setMaximum(2000000000);
    if (!quant->type) spinPDG2->setEnabled(false);
    else spinPDG2->setValue(quant->ratio.fPDGID2);
    name1 = new QLabel();
    name2 = new QLabel();
    name1->setText(QString(model->TPS()->GetNameFromPDG(spinPDG1->value()).c_str()));
    name2->setText(QString(model->TPS()->GetNameFromPDG(spinPDG2->value()).c_str()));
    layGrid->addWidget(labPDG1, 0, 0);
    layGrid->addWidget(spinPDG1, 0, 1);
    layGrid->addWidget(name1, 0, 2);
    layGrid->addWidget(labPDG2, 1, 0);
    layGrid->addWidget(spinPDG2, 1, 1);
    layGrid->addWidget(name2, 1, 2);

    connect(spinPDG1, SIGNAL(valueChanged(int)), this, SLOT(PDGChanged1(int)));
    connect(spinPDG2, SIGNAL(valueChanged(int)), this, SLOT(PDGChanged2(int)));

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
		spinPDG2->setEnabled(false);
		comboFeeddown2->setEnabled(false);
	}
	else {
		spinPDG2->setEnabled(true);
		comboFeeddown2->setEnabled(true);
	}
}

void QuantityDialog::returnOK() {
  if (!checkType->isChecked()) {
		quant->type  = FittedQuantity::Multiplicity;
    quant->mult  = ExperimentMultiplicity(spinPDG1->value(), spinValue->value(), spinError->value());
		quant->mult.fFeedDown = static_cast<Feeddown::Type>(comboFeeddown1->currentIndex());
  }
	else {
		quant->type  = FittedQuantity::Ratio;
		quant->ratio = ExperimentRatio(spinPDG1->value(), spinPDG2->value(), spinValue->value(), spinError->value());
		quant->ratio.fFeedDown1 = static_cast<Feeddown::Type>(comboFeeddown1->currentIndex());
		quant->ratio.fFeedDown2 = static_cast<Feeddown::Type>(comboFeeddown2->currentIndex());
	}
  this->accept();
}

void QuantityDialog::returnCancel() {
    this->reject();
}

void QuantityDialog::PDGChanged1(int pdgid) {
    name1->setText(QString(model->TPS()->GetNameFromPDG(pdgid).c_str()));
}

void QuantityDialog::PDGChanged2(int pdgid) {
    name2->setText(QString(model->TPS()->GetNameFromPDG(pdgid).c_str()));
}
