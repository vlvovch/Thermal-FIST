/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef QUANTITYDIALOG_H
#define QUANTITYDIALOG_H

#include <QDialog>
#include <QPushButton>
#include <QSpinBox>
#include <QDoubleSpinBox>
#include <QCheckBox>
#include <QDialogButtonBox>
#include <QLabel>
#include <QComboBox>
#include <QLineEdit>

#include "HRGBase/ThermalModelBase.h"
#include "HRGFit/ThermalModelFit.h"
#include "tablemodel.h"


class QuantityDialog : public QDialog
{
    Q_OBJECT

    thermalfist::ThermalModelBase *model;
    thermalfist::FittedQuantity *quant;

    int pid;

    QCheckBox *checkType;
    QLineEdit *lePDG1;
    QLineEdit *lePDG2;
    QLabel *name1, *name2;
    QDoubleSpinBox *spinValue;
    QDoubleSpinBox *spinError;

    QComboBox *comboFeeddown1;
    QComboBox *comboFeeddown2;

    QPushButton *buttonOK;
    QPushButton *buttonCancel;

    QDialogButtonBox *buttonBox;

    QString GetParticleInfo();
public:
    explicit  QuantityDialog(QWidget *parent = 0, thermalfist::ThermalModelBase *mod = NULL, thermalfist::FittedQuantity *quantity = NULL);

signals:

public slots:
    void changeType(bool flag);
    void returnOK();
    void returnCancel();
    void PDGChanged1();
    void PDGChanged2();
};

#endif // QUANTITYDIALOG_H
