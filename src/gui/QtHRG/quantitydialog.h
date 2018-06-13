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

#include "tablemodel.h"

class ThermalModelBase;
class FittedQuantity;

class QuantityDialog : public QDialog
{
    Q_OBJECT

    //ThermalParticle *fParticle;
    //ThermalParticleSystem *fTPS;

    ThermalModelBase *model;
    FittedQuantity *quant;

    int pid;

    QCheckBox *checkType;
    QSpinBox *spinPDG1;
    QSpinBox *spinPDG2;
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
    explicit  QuantityDialog(QWidget *parent = 0, ThermalModelBase *mod = NULL, FittedQuantity *quantity = NULL);

signals:

public slots:
    void changeType(bool flag);
    void returnOK();
    void returnCancel();
    void PDGChanged1(int pdgid);
    void PDGChanged2(int pdgid);
};

#endif // QUANTITYDIALOG_H
