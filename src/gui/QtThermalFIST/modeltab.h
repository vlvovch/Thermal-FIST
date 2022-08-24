/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef MODELTAB_H
#define MODELTAB_H

#include <QWidget>
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

#include "HRGBase/ThermalModelBase.h"
#include "HRGFit/ThermalModelFit.h"
#include "BaseStructures.h"
#include "configwidgets.h"


class TableModel;
class FitToExperimentTab;

class ModelTab : public QWidget
{
    Q_OBJECT

    QTableView *tableParticles;

    QCheckBox *checkOnlyStable;
    QPushButton *buttonResults;
    QPushButton *buttonCorrelations;
    QPushButton *labelValid;
    QLabel *labelHint;

    QLabel *labelmuS, *labelmuC;
    QLabel *labelgammaS, *labelgammaC;
    QDoubleSpinBox *spinTemperature, *spinmuB, *spingammaq, *spingammaS, *spingammaC, *spinmuS, *spinmuQ, *spinmuC, *spinVolumeR;
    QDoubleSpinBox *spinVolumeRSC;
    QLabel *labelVolumeVal;
    QLabel *labelB, *labelQ, *labelS, *labelC;
    QSpinBox *spinB, *spinS, *spinQ, *spinC;

    QCheckBox *checkFluctuations;
    QCheckBox *checkMuInitials;

    QCheckBox *checkOMP;

    QPushButton *buttonCalculate;
    QPushButton *buttonWriteToFile;
    QPushButton *buttonBenchmark;
    QPushButton *buttonCalculateFitted;

    TableModel *myModel;

    thermalfist::ThermalModelBase *model;

    ChargesFluctuations flucts;

    QTextEdit *teDebug;

    const FitToExperimentTab *tabFit;

    ModelConfigWidget *configWidget;

    int getCurrentRow();

public:
    ModelTab(QWidget *parent = 0, thermalfist::ThermalModelBase *model=NULL);
    ~ModelTab();
    ThermalModelConfig getConfig();
    ThermalModelConfig getConfigFromFit(thermalfist::ThermalModelFit *fit, const ThermalModelConfig & configfit);
    void setFitTab(const FitToExperimentTab *tab) { tabFit = tab; }
    void updateControlsWithConfig(const ThermalModelConfig & config);
private slots:
    void changedRow();
    void performCalculation(const ThermalModelConfig & config);
    void calculate();
    void calculateFitted();
    void writetofile();
    void benchmark();
    void particleInfoDoubleClick(const QModelIndex &);
    void switchStability(bool);
    void showResults();
    void showCorrelations();
    void setModel(thermalfist::ThermalModelBase *model);
    void modelChanged();
    void changeVolumeRSC(double);
    void showValidityCheckLog();
    void computeHigherOrderFluctuations();
public:
    void updateModel();
    void resetTPS();
    void updateFontSizes();
};

#endif // MODELTAB_H
