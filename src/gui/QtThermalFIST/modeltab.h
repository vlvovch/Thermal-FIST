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
//#include "fittoexperimenttab.h"

class TableModel;
class FitToExperimentTab;

class ModelTab : public QWidget
{
    Q_OBJECT

    QTableView *tableParticles;

		QRadioButton *radIdeal, *radEVD, *radEVCRS, *radQVDW;
		QRadioButton *radGCE, *radCE, *radSCE;

    QRadioButton *radioBoltz, *radioQuant;
		QCheckBox *CBBoseOnly, *CBPionsOnly;
		QCheckBox *CBQuadratures;

    QCheckBox *checkOnlyStable;
    QPushButton *buttonResults;
		QPushButton *labelValid;
    QLabel *labelHint;

		QLabel *labelmuS, *labelmuC;
    QDoubleSpinBox *spinTemperature, *spinmuB, *spingammaq, *spingammaS, *spinmuS, *spinmuQ, *spinmuC, *spinVolumeR;
		QDoubleSpinBox *spinVolumeRSC;
		QLabel *labelB, *labelQ, *labelS, *labelC;
    QSpinBox *spinB, *spinS, *spinQ, *spinC;
    QDoubleSpinBox *spinQBRatio;
    QDoubleSpinBox *spinRadius;
		

		QCheckBox *checkFixMuQ, *checkFixMuS, *checkFixMuC;

		QComboBox *comboWidth;
    QCheckBox *checkBratio;
		QCheckBox *checkFluctuations;

    QCheckBox *checkOMP;

    QRadioButton *radioUniform, *radioBaglike, *radioMesons, *radioCustomEV;
		QString strEVPath;


    QPushButton *buttonCalculate;
    QPushButton *buttonWriteToFile;
    QPushButton *buttonBenchmark;
    QPushButton *buttonCalculateFitted;

    TableModel *myModel;

    thermalfist::ThermalModelBase *model;

		ChargesFluctuations flucts;

    QTextEdit *teDebug;

    const FitToExperimentTab *tabFit;

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
		void loadEVFromFile();
    void showResults();
    void setModel(thermalfist::ThermalModelBase *model);
    void modelChanged();
		void changeVolumeRSC(double);
		void showValidityCheckLog();
		void computeHigherOrderFluctuations();
public:
    void updateModel();
    void resetTPS();
};

#endif // MODELTAB_H
