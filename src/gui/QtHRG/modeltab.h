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

#include "BaseStructures.h"


class TableModel;
class ThermalModelBase;
//class ThermalParticleSystem;

class ModelTab : public QWidget
{
    Q_OBJECT

    QTableView *tableParticles;

    //QRadioButton *radioIdeal, *radioEVMF, *radioEVMulti, *radioCE, *radioIdealCanonStrangeness, *radioIdealCanonCharm, *radioVDWHRG;
		//QRadioButton *radioEVCanonStrangeness, *radioVDWCanonStrangeness;

		//QComboBox *comboModel;

		QRadioButton *radIdeal, *radEVD, *radEVCRS, *radQVDW;
		QRadioButton *radGCE, *radCE, *radSCE;

    QRadioButton *radioBoltz, *radioQuant;
		QCheckBox *CBBoseOnly, *CBPionsOnly;
		QCheckBox *CBQuadratures;

    QCheckBox *checkOnlyStable;
    QPushButton *buttonResults;
		//QLabel *labelValid;
		QPushButton *labelValid;

		QLabel *labelmuS, *labelmuC;
    QDoubleSpinBox *spinTemperature, *spinmuB, *spingammaq, *spingammaS, *spinmuS, *spinmuQ, *spinmuC, *spinVolumeR;
		QDoubleSpinBox *spinVolumeRSC;
    QSpinBox *spinB, *spinS, *spinQ;
    QDoubleSpinBox *spinQBRatio;
    QDoubleSpinBox *spinRadius;
		

		QCheckBox *checkFixMuQ, *checkFixMuS, *checkFixMuC;

    //QCheckBox *checkFiniteWidth;
		QComboBox *comboWidth;
    QCheckBox *checkBratio;
		QCheckBox *checkFluctuations;

    QCheckBox *checkOMP;

    QRadioButton *radioUniform, *radioBaglike, *radioMesons, *radioCustomEV;
		QString strEVPath;
    //QCheckBox *checkMesons;


    QPushButton *buttonCalculate;
    QPushButton *buttonWriteToFile;
    QPushButton *buttonBenchmark;

    TableModel *myModel;

    //ThermalParticleSystem *TPS;
    ThermalModelBase *model;

		ChargesFluctuations flucts;

    QTextEdit *teDebug;

    int getCurrentRow();

public:
    ModelTab(QWidget *parent = 0, ThermalModelBase *model=NULL);
    ~ModelTab();
		ThermalModelConfig getConfig();
private slots:
    void changedRow();
		void performCalculation(const ThermalModelConfig & config);
    void calculate();
    void writetofile();
    void benchmark();
    void particleInfoDoubleClick(const QModelIndex &);
    void switchStability(bool);
		void loadEVFromFile();
    void showResults();
    void setModel(ThermalModelBase *model);
    void modelChanged();
		void changeVolumeRSC(double);
		void showValidityCheckLog();
		void computeHigherOrderFluctuations();
public:
    void updateModel();
    void resetTPS();
};

#endif // MODELTAB_H
