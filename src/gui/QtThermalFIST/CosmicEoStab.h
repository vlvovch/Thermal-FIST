/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef COSMICEOSTAB_H
#define COSMICEOSTAB_H

#include <QWidget>
#include <QPushButton>
#include <QSpinBox>
#include <QDoubleSpinBox>
#include <QCheckBox>
#include <QComboBox>
#include <QRadioButton>
#include <QThread>
#include <QTimer>
#include <map>

#include "HRGBase/ThermalModelBase.h"
#include "qcustomplot.h"
#include "BaseStructures.h"
#include "configwidgets.h"
#include "HRGBase/xMath.h"
#include "CosmicEos/CosmicEoS.h"

#include "calculationtabledialog.h"

class CosmicEoSWorker : public QThread
{
  Q_OBJECT

  thermalfist::CosmicEoS *cosmos;
  double Tmin, Tmax, dT;
  std::vector<double> cParams;
  std::vector<double> chems;
  int mode;
  int *currentSize;
  int *stop;

  std::vector< ThermodynamicsCosmic > *paramsTD;
  std::vector< Thermodynamics > *paramsTDHRG;
  std::vector<double> *varvalues;

  void run() Q_DECL_OVERRIDE;

public:
  CosmicEoSWorker(thermalfist::CosmicEoS *mod = NULL,
      double Tmin = 10.,
      double Tmax = 200.,
      double dT = 2.,
      std::vector<double> cParamVals = {8.6e-11,0.,0.,0.,0.},
      std::vector< ThermodynamicsCosmic > *paramsTDo = NULL,
      std::vector< Thermodynamics > *paramsTDHRGo = NULL,
      std::vector<double> *varvalueso = NULL,
      int *currentSizeo = NULL,
      int *stopo = NULL,
      QObject * parent = 0) :
  QThread(parent), Tmin(Tmin), Tmax(Tmax), dT(dT), cParams(cParamVals) {
      cosmos = mod;
      paramsTD = paramsTDo;
      paramsTDHRG = paramsTDHRGo;
      varvalues = varvalueso;
      currentSize = currentSizeo;
      stop = stopo;

      chems = std::vector<double>({ 0.700, -1.e-7, -1.e-7, -1.e-7, -1.e-7 });

      cosmos->SetAsymmetries(cParamVals);
    }
signals:
  void calculated();
};

class CosmicEoSTab : public QWidget
{
    Q_OBJECT

    QComboBox *comboQuantity;
    QCheckBox *CBratio;
    QCheckBox *CBxAxis;
    QComboBox *comboQuantity2;

    QComboBox *comboParticle;
    QComboBox *comboParticle2;
    QLabel    *labelParticle, *labelParticle2;
    QComboBox *comboFeeddown;
    QComboBox *comboFeeddown2;
    QLabel    *labelFeeddown, *labelFeeddown2;

    QCheckBox *CBflipAxes;

    QLabel *labelmuB, *labelTMin, *labelTMax, *labeldT;

    QCustomPlot *plotDependence;

    std::vector<ThermodynamicsCosmic> paramsTD;
    std::vector<Thermodynamics> paramsTDHRG;
//    std::vector<ChargesFluctuations> paramsFl;
    std::vector<double> varvalues;

    int fCurrentSize;
    bool fRunning;
    int fStop;

    std::map<QString, int> parammap;
    std::vector<QString> paramnames;

    QDoubleSpinBox *spinTMin, *spinTMax, *spindT;
//    QDoubleSpinBox *spinmuB;
//    QLabel *labelConstr;

    QDoubleSpinBox *spinnB, *spinnQ, *spinnE, *spinnMu, *spinnTau;

    QPushButton *buttonCalculate;

    thermalfist::ThermalParticleSystem *TPS;
    thermalfist::ThermalModelBase *model;
    thermalfist::CosmicEoS *cosmos;

    QTimer *calcTimer;

    QString cpath;

    ModelConfigWidget *configWidget;

    // Table for the calculation results table widget
    CalculationTable calcTable;

    QPushButton *buttonEoSTable;

public:
    CosmicEoSTab(QWidget *parent = 0, thermalfist::ThermalModelBase *model=NULL);
    ~CosmicEoSTab();

    ThermalModelConfig getConfig();
signals:

public slots:
    void calculate();
    void replot();
    void finalize();
    void modelChanged();
    void resetTPS();
    void fillParticleLists();
    void contextMenuRequest(QPoint pos);
    void saveAsPdf();
    void saveAsPng();
    void saveAsAscii();
    // saveAs type: 0 - pdf, 1 - png, 2 - ascii data
    void saveAs(int type);

    void showEoSTable();
private:
  std::vector<double> getValues(int index, int num = 0);
  std::vector<double> getValuesRatio(int index, int index2);
  QString getParameterName() const;
  void recomputeCalcTable();
};


#endif // CosmicEoStab_H
