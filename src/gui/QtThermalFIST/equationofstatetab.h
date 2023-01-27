/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef EQUATIONOFSTATETAB_H
#define EQUATIONOFSTATETAB_H

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


class EoSWorker : public QThread
{
  Q_OBJECT

  thermalfist::ThermalModelBase *model;
  double Tmin, Tmax, dT;
  double cParam;
  int mode;
  int *currentSize;
  int *stop;

  std::vector< Thermodynamics > *paramsTD;
  std::vector< ChargesFluctuations > *paramsFl;
  std::vector<double> *varvalues;

  void run() Q_DECL_OVERRIDE;

public:
  EoSWorker(thermalfist::ThermalModelBase *mod = NULL, 
      double Tmin = 100., 
      double Tmax = 200.,
      double dT = 5.,
      double cParamVal = 0.,
      int mmode = 0,
      std::vector< Thermodynamics > *paramsTDo = NULL,
      std::vector< ChargesFluctuations > *paramsFlo = NULL,
      std::vector<double> *varvalueso = NULL,
      int *currentSizeo = NULL,
      int *stopo = NULL,
      QObject * parent = 0) :
  QThread(parent), Tmin(Tmin), Tmax(Tmax), dT(dT), cParam(cParamVal), mode(mmode) {
      model = mod;
      paramsTD = paramsTDo;
      paramsFl = paramsFlo;
      varvalues = varvalueso;
      currentSize = currentSizeo;
      stop = stopo;

      if (mode == 0) {
        mod->SetBaryonChemicalPotential(cParam * 1.e-3);
        mod->FillChemicalPotentials();
      }
      else if (mode == 1) {
        mod->SetBaryonChemicalPotential(0.);
        mod->FillChemicalPotentials();
      }
      else if (mode == 2) {
        mod->SetTemperature(cParam * 1.e-3);
      }
    }
signals:
  void calculated();
};

class EquationOfStateTab : public QWidget
{
    Q_OBJECT

    QComboBox *comboQuantity;
    QCheckBox *CBratio;
    QComboBox *comboQuantity2;

    QComboBox *comboParticle;
    QComboBox *comboParticle2;
    QLabel    *labelParticle, *labelParticle2;
    QComboBox *comboFeeddown;
    QComboBox *comboFeeddown2;
    QLabel    *labelFeeddown, *labelFeeddown2;

    QLabel *labelmuB, *labelTMin, *labelTMax, *labeldT;

    QCustomPlot *plotDependence;

    std::vector<Thermodynamics> paramsTD;
    std::vector<ChargesFluctuations> paramsFl;
    std::vector<double> varvalues;

    int fCurrentSize;
    bool fRunning;
    int fStop;

    std::map<QString, int> parammap;
    std::vector<QString> paramnames;

    QDoubleSpinBox *spinTMin, *spinTMax, *spindT;
    QDoubleSpinBox *spinmuB;
    QLabel *labelConstr;

    QPushButton *buttonCalculate;

    thermalfist::ThermalParticleSystem *TPS;
    thermalfist::ThermalModelBase *model;

    QTimer *calcTimer;

    QString cpath;

    ModelConfigWidget *configWidget;


    QComboBox* comboMode; // 0 - const muB, 1 - const muB/T, 2 - const T

    // Wuppertal-Budapest lattice data
    QVector< QVector<double> > dataWBx, dataWBy, dataWByerrp, dataWByerrm;
    std::map<QString, int> mapWB;

    // HotQCD lattice data
    QVector< QVector<double> > dataHotQCDx, dataHotQCDy, dataHotQCDyerrp, dataHotQCDyerrm;
    std::map<QString, int> mapHotQCD;

public:
    EquationOfStateTab(QWidget *parent = 0, thermalfist::ThermalModelBase *model=NULL);
    ~EquationOfStateTab();

    ThermalModelConfig getConfig();
signals:

public slots:
    void calculate();
    void replot();
    void finalize();
    void modelChanged();
    void resetTPS();
    void plotLatticeData();
    void fillParticleLists();
    void contextMenuRequest(QPoint pos);
    void saveAsPdf();
    void saveAsPng();
    void saveAsAscii();
    // saveAs type: 0 - pdf, 1 - png, 2 - ascii data
    void saveAs(int type);

private:
  std::vector<double> getValues(int index, int num = 0);
  std::vector<double> getValuesRatio(int index, int index2);
  void readLatticeData();
};


#endif // EQUATIONOFSTATETAB_H
