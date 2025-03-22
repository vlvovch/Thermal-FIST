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

#include "calculationtabledialog.h"


class EoSWorker : public QThread
{
  Q_OBJECT

  thermalfist::ThermalModelBase *model;
  ThermalModelConfig config;
  double Tmin, Tmax, dT;
  std::vector<double> cParams;
  int mode;
  int *currentSize;
  int *stop;

  std::vector< Thermodynamics > *paramsTD;
  std::vector< ChargesFluctuations > *paramsFl;
  std::vector<double> *varvalues;

  void run() Q_DECL_OVERRIDE;

public:
  EoSWorker(thermalfist::ThermalModelBase *mod = NULL,
      const ThermalModelConfig& cconfig = ThermalModelConfig(),
      double Tmin = 100., 
      double Tmax = 200.,
      double dT = 5.,
      const std::vector<double>& cParamVals = {0.,0.,0.},
      int mmode = 0,
      std::vector< Thermodynamics > *paramsTDo = NULL,
      std::vector< ChargesFluctuations > *paramsFlo = NULL,
      std::vector<double> *varvalueso = NULL,
      int *currentSizeo = NULL,
      int *stopo = NULL,
      QObject * parent = 0) :
  QThread(parent), Tmin(Tmin), Tmax(Tmax), dT(dT), cParams(cParamVals), mode(mmode) {
      model = mod;
      config = cconfig;
      paramsTD = paramsTDo;
      paramsFl = paramsFlo;
      varvalues = varvalueso;
      currentSize = currentSizeo;
      stop = stopo;

      if (mode == 0) {
        mod->SetBaryonChemicalPotential(cParams[0] * 1.e-3);
        mod->SetElectricChemicalPotential(cParams[1] * 1.e-3);
        mod->SetStrangenessChemicalPotential(cParams[2] * 1.e-3);
        mod->FillChemicalPotentials();
      }
      else if (mode == 1) {
        mod->SetBaryonChemicalPotential(0.);
        mod->FillChemicalPotentials();
      }
      else if (mode == 2) {
        mod->SetTemperature(cParams[0] * 1.e-3);
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
    QCheckBox *CBxAxis;
    QComboBox *comboQuantity2;

    QComboBox *comboParticle;
    QComboBox *comboParticle2;
    QLabel    *labelParticle, *labelParticle2;
    QComboBox *comboFeeddown;
    QComboBox *comboFeeddown2;
    QLabel    *labelFeeddown, *labelFeeddown2;

    QCheckBox *CBflipAxes;

    QLabel *labelTMin, *labelTMax, *labeldT;
    QLabel *labelEMin, *labelEMax, *labeldE;
    QLabel *labelmuB, *labelmuQ, *labelmuS;
    QLabel *labelTaux;

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
    QDoubleSpinBox *spinmuB, *spinmuQ, *spinmuS;
    QDoubleSpinBox *spinEMin, *spinEMax, *spindE;
    QDoubleSpinBox *spinTaux;
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

    // Table for the calculation results table widget
    CalculationTable calcTable;

    QPushButton *buttonEoSTable;

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

    void showEoSTable();

private:
  std::vector<double> getValues(int index, int num = 0);
  std::vector<double> getValuesRatio(int index, int index2);
  QString getParameterName() const;
  void readLatticeData();
  void recomputeCalcTable();
};


#endif // EQUATIONOFSTATETAB_H
