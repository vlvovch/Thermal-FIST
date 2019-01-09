/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef THERMALFITPLOTS_H
#define THERMALFITPLOTS_H

#include "qcustomplot.h"
#include "HRGFit/ThermalModelFit.h"

std::string BeautifyName(std::string name);

class PlotDialog : public QDialog
{
  Q_OBJECT

  QCustomPlot *plot;

public:
  explicit  PlotDialog(QWidget *parent = 0, QCustomPlot *plotin = 0);

private slots:
  void contextMenuRequest(QPoint pos);
  void saveAsPdf();
  void saveAsPng();
  void resetScales();

public slots :
  //void replot();
};

class YieldsPlot : public QCustomPlot
{
  Q_OBJECT

  QVector<double> dataValues, dataErrors, modelValues, dataX;
  QVector<QString> names;
  QVector<QCPItemText*> labels;
public:
  explicit YieldsPlot(QWidget *parent = 0, thermalfist::ThermalModelFit *fit = NULL);

protected:
  void resizeEvent(QResizeEvent *event) override;

};

class DeviationsPlot : public QCustomPlot
{
  Q_OBJECT

  QVector<double> deviationsValues, dataX;
  QVector<QString> names;
  QVector<QCPItemText*> labels;
public:
  explicit DeviationsPlot(QWidget *parent = 0, thermalfist::ThermalModelFit *fit = NULL);

protected:
  void resizeEvent(QResizeEvent *event) override;
};

class DataModelPlot : public QCustomPlot
{
  Q_OBJECT

  QVector<double> datamodelValues, datamodelErrors, dataX;
  QVector<QString> names;
  QVector<QCPItemText*> labels;
public:
  explicit DataModelPlot(QWidget *parent = 0, thermalfist::ThermalModelFit *fit = NULL);

protected:
  void resizeEvent(QResizeEvent *event) override;
};

class DataVsModelPlot : public QCustomPlot
{
  Q_OBJECT

  QVector<int> indmap, pdgs;
  QVector<double> dataValues, dataErrors, modelValues;
  QVector<QString> names;
  QVector<QCPItemText*> labels;
public:
  explicit DataVsModelPlot(QWidget *parent = 0, thermalfist::ThermalModelFit *fit = NULL);

protected:
  void resizeEvent(QResizeEvent *event) override;
};


#endif // THERMALFITPLOTS_H
