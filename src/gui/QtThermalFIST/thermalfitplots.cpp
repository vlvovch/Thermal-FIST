/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "thermalfitplots.h"

#include <QDebug>

using namespace thermalfist;


PlotDialog::PlotDialog(QWidget * parent, QCustomPlot * plotin) : QDialog(parent), plot(plotin)
{
  QVBoxLayout *layout = new QVBoxLayout;
  layout->addWidget(plot);

  plot->setContextMenuPolicy(Qt::CustomContextMenu);
  connect(plot, SIGNAL(customContextMenuRequested(QPoint)), this, SLOT(contextMenuRequest(QPoint)));

  plot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);

  setLayout(layout);

  this->setWindowTitle(tr("Thermal fit result"));

  plot->replot();
}

void PlotDialog::saveAsPdf()
{
  QString listpathprefix = QApplication::applicationDirPath() + "/" + plot->objectName() + ".pdf";
  QString path = QFileDialog::getSaveFileName(this, tr("Save plot as pdf"), listpathprefix, "*.pdf");
  if (path.length()>0)
  {
    plot->savePdf(path, plot->width(), plot->height());
  }
}

void PlotDialog::saveAsPng()
{
  QString listpathprefix = QApplication::applicationDirPath() + "/" + plot->objectName() + ".png";
  QString path = QFileDialog::getSaveFileName(this, tr("Save plot as png"), listpathprefix, "*.png");
  if (path.length()>0)
  {
    plot->savePng(path, plot->width(), plot->height());
  }
}

void PlotDialog::resetScales()
{
  plot->rescaleAxes();
  plot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
  plot->replot();
}


void PlotDialog::contextMenuRequest(QPoint pos)
{
  QMenu *menu = new QMenu(this);
  menu->setAttribute(Qt::WA_DeleteOnClose);

  menu->addAction("Save as pdf", this, SLOT(saveAsPdf()));
  menu->addAction("Save as png", this, SLOT(saveAsPng()));
  menu->addAction("Reset scales", this, SLOT(resetScales()));

  menu->popup(plot->mapToGlobal(pos));
}

YieldsPlot::YieldsPlot(QWidget * parent, thermalfist::ThermalModelFit * fit) :
  QCustomPlot(parent)
{
  this->setObjectName(tr("FitYields"));
  
  int Npoints = fit->FittedQuantities().size();

  double tmin = 1.e12, tmax = 0.;
  
  for (int i = 0; i < Npoints; ++i) {
    dataX.push_back(0.5 + i);
    dataValues.push_back(fit->FittedQuantities()[i].Value());
    dataErrors.push_back(fit->FittedQuantities()[i].ValueError());
    if (i < fit->ModelDataSize())
      modelValues.push_back(fit->ModelData(i));
    if (fit->FittedQuantities()[i].type == 0) {
      names.push_back(BeautifyName(fit->model()->TPS()->GetNameFromPDG(fit->FittedQuantities()[i].mult.fPDGID)).c_str());
    }
    else {
      names.push_back((BeautifyName(fit->model()->TPS()->GetNameFromPDG(fit->FittedQuantities()[i].ratio.fPDGID1))
      + "/" + BeautifyName(fit->model()->TPS()->GetNameFromPDG(fit->FittedQuantities()[i].ratio.fPDGID2))).c_str());
    }
    tmin = std::min(tmin, dataValues[i]);
    //tmin = std::min(tmin, modelValues[i]);
    tmax = std::max(tmax, dataValues[i]);
  }
  
  this->xAxis->setTickLabels(false);
  //this->yAxis->setTickLabels(false);
  this->xAxis2->setVisible(true);
  this->xAxis2->setTickLabels(false);
  this->yAxis2->setVisible(true);
  this->yAxis2->setTickLabels(false);

  this->xAxis->setPadding((40 * height()) / 600);

  this->xAxis->setRange(0., 1. * Npoints);
  this->xAxis2->setRange(0., 1. * Npoints);

  this->yAxis->grid()->setVisible(false);
  this->yAxis2->grid()->setVisible(false);

  QSharedPointer<QCPAxisTickerText> textTicker(new QCPAxisTickerText);
  for (int i = 1; i < Npoints; ++i)
    textTicker->addTick(i, names[i]);

  this->xAxis->setTicker(textTicker);
  this->xAxis2->setTicker(textTicker);

  this->yAxis->setScaleType(QCPAxis::stLogarithmic);
  QSharedPointer<QCPAxisTickerLog> ticker(new QCPAxisTickerLog);
  ticker->setSubTickCount(4);
  this->yAxis->setTicker(ticker);
  this->yAxis->setRange(0.5 * tmin, 2. * tmax);
  this->yAxis2->setRange(0.5 * tmin, 2. * tmax);

  this->yAxis2->setScaleType(QCPAxis::stLogarithmic);
  this->yAxis2->setTicker(ticker);

  this->addGraph();
  this->graph(0)->setLineStyle(QCPGraph::lsNone);
  this->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, 7));
  this->graph(0)->setPen(QPen(Qt::red));
  this->graph(0)->setName("Data");
  this->graph(0)->data()->clear();
  this->graph(0)->setData(dataX, dataValues);

  QCPErrorBars *errorBars = new QCPErrorBars(this->xAxis, this->yAxis);
  errorBars->removeFromLegend();
  errorBars->setErrorType(QCPErrorBars::etValueError);
  errorBars->setPen(QPen(Qt::red));
  errorBars->setSymbolGap(0);
  errorBars->data().clear();
  errorBars->setData(dataErrors);
  errorBars->setDataPlottable(this->graph(0));

  for (int i = 0; i < fit->ModelDataSize(); ++i) {
    this->addGraph();
    this->graph(i + 1)->setName("Model");
    this->graph(i + 1)->setPen(QPen(Qt::black, 2, Qt::SolidLine));
    this->graph(i + 1)->data()->clear();
    this->graph(i + 1)->addData(0.5 + i - 0.25, fit->ModelData(i));
    this->graph(i + 1)->addData(0.5 + i + 0.25, fit->ModelData(i));
    if (i > 0)
      this->graph(i + 1)->removeFromLegend();
  }

  for (int i = 0; i < dataValues.size(); ++i) {
    labels.push_back(new QCPItemText(this));
    QCPItemText *textLabel = labels[i];
    textLabel->setLayer("overlay");
    textLabel->setClipToAxisRect(false);
    textLabel->position->setTypeX(QCPItemPosition::ptPlotCoords);
    textLabel->position->setTypeY(QCPItemPosition::ptAxisRectRatio);
    textLabel->position->setCoords(0.5 + i, 1. + 0.1 * 200 / (double)(height()));
    if (!names[i].startsWith("anti-")) {
      textLabel->setText(names[i]);
    }
    else {
      QFont tfont = textLabel->font();
      tfont.setOverline(true);
      textLabel->setFont(tfont);
      textLabel->setText(names[i].right(names[i].size() - 5));
    }
  }

  this->legend->setFont(QFont("Arial", 12));
  this->legend->setVisible(true);

  this->plotLayout()->insertRow(0);
  this->plotLayout()->addElement(0, 0, new QCPTextElement(this, "Yields", QFont("sans", font().pointSize() + 2, QFont::Bold)));
}

void YieldsPlot::resizeEvent(QResizeEvent *event) 
{
  for (int i = 0; i < labels.size(); ++i) {
    QCPItemText *textLabel = labels[i];
    textLabel->position->setCoords(0.5 + i, 1. + 0.1 * 200 / (double)(height()));
  }
  QCustomPlot::resizeEvent(event);
}

DataModelPlot::DataModelPlot(QWidget * parent, thermalfist::ThermalModelFit * fit)
{
  this->setObjectName(tr("FitDataModel"));
  
  int Npoints = fit->FittedQuantities().size();

  if (Npoints != fit->ModelDataSize())
    this->close();

  double tmin = 1.e12, tmax = 0.;

  for (int i = 0; i < Npoints; ++i) {
    dataX.push_back(0.5 + i);
    datamodelValues.push_back(fit->FittedQuantities()[i].Value() / fit->ModelData(i));
    datamodelErrors.push_back(fit->FittedQuantities()[i].ValueError() / fit->ModelData(i));
    
    if (fit->FittedQuantities()[i].type == 0) {
      names.push_back(BeautifyName(fit->model()->TPS()->GetNameFromPDG(fit->FittedQuantities()[i].mult.fPDGID)).c_str());
    }
    else {
      names.push_back((BeautifyName(fit->model()->TPS()->GetNameFromPDG(fit->FittedQuantities()[i].ratio.fPDGID1))
        + "/" + BeautifyName(fit->model()->TPS()->GetNameFromPDG(fit->FittedQuantities()[i].ratio.fPDGID2))).c_str());
    }

    tmin = std::min(tmin, datamodelValues[i] - 1.5 * datamodelErrors[i]);
    tmax = std::max(tmax, datamodelValues[i] + 1.5 * datamodelErrors[i]);
  }

  tmin = std::max(0., tmin);
  tmax = std::min(3., tmax);

  this->xAxis->setTickLabels(false);
  //this->yAxis->setTickLabels(false);
  this->xAxis2->setVisible(true);
  this->xAxis2->setTickLabels(false);
  this->yAxis2->setVisible(true);
  this->yAxis2->setTickLabels(false);

  //this->xAxis->setPadding((40 * height()) / 200);
  this->xAxis->setPadding(40);

  this->xAxis->setRange(0., 1. * Npoints);
  this->xAxis2->setRange(0., 1. * Npoints);

  this->yAxis->grid()->setVisible(false);
  this->yAxis2->grid()->setVisible(false);

  QSharedPointer<QCPAxisTickerText> textTicker(new QCPAxisTickerText);
  for (int i = 1; i < Npoints; ++i)
    textTicker->addTick(i, names[i]);

  this->xAxis->setTicker(textTicker);
  this->xAxis2->setTicker(textTicker);

  this->yAxis->setScaleType(QCPAxis::stLinear);
  QSharedPointer<QCPAxisTicker> ticker(new QCPAxisTicker);
  this->yAxis->setTicker(ticker);
  this->yAxis->setRange(tmin, tmax);
  this->yAxis2->setRange(tmin, tmax);

  this->yAxis2->setScaleType(QCPAxis::stLinear);
  this->yAxis2->setTicker(ticker);

  this->addGraph();
  this->graph(0)->setLineStyle(QCPGraph::lsNone);
  this->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, 7));
  this->graph(0)->setPen(QPen(Qt::red));
  this->graph(0)->setName("Data/Model");
  this->graph(0)->data()->clear();
  this->graph(0)->setData(dataX, datamodelValues);

  QCPErrorBars *errorBars = new QCPErrorBars(this->xAxis, this->yAxis);
  errorBars->removeFromLegend();
  errorBars->setErrorType(QCPErrorBars::etValueError);
  errorBars->setPen(QPen(Qt::red));
  errorBars->setSymbolGap(0);
  errorBars->data().clear();
  errorBars->setData(datamodelErrors);
  errorBars->setDataPlottable(this->graph(0));

  for (int i = 0; i < datamodelValues.size(); ++i) {
    labels.push_back(new QCPItemText(this));
    QCPItemText *textLabel = labels[i];
    textLabel->setLayer("overlay");
    textLabel->setClipToAxisRect(false);
    textLabel->position->setTypeX(QCPItemPosition::ptPlotCoords);
    textLabel->position->setTypeY(QCPItemPosition::ptAxisRectRatio);
    textLabel->position->setCoords(0.5 + i, 1. + 0.1 * 200 / (double)(height()));
    if (!names[i].startsWith("anti-")) {
      textLabel->setText(names[i]);
    }
    else {
      QFont tfont = textLabel->font();
      tfont.setOverline(true);
      textLabel->setFont(tfont);
      textLabel->setText(names[i].right(names[i].size() - 5));
    }
  }

  this->addGraph();
  this->graph(1)->setName("Unity");
  this->graph(1)->setPen(QPen(Qt::black, 2, Qt::SolidLine));
  this->graph(1)->data()->clear();
  this->graph(1)->addData(0., 1.);
  this->graph(1)->addData(1. * Npoints, 1.);
  this->graph(1)->removeFromLegend();

  this->plotLayout()->insertRow(0);
  this->plotLayout()->addElement(0, 0, new QCPTextElement(this, "Data/Model", QFont("sans", font().pointSize() + 2, QFont::Bold)));
  //this->legend->setFont(QFont("Arial", 12));
  //this->legend->setVisible(true);
}

void DataModelPlot::resizeEvent(QResizeEvent * event)
{
  for (int i = 0; i < labels.size(); ++i) {
    QCPItemText *textLabel = labels[i];
    textLabel->position->setCoords(0.5 + i, 1. + 0.1 * 200 / (double)(height()));
  }
  QCustomPlot::resizeEvent(event);
}

DataVsModelPlot::DataVsModelPlot(QWidget * parent, thermalfist::ThermalModelFit * fit)
{
  this->setObjectName(tr("FitDataVsModel"));
  
  int Npoints = fit->FittedQuantities().size();

  if (Npoints != fit->ModelDataSize())
    this->close();

  QVector< QPair<double, int> > sorter;
  for (int i = 0; i < Npoints; ++i) {
    sorter.push_back(QPair<double,int>(fit->ModelData(i),i));
  }

  qSort(sorter.begin(), sorter.end());

  for (int i = 0; i < Npoints; ++i) {
    indmap.push_back(sorter[i].second);
  }

  double tmin = 1.e12, tmax = 0.;

  for (int i = 0; i < Npoints; ++i) {
    dataValues.push_back(fit->FittedQuantities()[indmap[i]].Value());
    dataErrors.push_back(fit->FittedQuantities()[indmap[i]].ValueError());
    modelValues.push_back(fit->ModelData(indmap[i]));
    if (fit->FittedQuantities()[i].type == 0) {
      names.push_back(BeautifyName(fit->model()->TPS()->GetNameFromPDG(fit->FittedQuantities()[indmap[i]].mult.fPDGID)).c_str());
      pdgs.push_back(fit->FittedQuantities()[indmap[i]].mult.fPDGID);
    }
    else {
      names.push_back((BeautifyName(fit->model()->TPS()->GetNameFromPDG(fit->FittedQuantities()[indmap[i]].ratio.fPDGID1))
        + "/" + BeautifyName(fit->model()->TPS()->GetNameFromPDG(fit->FittedQuantities()[indmap[i]].ratio.fPDGID2))).c_str());
      pdgs.push_back(fit->FittedQuantities()[indmap[i]].ratio.fPDGID1);
    }
    tmin = std::min(tmin, dataValues[i]);
    tmax = std::max(tmax, dataValues[i]);
  }

  this->xAxis2->setVisible(true);
  this->xAxis2->setTickLabels(false);
  this->yAxis2->setVisible(true);
  this->yAxis2->setTickLabels(false);

  //this->xAxis->setPadding((40 * height()) / 600);

  tmin *= 0.25;
  tmax *= 4.;

  this->xAxis->setScaleType(QCPAxis::stLogarithmic);
  QSharedPointer<QCPAxisTickerLog> tickerx(new QCPAxisTickerLog);
  tickerx->setSubTickCount(4);
  this->xAxis->setTicker(tickerx);
  this->xAxis->setRange(tmin, tmax);
  this->xAxis2->setRange(tmin, tmax);

  this->xAxis2->setScaleType(QCPAxis::stLogarithmic);
  this->xAxis2->setTicker(tickerx);

  this->yAxis->grid()->setVisible(false);
  this->yAxis2->grid()->setVisible(false);

  this->yAxis->setScaleType(QCPAxis::stLogarithmic);
  QSharedPointer<QCPAxisTickerLog> ticker(new QCPAxisTickerLog);
  ticker->setSubTickCount(4);
  this->yAxis->setTicker(ticker);
  this->yAxis->setRange(tmin, tmax);
  this->yAxis2->setRange(tmin, tmax);

  this->yAxis2->setScaleType(QCPAxis::stLogarithmic);
  this->yAxis2->setTicker(ticker);

  this->xAxis->setLabel("Model");
  this->xAxis->setLabelFont(QFont("Arial", font().pointSize() + 4));
  this->yAxis->setLabel("Data");
  this->yAxis->setLabelFont(QFont("Arial", font().pointSize() + 4));

  this->addGraph();
  this->graph(0)->setLineStyle(QCPGraph::lsNone);
  this->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, 7));
  this->graph(0)->setPen(QPen(Qt::red));
  this->graph(0)->setName("Data");
  this->graph(0)->data()->clear();
  this->graph(0)->setData(modelValues, dataValues);

  QCPErrorBars *errorBars = new QCPErrorBars(this->xAxis, this->yAxis);
  errorBars->removeFromLegend();
  errorBars->setErrorType(QCPErrorBars::etValueError);
  errorBars->setPen(QPen(Qt::red));
  errorBars->setSymbolGap(0);
  errorBars->data().clear();
  errorBars->setData(dataErrors);
  errorBars->setDataPlottable(this->graph(0));
  //for (int i = 0; i < dataErrors.size(); ++i)
    //qDebug() << dataErrors[i];

  this->addGraph();
  this->graph(1)->setName("Unity");
  this->graph(1)->setPen(QPen(Qt::black, 2, Qt::SolidLine));
  this->graph(1)->data()->clear();
  this->graph(1)->addData(tmin, tmin);
  this->graph(1)->addData(tmax, tmax);
  this->graph(1)->removeFromLegend();

  for (int i = 0; i < dataValues.size(); ++i) {
    labels.push_back(new QCPItemText(this));
    QCPItemText *textLabel = labels[i];
    textLabel->setLayer("overlay");
    textLabel->setClipToAxisRect(false);
    textLabel->position->setTypeX(QCPItemPosition::ptPlotCoords);
    textLabel->position->setTypeY(QCPItemPosition::ptPlotCoords);
    double xcoord = modelValues[i];
    if (pdgs[i] < 0) {
      xcoord = pow(10., log10(xcoord) + log10(this->xAxis->range().upper / this->xAxis->range().lower) * 0.03 * 600 / width());
    }
    else {
      xcoord = pow(10., log10(xcoord) - log10(this->xAxis->range().upper / this->xAxis->range().lower) * 0.03 * 600 / width());
    }
    textLabel->position->setCoords(xcoord, pow(10., log10(dataValues[i]) + log10(this->yAxis->range().upper/this->yAxis->range().lower) * 0.04 * 600 / height()));
    if (!names[i].startsWith("anti-")) {
      textLabel->setText(names[i]);
    }
    else {
      QFont tfont = textLabel->font();
      tfont.setOverline(true);
      textLabel->setFont(tfont);
      textLabel->setText(names[i].right(names[i].size() - 5));
    }
  }
}

void DataVsModelPlot::resizeEvent(QResizeEvent * event)
{
  for (int i = 0; i < labels.size(); ++i) {
    QCPItemText *textLabel = labels[i];
    double xcoord = modelValues[i];
    if (pdgs[i] < 0) {
      xcoord = pow(10., log10(xcoord) - log10(this->xAxis->range().upper / this->xAxis->range().lower) * 0.03 * 600 / width());
    }
    else {
      xcoord = pow(10., log10(xcoord) + log10(this->xAxis->range().upper / this->xAxis->range().lower) * 0.03 * 600 / width());
    }
    textLabel->position->setCoords(xcoord, pow(10., log10(dataValues[i]) + log10(this->yAxis->range().upper / this->yAxis->range().lower) * 0.04 * 600 / height()));
  }
  QCustomPlot::resizeEvent(event);
}

std::string BeautifyName(std::string name)
{
  std::string ret = "";
  std::vector< std::string > strs_from;
  std::vector< std::string > strs_to;
  strs_from.push_back("π");
  strs_to.push_back("pi");
  strs_from.push_back("μ");
  strs_to.push_back("mu");
  strs_from.push_back("γ");
  strs_to.push_back("gamma");
  //strs_from.push_back("a");
  //strs_to.push_back("a");
  strs_from.push_back("ρ");
  strs_to.push_back("rho");
  //strs_from.push_back("f");
  //strs_to.push_back("f");
  strs_from.push_back("ω");
  strs_to.push_back("omega");
  strs_from.push_back("χ");
  strs_to.push_back("chi");
  strs_from.push_back("*");
  strs_to.push_back("*");
  strs_from.push_back("J");
  strs_to.push_back("J");
  strs_from.push_back("Ψ");
  strs_to.push_back("Psi");
  strs_from.push_back("φ");
  strs_to.push_back("phi");
  strs_from.push_back("η");
  strs_to.push_back("eta");
  strs_from.push_back("χ");
  strs_to.push_back("chi");
  strs_from.push_back("⁰");
  strs_to.push_back("0");
  strs_from.push_back("₀");
  strs_to.push_back("0");
  strs_from.push_back("Υ");
  strs_to.push_back("Y");
  strs_from.push_back("N");
  strs_to.push_back("N");
  strs_from.push_back("Δ");
  strs_to.push_back("Delta");
  strs_from.push_back("Σ");
  strs_to.push_back("Sigma");
  strs_from.push_back("Σ");
  strs_to.push_back("Sigma");
  strs_from.push_back("Λ");
  strs_to.push_back("Lambda");
  strs_from.push_back("Ξ");
  strs_to.push_back("Ksi");
  strs_from.push_back("Ω");
  strs_to.push_back("Omega");
  strs_from.push_back("₁");
  strs_to.push_back("1");
  strs_from.push_back("₂");
  strs_to.push_back("2");
  strs_from.push_back("₃");
  strs_to.push_back("3");
  strs_from.push_back("₄");
  strs_to.push_back("4");
  strs_from.push_back("⁺");
  strs_to.push_back("+");
  strs_from.push_back("⁻");
  strs_to.push_back("-");
  strs_from.push_back("⁺⁺");
  strs_to.push_back("++");
  strs_from.push_back("σ");
  strs_to.push_back("sigma");
  strs_from.push_back("K");
  strs_to.push_back("K");
  strs_from.push_back("e");
  strs_to.push_back("e");

  // Light nuclei
  strs_from.push_back("d");
  strs_to.push_back("Deuteron");
  strs_from.push_back("Λn");
  strs_to.push_back("LambdaNeutron");
  strs_from.push_back("Λp");
  strs_to.push_back("LambdaProton");
  strs_from.push_back("ΛΛ");
  strs_to.push_back("DiLambda");
  strs_from.push_back("³He");
  strs_to.push_back("HE3");
  strs_from.push_back("³H");
  strs_to.push_back("Triton");
  strs_from.push_back("³HΛ");
  strs_to.push_back("HyperTriton");
  strs_from.push_back("⁴He");
  strs_to.push_back("Alpha");
  
  for (int ind = 0; ind < name.size(); ++ind) {
    bool fl = false;

    if (ind == 0 && name.size() >= 5 && name.substr(0, 5) == "anti-") {
      fl = true;
      ret += name.substr(0, 5);
      ind += 4;
      continue;
    }

    if (name[ind] == '(') {
      fl = true;
      ret += name.substr(ind);
      break;
    }

    for (int j = 0; j < strs_to.size() && !fl; ++j) {
      if (ind + strs_to[j].size() <= name.size()
        && strs_to[j] == name.substr(ind, strs_to[j].size())) {
        //qDebug() << QString(strs_to[j].c_str()) << " " << QString(strs_from[j].c_str());
        fl = true;
        ret += strs_from[j];
        ind += strs_to[j].size() - 1;
        break;
      }
    }
    if (!fl)
      ret += name[ind];
  }
  return ret;
}

DeviationsPlot::DeviationsPlot(QWidget * parent, thermalfist::ThermalModelFit * fit)
{
  this->setObjectName(tr("FitDeviations"));

  int Npoints = fit->FittedQuantities().size();

  if (Npoints != fit->ModelDataSize())
    this->close();

  double tmax = 0.;

  for (int i = 0; i < Npoints; ++i) {
    dataX.push_back(0.5 + i);
    deviationsValues.push_back((fit->ModelData(i) - fit->FittedQuantities()[i].Value()) / fit->FittedQuantities()[i].ValueError());

    if (fit->FittedQuantities()[i].type == 0) {
      names.push_back(BeautifyName(fit->model()->TPS()->GetNameFromPDG(fit->FittedQuantities()[i].mult.fPDGID)).c_str());
    }
    else {
      names.push_back((BeautifyName(fit->model()->TPS()->GetNameFromPDG(fit->FittedQuantities()[i].ratio.fPDGID1))
        + "/" + BeautifyName(fit->model()->TPS()->GetNameFromPDG(fit->FittedQuantities()[i].ratio.fPDGID2))).c_str());
    }

    tmax = std::max(tmax, fabs(deviationsValues[i]));
  }

  tmax = std::min(3., static_cast<int>(tmax) + 1.);

  this->xAxis->setTickLabels(false);
  //this->yAxis->setTickLabels(false);
  this->xAxis2->setVisible(true);
  this->xAxis2->setTickLabels(false);
  this->yAxis2->setVisible(true);
  this->yAxis2->setTickLabels(false);

  //this->xAxis->setPadding((40 * height()) / 200);
  this->xAxis->setPadding(40);

  this->xAxis->setRange(0., 1. * Npoints);
  this->xAxis2->setRange(0., 1. * Npoints);

  //this->yAxis->grid()->setVisible(false);
  //this->yAxis2->grid()->setVisible(false);

  QSharedPointer<QCPAxisTickerText> textTicker(new QCPAxisTickerText);
  for (int i = 1; i < Npoints; ++i)
    textTicker->addTick(i, names[i]);

  this->xAxis->setTicker(textTicker);
  this->xAxis2->setTicker(textTicker);

  QSharedPointer<QCPAxisTickerText> textTickerY(new QCPAxisTickerText);
  for (int i = -3; i <= 3; ++i)
    textTickerY->addTick(i, QString::number(i));

  this->yAxis->setTicker(textTickerY);
  this->yAxis->setRange(-tmax, tmax);
  this->yAxis2->setTicker(textTickerY);
  this->yAxis2->setRange(-tmax, tmax);

  for (int i = 0; i < fit->ModelDataSize(); ++i) {
    this->addGraph();
    this->graph(i)->setName("Model");
    this->graph(i)->setPen(QPen(Qt::black, 2, Qt::SolidLine));
    this->graph(i)->data()->clear();
    this->graph(i)->addData(0.5 + i - 0.25, deviationsValues[i]);
    this->graph(i)->addData(0.5 + i + 0.25, deviationsValues[i]);
    if (i > 0)
      this->graph(i)->removeFromLegend();
  }

  for (int i = 0; i < deviationsValues.size(); ++i) {
    labels.push_back(new QCPItemText(this));
    QCPItemText *textLabel = labels[i];
    textLabel->setLayer("overlay");
    textLabel->setClipToAxisRect(false);
    textLabel->position->setTypeX(QCPItemPosition::ptPlotCoords);
    textLabel->position->setTypeY(QCPItemPosition::ptAxisRectRatio);
    textLabel->position->setCoords(0.5 + i, 1. + 0.1 * 200 / (double)(height()));
    if (!names[i].startsWith("anti-")) {
      textLabel->setText(names[i]);
    }
    else {
      QFont tfont = textLabel->font();
      tfont.setOverline(true);
      textLabel->setFont(tfont);
      textLabel->setText(names[i].right(names[i].size() - 5));
    }
  }

  this->plotLayout()->insertRow(0);
  this->plotLayout()->addElement(0, 0, new QCPTextElement(this, "(Model-Data)/σ", QFont("sans", font().pointSize() + 2, QFont::Bold)));
}

void DeviationsPlot::resizeEvent(QResizeEvent * event)
{
  for (int i = 0; i < labels.size(); ++i) {
    QCPItemText *textLabel = labels[i];
    textLabel->position->setCoords(0.5 + i, 1. + 0.1 * 200 / (double)(height()));
  }
  QCustomPlot::resizeEvent(event);
}
