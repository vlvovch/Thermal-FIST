/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "spectramodel.h"

#include <QDebug>

#include "HRGBase/ThermalModelBase.h"

using namespace thermalfist;


const int SpectraModel::columnNumber = 9;

SpectraModel::SpectraModel(QObject *parent)
  :QAbstractTableModel(parent)
{
  RowToParticle.resize(0);
}

int SpectraModel::rowCount(const QModelIndex & /*parent*/) const
{
  if (RowToParticle.size() > 0) return RowToParticle.size() + spectra->fNetParticles.size() + spectra->fNetCharges.size() + spectra->fTotalCharges.size() + spectra->fPositiveCharges.size() + spectra->fNegativeCharges.size();
  return RowToParticle.size();
}

int SpectraModel::columnCount(const QModelIndex & /*parent*/) const
{
  return columnNumber;
}

QVariant SpectraModel::data(const QModelIndex &index, int role) const
{
  int row = index.row();
  int col = index.column();
  switch (role) {
  case Qt::DisplayRole:
    if (col == 0) {
      if (row < spectra->fNames.size()) return QString(spectra->fNames[RowToParticle[row]].c_str());
      if (row < spectra->fNames.size() + spectra->fNetParticles.size())
        return spectra->fNetParticles[row - spectra->fNames.size()].GetName().c_str();
      if (row == spectra->fNames.size() + spectra->fNetParticles.size()) return QString("net-baryon");
      if (row == spectra->fNames.size() + spectra->fNetParticles.size() + 1) return QString("net-charge");
      if (row == spectra->fNames.size() + spectra->fNetParticles.size() + 2) return QString("net-strangeness");
      if (row == spectra->fNames.size() + spectra->fNetParticles.size() + 3) return QString("net-charm");
      if (row == spectra->fNames.size() + spectra->fNetParticles.size() + spectra->fNetCharges.size()) return QString("baryon hadrons");
      if (row == spectra->fNames.size() + spectra->fNetParticles.size() + spectra->fNetCharges.size() + 1) return QString("charge hadrons");
      if (row == spectra->fNames.size() + spectra->fNetParticles.size() + spectra->fNetCharges.size() + 2) return QString("strange hadrons");
      if (row == spectra->fNames.size() + spectra->fNetParticles.size() + spectra->fNetCharges.size() + 3) return QString("charm hadrons"); 
      if (row == spectra->fNames.size() + spectra->fNetParticles.size() + spectra->fNetCharges.size() + spectra->fTotalCharges.size()) return QString("baryon+ hadrons");
      if (row == spectra->fNames.size() + spectra->fNetParticles.size() + spectra->fNetCharges.size() + spectra->fTotalCharges.size() + 1) return QString("charge+ hadrons");
      if (row == spectra->fNames.size() + spectra->fNetParticles.size() + spectra->fNetCharges.size() + spectra->fTotalCharges.size() + 2) return QString("strange+ hadrons");
      if (row == spectra->fNames.size() + spectra->fNetParticles.size() + spectra->fNetCharges.size() + spectra->fTotalCharges.size() + 3) return QString("charm+ hadrons"); 
      if (row == spectra->fNames.size() + spectra->fNetParticles.size() + spectra->fNetCharges.size() + spectra->fTotalCharges.size() + spectra->fPositiveCharges.size()) return QString("baryon- hadrons");
      if (row == spectra->fNames.size() + spectra->fNetParticles.size() + spectra->fNetCharges.size() + spectra->fTotalCharges.size() + spectra->fPositiveCharges.size() + 1) return QString("charge- hadrons");
      if (row == spectra->fNames.size() + spectra->fNetParticles.size() + spectra->fNetCharges.size() + spectra->fTotalCharges.size() + spectra->fPositiveCharges.size() + 2) return QString("strange- hadrons");
      if (row == spectra->fNames.size() + spectra->fNetParticles.size() + spectra->fNetCharges.size() + spectra->fTotalCharges.size() + spectra->fPositiveCharges.size() + 3) return QString("charm- hadrons");
    }
    if (col == 1) {
      if (row < spectra->fParticles.size()) return spectra->fParticles[RowToParticle[row]].GetPDGID();
      else return "--";
    }
    if (col == 2) {
      if (row < spectra->fMasses.size()) return spectra->fMasses[RowToParticle[row]];
      else return "--";
    }
    {
      if (row < spectra->fParticles.size()) {
        if (spectra->fParticles[RowToParticle[row]].GetEvents() == 0) return QVariant();
        if (col == 3) return QString::number(spectra->fParticles[RowToParticle[row]].GetMean(), 'f', 3) + " ± " + QString::number(spectra->fParticles[RowToParticle[row]].GetMeanError(), 'f', 3);//QString("%1 ± %2").arg(spectra->fParticles[RowToParticle[row]].GetMean()).arg(spectra->fParticles[RowToParticle[row]].GetMeanError());
        if (col == 4) return QString::number(spectra->fParticles[RowToParticle[row]].GetVariance(), 'f', 3) + " ± " + QString::number(spectra->fParticles[RowToParticle[row]].GetVarianceError(), 'f', 3);
        if (col == 5) return QString::number(spectra->fParticles[RowToParticle[row]].GetScaledVariance(), 'f', 3) + " ± " + QString::number(spectra->fParticles[RowToParticle[row]].GetScaledVarianceError(), 'f', 3);
        if (col == 6) return QString::number(spectra->fParticles[RowToParticle[row]].GetSkewness(), 'f', 3) + " ± " + QString::number(spectra->fParticles[RowToParticle[row]].GetSkewnessError(), 'f', 3);
        if (col == 7) return QString::number(spectra->fParticles[RowToParticle[row]].GetKurtosis(), 'f', 3) + " ± " + QString::number(spectra->fParticles[RowToParticle[row]].GetKurtosisError(), 'f', 3);
      }
      else if (row < spectra->fNames.size() + spectra->fNetParticles.size()) {
        int tind = row - spectra->fParticles.size();
        if (spectra->fNetParticles[tind].n == 0.) return QVariant();
        if (col == 3) return QString::number(spectra->fNetParticles[tind].GetMean(), 'f', 3) + " ± " + QString::number(spectra->fNetParticles[tind].GetMeanError(), 'f', 3);//QString("%1 ± %2").arg(spectra->fParticles[RowToParticle[row]].GetMean()).arg(spectra->fParticles[RowToParticle[row]].GetMeanError());
        if (col == 4) return QString::number(spectra->fNetParticles[tind].GetVariance(), 'f', 3) + " ± " + QString::number(spectra->fNetParticles[tind].GetVarianceError(), 'f', 3);
        if (col == 5) return QString::number(spectra->fNetParticles[tind].GetScaledVariance(), 'f', 3) + " ± " + QString::number(spectra->fNetParticles[tind].GetScaledVarianceError(), 'f', 3);
        if (col == 6) return QString::number(spectra->fNetParticles[tind].GetSkewness(), 'f', 3) + " ± " + QString::number(spectra->fNetParticles[tind].GetSkewnessError(), 'f', 3);
        if (col == 7) return QString::number(spectra->fNetParticles[tind].GetKurtosis(), 'f', 3) + " ± " + QString::number(spectra->fNetParticles[tind].GetKurtosisError(), 'f', 3);
      }
      else if (row < spectra->fNames.size() + spectra->fNetParticles.size() + spectra->fNetCharges.size()) {
        int tind = row - spectra->fParticles.size() - spectra->fNetParticles.size();
        if (spectra->fNetCharges[tind].n == 0.) return QVariant();
        if (col == 3) return QString::number(spectra->fNetCharges[tind].GetMean(), 'f', 3) + " ± " + QString::number(spectra->fNetCharges[tind].GetMeanError(), 'f', 3);//QString("%1 ± %2").arg(spectra->fParticles[RowToParticle[row]].GetMean()).arg(spectra->fParticles[RowToParticle[row]].GetMeanError());
        if (col == 4) return QString::number(spectra->fNetCharges[tind].GetVariance(), 'f', 3) + " ± " + QString::number(spectra->fNetCharges[tind].GetVarianceError(), 'f', 3);
        if (col == 5) return QString::number(spectra->fNetCharges[tind].GetScaledVariance(), 'f', 3) + " ± " + QString::number(spectra->fNetCharges[tind].GetScaledVarianceError(), 'f', 3);
        if (col == 6) return QString::number(spectra->fNetCharges[tind].GetSkewness(), 'f', 3) + " ± " + QString::number(spectra->fNetCharges[tind].GetSkewnessError(), 'f', 3);
        if (col == 7) return QString::number(spectra->fNetCharges[tind].GetKurtosis(), 'f', 3) + " ± " + QString::number(spectra->fNetCharges[tind].GetKurtosisError(), 'f', 3);
      }
      else if (row < spectra->fNames.size() + spectra->fNetParticles.size() + spectra->fNetCharges.size() + spectra->fTotalCharges.size()) {
        int tind = row - spectra->fParticles.size() - spectra->fNetParticles.size() - spectra->fNetCharges.size();
        if (spectra->fTotalCharges[tind].n == 0.) return QVariant();
        if (col == 3) return QString::number(spectra->fTotalCharges[tind].GetMean(), 'f', 3) + " ± " + QString::number(spectra->fTotalCharges[tind].GetMeanError(), 'f', 3);//QString("%1 ± %2").arg(spectra->fParticles[RowToParticle[row]].GetMean()).arg(spectra->fParticles[RowToParticle[row]].GetMeanError());
        if (col == 4) return QString::number(spectra->fTotalCharges[tind].GetVariance(), 'f', 3) + " ± " + QString::number(spectra->fTotalCharges[tind].GetVarianceError(), 'f', 3);
        if (col == 5) return QString::number(spectra->fTotalCharges[tind].GetScaledVariance(), 'f', 3) + " ± " + QString::number(spectra->fTotalCharges[tind].GetScaledVarianceError(), 'f', 3);
        if (col == 6) return QString::number(spectra->fTotalCharges[tind].GetSkewness(), 'f', 3) + " ± " + QString::number(spectra->fTotalCharges[tind].GetSkewnessError(), 'f', 3);
        if (col == 7) return QString::number(spectra->fTotalCharges[tind].GetKurtosis(), 'f', 3) + " ± " + QString::number(spectra->fTotalCharges[tind].GetKurtosisError(), 'f', 3);
      }
      else if (row < spectra->fNames.size() + spectra->fNetParticles.size() + spectra->fNetCharges.size() + spectra->fTotalCharges.size() + spectra->fPositiveCharges.size()) {
        int tind = row - spectra->fParticles.size() - spectra->fNetParticles.size() - spectra->fNetCharges.size() - spectra->fTotalCharges.size();
        if (spectra->fPositiveCharges[tind].n == 0.) return QVariant();
        if (col == 3) return QString::number(spectra->fPositiveCharges[tind].GetMean(), 'f', 3) + " ± " + QString::number(spectra->fPositiveCharges[tind].GetMeanError(), 'f', 3);//QString("%1 ± %2").arg(spectra->fParticles[RowToParticle[row]].GetMean()).arg(spectra->fParticles[RowToParticle[row]].GetMeanError());
        if (col == 4) return QString::number(spectra->fPositiveCharges[tind].GetVariance(), 'f', 3) + " ± " + QString::number(spectra->fPositiveCharges[tind].GetVarianceError(), 'f', 3);
        if (col == 5) return QString::number(spectra->fPositiveCharges[tind].GetScaledVariance(), 'f', 3) + " ± " + QString::number(spectra->fPositiveCharges[tind].GetScaledVarianceError(), 'f', 3);
        if (col == 6) return QString::number(spectra->fPositiveCharges[tind].GetSkewness(), 'f', 3) + " ± " + QString::number(spectra->fPositiveCharges[tind].GetSkewnessError(), 'f', 3);
        if (col == 7) return QString::number(spectra->fPositiveCharges[tind].GetKurtosis(), 'f', 3) + " ± " + QString::number(spectra->fPositiveCharges[tind].GetKurtosisError(), 'f', 3);
      }
      else if (row < spectra->fNames.size() + spectra->fNetParticles.size() + spectra->fNetCharges.size() + spectra->fTotalCharges.size() + spectra->fPositiveCharges.size() + spectra->fNegativeCharges.size()) {
        int tind = row - spectra->fParticles.size() - spectra->fNetParticles.size() - spectra->fNetCharges.size() - spectra->fTotalCharges.size() - spectra->fPositiveCharges.size();
        if (spectra->fNegativeCharges[tind].n == 0.) return QVariant();
        if (col == 3) return QString::number(spectra->fNegativeCharges[tind].GetMean(), 'f', 3) + " ± " + QString::number(spectra->fNegativeCharges[tind].GetMeanError(), 'f', 3);//QString("%1 ± %2").arg(spectra->fParticles[RowToParticle[row]].GetMean()).arg(spectra->fParticles[RowToParticle[row]].GetMeanError());
        if (col == 4) return QString::number(spectra->fNegativeCharges[tind].GetVariance(), 'f', 3) + " ± " + QString::number(spectra->fNegativeCharges[tind].GetVarianceError(), 'f', 3);
        if (col == 5) return QString::number(spectra->fNegativeCharges[tind].GetScaledVariance(), 'f', 3) + " ± " + QString::number(spectra->fNegativeCharges[tind].GetScaledVarianceError(), 'f', 3);
        if (col == 6) return QString::number(spectra->fNegativeCharges[tind].GetSkewness(), 'f', 3) + " ± " + QString::number(spectra->fNegativeCharges[tind].GetSkewnessError(), 'f', 3);
        if (col == 7) return QString::number(spectra->fNegativeCharges[tind].GetKurtosis(), 'f', 3) + " ± " + QString::number(spectra->fNegativeCharges[tind].GetKurtosisError(), 'f', 3);
      }
    }
    if (col == 8 && row < spectra->fParticles.size()) 
      return QString::number(spectra->fParticles[RowToParticle[row]].GetMeanPt(), 'f', 3) + " ± " + QString::number(spectra->fParticles[RowToParticle[row]].GetMeanPtError(), 'f', 3);
    //if (col == 8 && row < spectra->fParticles.size() && spectra->fParticles[RowToParticle[row]].GetAcceptance()) return "*";
    break;
  case Qt::FontRole:
    //            if (model->TPS()->Particles()[RowToParticle[row]].IsStable()) //change font only for cell(0,0)
    //            {
    //                QFont boldFont;
    //                boldFont.setBold(true);
    //                return boldFont;
    //            }
    break;
    /*case Qt::BackgroundRole:

        if (row == 1 && col == 2)  //change background only for cell(1,2)
        {
            QBrush redBackground(Qt::red);
            return redBackground;
        }
        break;*/
  case Qt::TextAlignmentRole:

    //            if (col>=3 && col<=8) //add a checkbox to cell(1,0)
    //            {
    //                return Qt::AlignCenter;
    //            }
    if (col == 8) return Qt::AlignCenter;
    break;
  case Qt::CheckStateRole:

    //            if (col>=3 && col<=8) //add a checkbox to cell(1,0)
    //            {
    //                //return Qt::Unchecked;
    //            }
    break;
  }
  return QVariant();
}

QVariant SpectraModel::headerData(int section, Qt::Orientation orientation, int role) const
{
  if (role == Qt::DisplayRole)
  {
    if (orientation == Qt::Horizontal) {
      switch (section)
      {
      case 1:
        return tr("PDG ID");
      case 0:
        return tr("Name");
      case 2:
        return tr("m [GeV]");
      case 3:
        return tr("Multiplicity");
      case 4:
        return tr("Variance");
      case 5:
        return tr("Scaled variance");
      case 6:
        return tr("Skewness");
      case 7:
        return tr("Kurtosis");
      case 8:
        return tr("<pT> [GeV]");
      }
    }
    else if (orientation == Qt::Vertical) return section + 1;
  }
  return QVariant();
}


void SpectraModel::reset() {
  beginResetModel();
  RowToParticle.resize(0);
  for (int i = 0; i < spectra->fParticles.size(); ++i) {
    RowToParticle.push_back(i);
  }
  endResetModel();
}

void SpectraModel::setSpectra(ParticlesSpectra *spec) {
  beginResetModel();
  spectra = spec;
  RowToParticle.resize(0);
  if (spectra != NULL) {
    for (int i = 0; i < spectra->fParticles.size(); ++i) {
      RowToParticle.push_back(i);
    }
  }
  endResetModel();
}
