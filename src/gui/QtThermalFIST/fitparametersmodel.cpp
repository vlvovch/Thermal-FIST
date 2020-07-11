/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "fitparametersmodel.h"

#include <QDebug>

using namespace thermalfist;


FitParametersModel::FitParametersModel(QObject *parent, thermalfist::ThermalModelFitParameters *params)
  :QAbstractTableModel(parent), fParameters(params)
{
  columnNumber = 5;
}


int FitParametersModel::rowCount(const QModelIndex & /*parent*/) const
{
  return ThermalModelFitParameters::ParameterCount;
}

int FitParametersModel::columnCount(const QModelIndex & /*parent*/) const
{
  return columnNumber;
}

QVariant FitParametersModel::data(const QModelIndex &index, int role) const
{
  int row = index.row();
  int col = index.column();
  switch (role) {
  case Qt::DisplayRole:
  case Qt::EditRole:
  {
    if (row >= ThermalModelFitParameters::ParameterCount) return QVariant();

    FitParameter param = fParameters->GetParameter(row);
    double mn = 1.;

    QString name = QString::fromStdString(param.name);
    if (name == "T" || name == "muB" || name == "muQ" || name == "muS"
      || name == "muC" || name == "Tkin") {
      name += " (MeV)";
      mn = 1.e3;
      if (name.startsWith("mu"))
        name = "μ" + name.right(name.size() - 2);
    }

    if (name == "R" || name == "Rc")
      name += " (fm)";

    if (name.startsWith("gamma")) {
      name = "γ" + name.right(name.size() - 5);
    }

    if (col == 0) {
      return name;
    }
    else if (col == 2) {
      return mn * param.value;
    }
    else if (col == 3) {
      return mn * param.xmin;
    }
    else if (col == 4) {
      return mn * param.xmax;
    }
    else return QVariant();
    break;
  }
  case Qt::TextAlignmentRole:

    if (col == 1) //add a checkbox to cell(1,0)
    {
      return Qt::AlignCenter;
    }
    break;
  case Qt::CheckStateRole:
    if (row >= ThermalModelFitParameters::ParameterCount) return QVariant();
    if (col == 1) //add a checkbox to cell(1,0)
    {
      if (fParameters->GetParameter(row).toFit)
        return Qt::Checked;
      else
        return Qt::Unchecked;
    }
  }
  return QVariant();
}

QVariant FitParametersModel::headerData(int section, Qt::Orientation orientation, int role) const
{
  if (role == Qt::DisplayRole)
  {
    if (orientation == Qt::Horizontal) {
      {
        if (section == 0) return tr("Parameter");
        if (section == 1) return tr("Fit?");
        if (section == 2) return tr("Initial value");
        if (section == 3) return tr("Min value");
        if (section == 4) return tr("Max value");
      }
    }
    else if (orientation == Qt::Vertical) return section + 1;
  }
  return QVariant();
}

bool FitParametersModel::setData(const QModelIndex & index, const QVariant & value, int role)
{
  int row = index.row();
  int col = index.column();
  if (role == Qt::EditRole) {
    FitParameter &param = fParameters->GetParameter(row);
    double mn = 1.;

    std::string name = param.name;
    if (name == "T" || name == "muB" || name == "muQ" || name == "muS"
      || name == "muC" || name == "Tkin") {
      mn = 1.e3;
    }

    if (col == 2) {
      param.value = value.toDouble() / mn;
      return true;
    }
    else if (col == 3) {
      param.xmin = value.toDouble() / mn;
      return true;
    }
    else if (col == 4) {
      param.xmax = value.toDouble() / mn;
      return true;
    }
    else return false;
  }
  else if (role == Qt::CheckStateRole) {
    if (col == 1) {
      fParameters->GetParameter(row).toFit = value.toBool();
      return true;
    }
    else return false;
  }
  else return true;
}
