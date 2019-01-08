/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef FITPARAMETERSMODEL_H
#define FITPARAMETERSMODEL_H

#include <QAbstractTableModel>

#include "HRGFit/ThermalModelFit.h"

class FitParametersModel : public QAbstractTableModel
{
  Q_OBJECT
  thermalfist::ThermalModelFitParameters *fParameters;
  int columnNumber;
public:
  FitParametersModel(QObject *parent, thermalfist::ThermalModelFitParameters *params = NULL);
  int rowCount(const QModelIndex &parent = QModelIndex()) const;
  int columnCount(const QModelIndex &parent = QModelIndex()) const;
  QVariant data(const QModelIndex &index, int role = Qt::DisplayRole) const;
  QVariant headerData(int section, Qt::Orientation orientation, int role) const;
  Qt::ItemFlags flags ( const QModelIndex & index ) const {
      if (!index.isValid()) {
          return Qt::ItemIsEnabled;
      }
      if (index.column() > 1) {
          return QAbstractItemModel::flags(index) | Qt::ItemIsEditable;
      }
      if (index.column() == 1) {
        return QAbstractItemModel::flags(index) | Qt::ItemIsUserCheckable;
      }
      return QAbstractItemModel::flags(index);
  }
  bool setData(const QModelIndex & index, const QVariant & value, int role);
  void setParameters(thermalfist::ThermalModelFitParameters *params) { fParameters = params; }
};

#endif // FITPARAMETERSMODEL_H