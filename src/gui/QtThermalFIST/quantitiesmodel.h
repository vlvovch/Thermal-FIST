/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef QUANTITIESMODEL_H
#define QUANTITIESMODEL_H

#include <QAbstractTableModel>

#include "HRGBase/ThermalModelBase.h"
#include "HRGFit/ThermalModelFit.h"



class QuantitiesModel : public QAbstractTableModel
{
  Q_OBJECT
  std::vector<thermalfist::FittedQuantity> *quantities;
  thermalfist::ThermalModelFit *fitcopy;
  thermalfist::ThermalParticleSystem *TPS;
public:
  QuantitiesModel(QObject *parent, std::vector<thermalfist::FittedQuantity>  *quants = NULL, thermalfist::ThermalModelFit *fitop = NULL);
  int rowCount(const QModelIndex &parent = QModelIndex()) const;
  int columnCount(const QModelIndex &parent = QModelIndex()) const;
  QVariant data(const QModelIndex &index, int role = Qt::DisplayRole) const;
  QVariant headerData(int section, Qt::Orientation orientation, int role) const;
  Qt::ItemFlags flags(const QModelIndex & index) const;
  bool setData(const QModelIndex & index, const QVariant & value, int role);
  void addQuantity() {
    beginInsertRows(QModelIndex(), quantities->size(), quantities->size());
    quantities->push_back(thermalfist::FittedQuantity());
    endInsertRows();
  }
  void removeQuantity(int number) {
    if (number >= 0 && number < quantities->size()) {
      beginRemoveRows(QModelIndex(), number, number);
      quantities->erase(quantities->begin() + number);
      endRemoveRows();
    }
  }

  void setModel(thermalfist::ThermalModelFit *fitop);
  void setQuantities(std::vector<thermalfist::FittedQuantity>  *quants);

  void updateAll() {
    QModelIndex start_ix = createIndex(0, 0);
    QModelIndex end_ix = createIndex(quantities->size() - 1, 5 - 1);
    emit(dataChanged(start_ix, end_ix));
  }

};


#endif
