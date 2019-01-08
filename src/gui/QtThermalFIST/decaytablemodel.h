/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef DECAYTABLEMODEL_H
#define DECAYTABLEMODEL_H

#include <QAbstractTableModel>

#include "HRGBase/ThermalModelBase.h"

class DecayTableModel : public QAbstractTableModel
{
  Q_OBJECT
  thermalfist::ThermalParticle *fParticle;
  thermalfist::ThermalParticleSystem *fTPS;
  int columnNumber;
  double bratioSum;
public:
  DecayTableModel(QObject *parent, thermalfist::ThermalParticle *Particle = NULL, thermalfist::ThermalParticleSystem *TPS = NULL);
  int rowCount(const QModelIndex &parent = QModelIndex()) const;
  int columnCount(const QModelIndex &parent = QModelIndex()) const;
  QVariant data(const QModelIndex &index, int role = Qt::DisplayRole) const;
  QVariant headerData(int section, Qt::Orientation orientation, int role) const;
  /*Qt::ItemFlags flags ( const QModelIndex & index ) const {
      if (!index.isValid()) {
          return Qt::ItemIsEnabled;
      }
      if (index.column()!=1) {
          return QAbstractItemModel::flags(index) | Qt::ItemIsEditable;
      }
      return QAbstractItemModel::flags(index);
  }*/
  bool setData(const QModelIndex & index, const QVariant & value, int role);
  void addDecay() {
    beginInsertRows(QModelIndex(), fParticle->Decays().size(), fParticle->Decays().size());
    thermalfist::ThermalParticle::ParticleDecaysVector decays = fParticle->Decays();
    decays.push_back(thermalfist::ParticleDecayChannel());
    fParticle->SetDecays(decays);
    endInsertRows();
  }
  void removeDecay(int number) {
    if (number >= 0 && number < fParticle->Decays().size()) {
      beginRemoveRows(QModelIndex(), number, number);
      thermalfist::ThermalParticle::ParticleDecaysVector decays = fParticle->Decays();
      decays.erase(decays.begin() + number);
      fParticle->SetDecays(decays);
      endRemoveRows();
    }
  }
  void addColumn() {
    beginInsertColumns(QModelIndex(), columnNumber, columnNumber);
    columnNumber++;
    endInsertColumns();
  }
  void removeColumn() {
    int tcolumnNumber = 3;

    bratioSum = 0.;
    for (unsigned int i = 0; i < fParticle->Decays().size(); i++) {
      if (fParticle->Decays()[i].mDaughters.size() + 2 > tcolumnNumber) tcolumnNumber = 2 + fParticle->Decays()[i].mDaughters.size();
      bratioSum += fParticle->Decays()[i].mBratio;
    }
    if (tcolumnNumber < columnNumber) {
      beginRemoveColumns(QModelIndex(), columnNumber - 1, columnNumber - 1);
      columnNumber--;
      endRemoveColumns();
    }
  }
};

#endif