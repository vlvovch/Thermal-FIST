/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "decaytablemodel.h"

#include <QDebug>

using namespace thermalfist;


DecayTableModel::DecayTableModel(QObject *parent, ThermalParticle *Particle, ThermalParticleSystem *TPS)
  :QAbstractTableModel(parent), fParticle(Particle), fTPS(TPS)
{
  columnNumber = 3;
  bratioSum = 0.;

  for (unsigned int i = 0; i < Particle->Decays().size(); i++) {
    if (Particle->Decays()[i].mDaughters.size() + 2 > columnNumber) columnNumber = 2 + Particle->Decays()[i].mDaughters.size();
    bratioSum += Particle->Decays()[i].mBratio;
  }
}


int DecayTableModel::rowCount(const QModelIndex & /*parent*/) const
{
  return fParticle->Decays().size();
}

int DecayTableModel::columnCount(const QModelIndex & /*parent*/) const
{
  return columnNumber;
}

QVariant DecayTableModel::data(const QModelIndex &index, int role) const
{
  int row = index.row();
  int col = index.column();
  ParticleDecayChannel decay;
  switch (role) {
  case Qt::DisplayRole:
    if (row >= fParticle->Decays().size()) return QVariant();
    decay = fParticle->Decays()[row];
    if (col == 0) return decay.mBratio * 100.;
    else if (col == 1) return decay.mBratio / bratioSum * 100.;
    else if ((col - 2) < decay.mDaughters.size()) {
      int PDGID = decay.mDaughters[col - 2];
      return QString("%1   %2").arg(QString::number(PDGID), QString(fTPS->GetNameFromPDG(PDGID).c_str()));
      //if (fTPS->PdgToId(PDGID) == -1) return QString("%1   %2").arg(QString::number(PDGID), "???");
      //else return QString("%1   %2").arg(QString::number(PDGID), QString(fTPS->ParticleByPDG(PDGID).Name().c_str()));
    }
    else return QVariant();
    break;
  case Qt::EditRole:
    if (row >= fParticle->Decays().size()) return QVariant();
    decay = fParticle->Decays()[row];
    if (col == 0) return decay.mBratio * 100.;
    else if (col == 1) return decay.mBratio / bratioSum * 100.;
    else if ((col - 2) < decay.mDaughters.size()) {
      return decay.mDaughters[col - 2];
    }
    else return 0;
    break;
  }
  return QVariant();
}

QVariant DecayTableModel::headerData(int section, Qt::Orientation orientation, int role) const
{
  if (role == Qt::DisplayRole)
  {
    if (orientation == Qt::Horizontal) {
      {
        if (section == 0) return tr("Branching ratio (%)");
        if (section == 1) return tr("Normalized (%)");
        if (section >= 2) return QString(tr("Daughter %1")).arg(section - 1);
      }
    }
    else if (orientation == Qt::Vertical) return section + 1;
  }
  return QVariant();
}

bool DecayTableModel::setData(const QModelIndex & index, const QVariant & value, int role)
{
  int row = index.row();
  int col = index.column();
  if (role == Qt::EditRole) {
    if (col == 0) {
      fParticle->Decays()[row].mBratio = value.toDouble() / 100.;
      return true;
    }
    else if (col != 1 && (col - 2) < fParticle->Decays()[row].mDaughters.size()) {
      if (value.toInt() != 0)
        fParticle->Decays()[row].mDaughters[col - 2] = value.toInt();
      else fParticle->Decays()[row].mDaughters.erase(fParticle->Decays()[row].mDaughters.begin() + (int)(col - 2));
      return true;
    }
    else if (col != 1 && value.toInt() != 0) {
      fParticle->Decays()[row].mDaughters.push_back(value.toInt());
      return true;
    }
    else return false;
  }
  else return false;
}
