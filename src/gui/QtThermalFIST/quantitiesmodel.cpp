/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "quantitiesmodel.h"

#include <QDebug>
#include <QFont>
#include <QDir>
#include <fstream>

using namespace thermalfist;


QuantitiesModel::QuantitiesModel(QObject *parent, std::vector<FittedQuantity>  *quants, ThermalModelFit *fitop)
  : quantities(quants), fitcopy(fitop)
{
}


int QuantitiesModel::rowCount(const QModelIndex & /*parent*/) const
{
  return quantities->size();
}

int QuantitiesModel::columnCount(const QModelIndex & /*parent*/) const
{
  return 8;
}

QVariant QuantitiesModel::data(const QModelIndex &index, int role) const
{
  int row = index.row();
  int col = index.column();
  FittedQuantity quantity;
  switch (role) {
  case Qt::DisplayRole:
    if (row >= quantities->size()) return QVariant();
    quantity = (*quantities)[row];
    if (col == 0) {
      if (quantity.type == FittedQuantity::Multiplicity) {
        if (fitcopy != NULL)
          return QString::fromStdString(fitcopy->model()->TPS()->GetNameFromPDG(quantity.mult.fPDGID));
        else
          return quantity.mult.fPDGID;
      }
      else {
        QString name1, name2;
        if (fitcopy != NULL) {
          name1 = QString::fromStdString(fitcopy->model()->TPS()->GetNameFromPDG(quantity.ratio.fPDGID1));
          name2 = QString::fromStdString(fitcopy->model()->TPS()->GetNameFromPDG(quantity.ratio.fPDGID2));
        }
        else {
          name1 = QString::number(quantity.ratio.fPDGID1);
          name2 = QString::number(quantity.ratio.fPDGID2);
        }

        return name1 + "/" + name2;
      }
    }
    else if (col == 2) {
      if (quantity.type == FittedQuantity::Multiplicity) return quantity.mult.fValue;
      else return quantity.ratio.fValue;
    }
    else if (col == 3) {
      if (quantity.type == FittedQuantity::Multiplicity) return quantity.mult.fError;
      else return quantity.ratio.fError;
    }
    else if (col == 4 && fitcopy != NULL && row < fitcopy->ModelDataSize()) {
      if (quantity.type == FittedQuantity::Multiplicity) {
        return fitcopy->ModelData(row);
      }
      else {
        return fitcopy->ModelData(row);
      }
    }
    else if (col == 5 && fitcopy != NULL && row < fitcopy->ModelDataSize()) {
      if (quantity.type == FittedQuantity::Multiplicity) {

        double val = fitcopy->ModelData(row);

        double valexp = quantity.mult.fValue;
        double err = quantity.mult.fError;

        return (val - valexp) / err;
      }
      else {
        double val = fitcopy->ModelData(row);

        double valexp = quantity.ratio.fValue;
        double err = quantity.ratio.fError;

        return (val - valexp) / err;
      }
    }
    else if (col == 6 && fitcopy != NULL && row < fitcopy->ModelDataSize()) {
      if (quantity.type == FittedQuantity::Multiplicity) {
        double val = fitcopy->ModelData(row);

        double valexp = quantity.mult.fValue;
        double err = quantity.mult.fError;

        QString str1 = QString::number(valexp / val, 'f', 3);
        QString str2 = QString::number(err / val, 'f', 3);

        return str1 + " ± " + str2;
      }
      else {
        double val = fitcopy->ModelData(row);

        double valexp = quantity.ratio.fValue;
        double err = quantity.ratio.fError;

        QString str1 = QString::number(valexp / val, 'f', 3);
        QString str2 = QString::number(err / val, 'f', 3);

        return str1 + " ± " + str2;
      }
    }
    else if (col == 7) {
      if (quantity.type == FittedQuantity::Multiplicity) {
        if (quantity.mult.fFeedDown == 0) return QString(tr("Primordial"));
        if (quantity.mult.fFeedDown == 1) return QString(tr("Stability flags"));
        if (quantity.mult.fFeedDown == 2) return QString(tr("Strong+EM+weak decays"));
        if (quantity.mult.fFeedDown == 3) return QString(tr("Strong+EM decays"));
        if (quantity.mult.fFeedDown == 4) return QString(tr("Strong decays"));
      }
      else {
        QString str1, str2;
        if (quantity.ratio.fFeedDown1 == 0) str1 = QString(tr("Primordial"));
        if (quantity.ratio.fFeedDown1 == 1) str1 = QString(tr("Stability flags"));
        if (quantity.ratio.fFeedDown1 == 2) str1 = QString(tr("Strong+EM+weak decays"));
        if (quantity.ratio.fFeedDown1 == 3) str1 = QString(tr("Strong+EM decays"));
        if (quantity.ratio.fFeedDown1 == 4) str1 = QString(tr("Strong decays"));
        if (quantity.ratio.fFeedDown2 == 0) str2 = QString(tr("Primordial"));
        if (quantity.ratio.fFeedDown2 == 1) str2 = QString(tr("Stability flags"));
        if (quantity.ratio.fFeedDown2 == 2) str2 = QString(tr("Strong+EM+weak decays"));
        if (quantity.ratio.fFeedDown2 == 3) str2 = QString(tr("Strong+EM decays"));
        if (quantity.ratio.fFeedDown2 == 4) str2 = QString(tr("Strong decays"));
        return str1 + "/" + str2;
      }
      return QVariant();
    }
    else return QVariant();
    break;
    /*case Qt::EditRole:
        if (row>=fParticle->Decays().size()) return QVariant();
        decay = fParticle->Decays()[row];
        if (col==0) return decay.mBratio * 100.;
        else if (col==1) return decay.mBratio / bratioSum * 100.;
        else if ((col-2)<decay.mDaughters.size()){
            return decay.mDaughters[col-2];
        }
        else return 0;
        break;
    case Qt::FontRole:
        if (TPS.fParticles[RowToParticle[row]].IsStable()) //change font only for cell(0,0)
        {
            QFont boldFont;
            boldFont.setBold(true);
            return boldFont;
        }
        break;
    case Qt::BackgroundRole:

        if (row == 1 && col == 2)  //change background only for cell(1,2)
        {
            QBrush redBackground(Qt::red);
            return redBackground;
        }
        break;*/
  case Qt::TextAlignmentRole:

    if (col == 1) //add a checkbox to cell(1,0)
    {
      return Qt::AlignCenter;
    }
    break;
  case Qt::CheckStateRole:
    if (row >= quantities->size()) return QVariant();
    quantity = (*quantities)[row];
    if (col == 1) //add a checkbox to cell(1,0)
    {
      if (quantity.toFit)
        return Qt::Checked;
      else
        return Qt::Unchecked;
    }
  }
  return QVariant();
}

QVariant QuantitiesModel::headerData(int section, Qt::Orientation orientation, int role) const
{
  if (role == Qt::DisplayRole)
  {
    if (orientation == Qt::Horizontal) {
      {
        if (section == 0) return tr("Name");
        if (section == 1) return tr("Fit?");
        if (section == 2) return tr("Exp. value");
        if (section == 3) return tr("Exp. error");
        if (section == 4) return tr("Model value");
        if (section == 5) return tr("Deviation");
        if (section == 6) return tr("Data/Model");
        if (section == 7) return tr("Feeddown");
      }
    }
    else if (orientation == Qt::Vertical) return section + 1;
  }
  return QVariant();
}

Qt::ItemFlags QuantitiesModel::flags(const QModelIndex &index) const
{
  Qt::ItemFlags flags = QAbstractItemModel::flags(index);
  if (index.column() == 1)
  {
    flags |= Qt::ItemIsUserCheckable;
    flags |= Qt::ItemIsSelectable;
  }
  return flags;
}

bool QuantitiesModel::setData(const QModelIndex & index, const QVariant & value, int role)
{
  int row = index.row();
  int col = index.column();
  if (row >= quantities->size())
    return false;
  FittedQuantity& quantity = (*quantities)[row];
  if (role == Qt::CheckStateRole) {
    quantity.toFit = value.toBool();
  }
  emit dataChanged(index, index);
  return true;
}

void QuantitiesModel::setModel(ThermalModelFit *fitop) {
  beginResetModel();
  fitcopy = fitop;
  endResetModel();
}

void QuantitiesModel::setQuantities(std::vector<FittedQuantity>  *quants) {
  beginResetModel();
  quantities = quants;
  endResetModel();
}

