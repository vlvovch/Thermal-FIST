/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "tablemodel.h"

#include <QDebug>
#include <QFont>
#include <QDir>
#include <fstream>

#include "HRGBase/ThermalModelBase.h"

using namespace thermalfist;

const int TableModel::columnNumber = 18;

TableModel::TableModel(QObject *parent)
    :QAbstractTableModel(parent)
{
    RowToParticle.resize(0);
    RowToParticleAll.resize(0);
    RowToParticleStable.resize(0);
    showOnlyStable = false;
}

int TableModel::rowCount(const QModelIndex & /*parent*/) const
{
   return RowToParticle.size();
}

int TableModel::columnCount(const QModelIndex & /*parent*/) const
{
    return columnNumber;
}

QVariant TableModel::data(const QModelIndex &index, int role) const
{

    int row = index.row();
    int col = index.column();
    switch(role){
        case Qt::DisplayRole:
            if (col==0) return QString(model->TPS()->Particles()[RowToParticle[row]].Name().c_str());
            if (col==1) return model->TPS()->Particles()[RowToParticle[row]].PdgId();
            if (col==2) return model->TPS()->Particles()[RowToParticle[row]].Mass();
            if ((col==3 && model->TPS()->Particles()[RowToParticle[row]].IsStable()) ||
                (col==4 && model->TPS()->Particles()[RowToParticle[row]].IsNeutral()!=0)) return "*";
            if (col == 3 && !model->TPS()->Particles()[RowToParticle[row]].IsStable())
              return QString("%1 decays").arg(model->TPS()->Particles()[RowToParticle[row]].Decays().size());
            if (col==5) {
                if (model->TPS()->Particles()[RowToParticle[row]].BaryonCharge()>0) return "+" + QString::number(model->TPS()->Particles()[RowToParticle[row]].BaryonCharge());
                else if (model->TPS()->Particles()[RowToParticle[row]].BaryonCharge()<0) return QString::number(model->TPS()->Particles()[RowToParticle[row]].BaryonCharge());
            }
            if (col==6 && model->TPS()->Particles()[RowToParticle[row]].ElectricCharge()!=0) {
                if (model->TPS()->Particles()[RowToParticle[row]].ElectricCharge()>0) return "+" + QString::number(model->TPS()->Particles()[RowToParticle[row]].ElectricCharge());
                else return model->TPS()->Particles()[RowToParticle[row]].ElectricCharge();
            }
            if (col==7) {
                if (model->TPS()->Particles()[RowToParticle[row]].Strangeness()>0) return "+" + QString::number(model->TPS()->Particles()[RowToParticle[row]].Strangeness());
                else if (model->TPS()->Particles()[RowToParticle[row]].Strangeness()<0) return QString::number(model->TPS()->Particles()[RowToParticle[row]].Strangeness());
                else if (model->TPS()->Particles()[RowToParticle[row]].AbsoluteStrangeness()>1.e-6) return QString::number(model->TPS()->Particles()[RowToParticle[row]].AbsoluteStrangeness()) + " " + tr("(hidden)");
                else if (model->TPS()->Particles()[RowToParticle[row]].Name().substr(0,2)=="K0") return tr("(hidden)");
            }
            if (col==8) {
                if (model->TPS()->Particles()[RowToParticle[row]].Charm()>0) return "+" + QString::number(model->TPS()->Particles()[RowToParticle[row]].Charm());
                else if (model->TPS()->Particles()[RowToParticle[row]].Charm()<0) return QString::number(model->TPS()->Particles()[RowToParticle[row]].Charm());
                else if (model->TPS()->Particles()[RowToParticle[row]].AbsoluteCharm()>1.e-6) return QString::number(model->TPS()->Particles()[RowToParticle[row]].AbsoluteCharm()) + " " + tr("(hidden)");
            }

            if (model->IsCalculated()) {
                if (col==9) return model->Densities()[RowToParticle[row]];
                if (col==10) return model->Densities()[RowToParticle[row]] * model->Volume();
                if (col==11) return model->TotalDensities()[RowToParticle[row]] * model->Volume();
                if (model->IsFluctuationsCalculated()) {
                  if (col == 12) return model->ScaledVariancePrimordial(RowToParticle[row]);
                  if (col == 13) return model->ScaledVarianceTotal(RowToParticle[row]);
                  if (model->Ensemble() == ThermalModelBase::GCE && model->InteractionModel() == ThermalModelBase::Ideal) {
                    if (col == 14) return model->SkewnessPrimordial(RowToParticle[row]);
                    if (col == 15) return model->SkewnessTotal(RowToParticle[row]);
                    if (col == 16) return model->KurtosisPrimordial(RowToParticle[row]);
                    if (col == 17) return model->KurtosisTotal(RowToParticle[row]);
                  }
                }
            }
            break;
        case Qt::FontRole:
            if (model->TPS()->Particles()[RowToParticle[row]].IsStable()) //change font only for cell(0,0)
            {
                QFont boldFont;
                boldFont.setBold(true);
                return boldFont;
            }
            break;
        /*case Qt::BackgroundRole:

            if (row == 1 && col == 2)  //change background only for cell(1,2)
            {
                QBrush redBackground(Qt::red);
                return redBackground;
            }
            break;*/
        case Qt::TextAlignmentRole:

            if (col>=3 && col<=8) //add a checkbox to cell(1,0)
            {
                return Qt::AlignCenter;
            }
            break;
        case Qt::CheckStateRole:

            if (col>=3 && col<=8) //add a checkbox to cell(1,0)
            {
                //return Qt::Unchecked;
            }
    }
    return QVariant();
}

QVariant TableModel::headerData(int section, Qt::Orientation orientation, int role) const
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
                    return tr("Mass");
                case 3:
                    return tr("Stable?");
                case 5:
                    return tr("B");
                case 4:
                    return tr("Neutral?");
                case 6:
                    return tr("Q");
                case 7:
                    return tr("S");
                case 8:
                    return tr("C");
                case 9:
                    return tr("Prim. density");
                case 10:
                    return tr("Prim. multiplicity");
                case 11:
                    return tr("Total multiplicity");
                case 12:
                    return tr("Scaled variance prim.");
                case 13:
                    return tr("Scaled variance total");
                case 14:
                    return tr("Skewness prim.");
                case 15:
                    return tr("Skewness total");
                case 16:
                    return tr("Kurtosis prim.");
                case 17:
                    return tr("Kurtosis total");
            }
        }
        else if (orientation == Qt::Vertical) return section + 1;
    }
    return QVariant();
}


void TableModel::reset() {
    beginResetModel();
    RowToParticle.resize(0);
    RowToParticleAll.resize(0);
    RowToParticleStable.resize(0);
    for(int i=0;i<model->TPS()->Particles().size();++i) {
        RowToParticleAll.push_back(i);
        if (model->TPS()->Particles()[i].IsStable()) RowToParticleStable.push_back(i);
    }
    if (showOnlyStable) RowToParticle = RowToParticleStable;
    else RowToParticle = RowToParticleAll;
    endResetModel();
}

void TableModel::setModel(ThermalModelBase *mod) {
    beginResetModel();
    model = mod;
    RowToParticle.resize(0);
    RowToParticleAll.resize(0);
    RowToParticleStable.resize(0);
    for(int i=0;i<model->TPS()->Particles().size();++i) {
        RowToParticleAll.push_back(i);
        if (model->TPS()->Particles()[i].IsStable()) RowToParticleStable.push_back(i);
    }
    if (showOnlyStable) RowToParticle = RowToParticleStable;
    else RowToParticle = RowToParticleAll;
    endResetModel();
}

