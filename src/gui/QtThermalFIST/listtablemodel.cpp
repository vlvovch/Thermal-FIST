/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "listtablemodel.h"
#include <QDebug>
#include <QFont>
#include <QDir>
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace thermalfist;

const int ListTableModel::columnNumber = 11;

ListTableModel::ListTableModel(QObject *parent)
    :QAbstractTableModel(parent)
{
    RowToParticle.resize(0);
    //TPS = TPScopy = NULL;
    TPS = new thermalfist::ThermalParticleSystem();
    TPScopy = new thermalfist::ThermalParticleSystem();
}

ListTableModel::~ListTableModel()
{
  if (TPS != NULL)
    delete TPS;
  if (TPScopy != NULL)
    delete TPScopy;
}

int ListTableModel::rowCount(const QModelIndex & /*parent*/) const
{
   return RowToParticle.size();
}

int ListTableModel::columnCount(const QModelIndex & /*parent*/) const
{
    return columnNumber;
}

QVariant ListTableModel::data(const QModelIndex &index, int role) const
{
    /*if (role == Qt::DisplayRole)
    {
       return QString("Row%1, Column%2")
                   .arg(index.row() + 1)
                   .arg(index.column() +1);
    }*/
    int row = index.row();
    int col = index.column();
    const ThermalParticle& part = TPS->Particles()[RowToParticle[row]];
    switch(role){
        case Qt::DisplayRole:
            if (col == 0) return QString(part.Name().c_str());
            if (col == 1) return part.PdgId();
            if (col == 2) return part.Mass();
            if ((col == 3 && part.IsStable()) ||
                (col == 4 && part.IsNeutral() != 0)) return "*";//
                /*    ||
                (col==6 && TPS->fParticles[RowToParticle[row]].fC!=0) ||
                (col==7 && TPS->fParticles[RowToParticle[row]].fAbsS>1.e-6) ||
                (col==8 && TPS->fParticles[RowToParticle[row]].fAbsC>1.e-6)) return "*";*/
            if (col == 5 && part.BaryonCharge() != 0) {
              if (part.BaryonCharge()>0) return "+" + QString::number(part.BaryonCharge());
              else return part.BaryonCharge();
            }
            if (col == 6 && part.ElectricCharge() != 0) {
                if (part.ElectricCharge()>0) return "+" + QString::number(part.ElectricCharge());
                else return part.ElectricCharge();
            }
            if (col == 7) {
              if (part.Strangeness()>0) return "+" + QString::number(part.Strangeness());
              else if (part.Strangeness()<0) return QString::number(part.Strangeness());
              else if (part.AbsoluteStrangeness()>1.e-6) return QString::number(part.AbsoluteStrangeness()) + " " + tr("(hidden)");
              else if (part.Name().substr(0, 2) == "K0") return tr("(hidden)");
            }
            if (col == 8) {
              if (part.Charm()>0) return "+" + QString::number(part.Charm());
              else if (part.Charm()<0) return QString::number(part.Charm());
              else if (part.AbsoluteCharm()>1.e-6) return QString::number(part.AbsoluteCharm()) + " " + tr("(hidden)");
            }
            if (col==9) {
                double bratiosum = 0.;
                for(int i=0;i<part.Decays().size();++i)
                    bratiosum += part.Decays()[i].mBratio;
                return bratiosum;
            }
            if (col==10) {
                QString ret = "";

                //if (!TPS->CheckDecayChargesConservation(RowToParticle[row]))
                //    ret += tr("Charges not OK!");

                std::vector<int> check = TPS->CheckDecayChargesConservationVector(RowToParticle[row]);
                if (check[0] != 1)
                  ret += "B";
                if (check[1] != 1)
                  ret += "Q";
                if (check[2] != 1)
                  ret += "S";
                if (check[3] != 1)
                  ret += "C";

                for (int i = 1; i < ret.size(); ++i) {
                  ret.insert(i, ',');
                  i++;
                }

                if (ret != "")
                  ret += tr(" not conserved!");

                double bratiosum = 0.;
                for(int i=0;i<part.Decays().size();++i)
                    bratiosum += part.Decays()[i].mBratio;

                if (fabs(bratiosum - 1.) > 1.e-8) {
                    if (ret != "")
                        ret += " / ";
                    ret += "Not 100%!";
                }

                if (abs(part.BaryonCharge()) <= 1 && (part.PdgId() % 10) != static_cast<int>(part.Degeneracy() + 0.5)) {
                  if (ret != "")
                    ret += " / ";
                  ret += "PDG code/degeneracy mismatch!";
                }

                return ret;
            }
            break;
        case Qt::FontRole:
            if (part.IsStable()) //change font only for cell(0,0)
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

QVariant ListTableModel::headerData(int section, Qt::Orientation orientation, int role) const
{
    if (role == Qt::DisplayRole)
    {
        if (orientation == Qt::Horizontal) {
            switch (section)
            {
                case 1:
                    return QString("PDG ID");
                case 0:
                    return QString("Name");
                case 2:
                    return QString("Mass");
                case 3:
                    return QString("Stable?");
                case 5:
                    return QString("B");
                case 4:
                    return QString("Neutral?");
                case 6:
                    return QString("Q");
                case 7:
                    return QString("S");
                case 8:
                    return QString("C");
                case 9:
                    return QString("Sum of BRs");
                case 10:
                    return QString("Check decays");
            }
        }
        else if (orientation == Qt::Vertical) return section + 1;
    }
    return QVariant();
}

//void ListTableModel::LoadParticleTable(const QString & filename) {
//    beginResetModel();
//    TPS = ThermalParticleSystem(filename.toStdString());
//    TPScopy = TPS;
//    RowToParticle.resize(0);
//    for(int i=0;i<TPS->Particles().size();++i)
//        if (TPS->Particles()[i].PdgId() > 0) RowToParticle.push_back(i);
//    endResetModel();
//}

void ListTableModel::reset() {
    beginResetModel();
    *TPS = *TPScopy;
    RowToParticle.resize(0);
    for(int i=0;i<TPS->Particles().size();++i)
        if (TPS->Particles()[i].PdgId() > 0) RowToParticle.push_back(i);
    endResetModel();
}

void ListTableModel::addParticle(const ThermalParticle & part) {
    beginInsertRows(QModelIndex(), RowToParticle.size(), RowToParticle.size());
        TPS->AddParticle(part);
        RowToParticle.push_back(TPS->Particles().size()-1);
    endInsertRows();
}

void ListTableModel::removeParticle(int number) {
    beginRemoveRows(QModelIndex(), number, number);
        int index = RowToParticle[number];
        /*TPS->PDGtoID.erase(TPS->fParticles[index].fPDGID);
        //TPS->IDtoPDG.erase(TPS->IDtoPDG.begin() + index);
        TPS->fParticles[index].fPDGID = 0;
        //TPS->fParticles.erase(TPS->fParticles.begin() + index);*/
        TPS->Particle(index).SetPdgId(0);
        TPS->FillPdgMap();
        RowToParticle.erase(RowToParticle.begin() + number);
    endRemoveRows();
}

DecayEditorTableModel::DecayEditorTableModel(QObject *parent, ThermalParticle *Particle, ThermalParticleSystem *TPS)
    :QAbstractTableModel(parent), fParticle(Particle), fTPS(TPS)
{
    //beginResetModel();
    columnNumber = 3;
    bratioSum = 0.;

    //qDebug() << Particle->fDecays.size();
    for(unsigned int i = 0; i < Particle->Decays().size(); i++) {
        if (Particle->Decays()[i].mDaughters.size()+2>columnNumber) columnNumber = 2 + Particle->Decays()[i].mDaughters.size();
        bratioSum += Particle->Decays()[i].mBratio;
    }
    //endResetModel();
}


int DecayEditorTableModel::rowCount(const QModelIndex & /*parent*/) const
{
   return fParticle->Decays().size();
}

int DecayEditorTableModel::columnCount(const QModelIndex & /*parent*/) const
{
    return columnNumber;
}

QVariant DecayEditorTableModel::data(const QModelIndex &index, int role) const
{
    /*if (role == Qt::DisplayRole)
    {
       return QString("Row%1, Column%2")
                   .arg(index.row() + 1)
                   .arg(index.column() +1);
    }*/
    int row = index.row();
    int col = index.column();
    ParticleDecayChannel decay;
    switch(role){
        case Qt::DisplayRole:
            if (row>=fParticle->Decays().size()) return QVariant();
            decay = fParticle->Decays()[row];
            if (col==0) return QString::number(decay.mBratio * 100., 'f', 4);
            else if (col==1) return QString::number(decay.mBratio / bratioSum * 100., 'f', 4);
            else if ((col-2)<decay.mDaughters.size()){
                int PDGID = decay.mDaughters[col-2];
                return QString("%1   %2").arg(QString::number(PDGID), QString(fTPS->GetNameFromPDG(PDGID).c_str()));
                //if (fTPS->PdgToId(PDGID)==-1) return QString("%1   %2").arg(QString::number(PDGID), "???");
                //else return QString("%1   %2").arg(QString::number(PDGID), QString(fTPS->Particles()[fTPS->PdgToId(PDGID)].Name().c_str()));
            }
            else return QVariant();
            break;
        case Qt::EditRole:
            if (row>=fParticle->Decays().size()) return QVariant();
            decay = fParticle->Decays()[row];
            if (col==0) return QString::number(decay.mBratio * 100., 'f', 4);
            else if (col==1) return decay.mBratio / bratioSum * 100.;
            else if ((col-2)<decay.mDaughters.size()){
                return decay.mDaughters[col-2];
            }
            else return 0;
            break;
    }
    return QVariant();
}

QVariant DecayEditorTableModel::headerData(int section, Qt::Orientation orientation, int role) const
{
    if (role == Qt::DisplayRole)
    {
        if (orientation == Qt::Horizontal) {
            {
                if (section==0) return "Branching ratio (%)";
                if (section==1) return "Normalized (%)";
                if (section>=2) return QString("Daughter %1").arg(section-1);
            }
        }
        else if (orientation == Qt::Vertical) return section + 1;
    }
    return QVariant();
}

bool DecayEditorTableModel::setData(const QModelIndex & index, const QVariant & value, int role)
{
    int row = index.row();
    int col = index.column();
    if (role==Qt::EditRole) {
        if (col==0) {
            fParticle->Decays()[row].mBratio = value.toDouble() / 100.;

            bratioSum = 0.;
            for(unsigned int i = 0; i < fParticle->Decays().size(); i++) {
                bratioSum += fParticle->Decays()[i].mBratio;
            }

            return true;
        }
        else if (col!=1 && (col-2)<fParticle->Decays()[row].mDaughters.size()) {
            if (value.toInt()!=0)
                fParticle->Decays()[row].mDaughters[col-2] = value.toInt();
            else fParticle->Decays()[row].mDaughters.erase(fParticle->Decays()[row].mDaughters.begin() + (int)(col - 2));
            return true;
        }
        else if (col!=1 && value.toInt()!=0) {
            fParticle->Decays()[row].mDaughters.push_back(value.toInt());
            return true;
        }
        else return false;
    }
    else return false;
}
