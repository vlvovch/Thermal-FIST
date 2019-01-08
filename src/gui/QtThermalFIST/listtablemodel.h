/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef LISTTABLEMODEL_H
#define LISTTABLEMODEL_H

#include <QAbstractTableModel>
#include "HRGBase/ThermalParticleSystem.h"

class ListTableModel : public QAbstractTableModel
{
    Q_OBJECT
    thermalfist::ThermalParticleSystem *TPS, *TPScopy;
    std::vector<int> RowToParticle;
public:
  ListTableModel(QObject *parent);
  ~ListTableModel();
    int rowCount(const QModelIndex &parent = QModelIndex()) const ;
    int columnCount(const QModelIndex &parent = QModelIndex()) const;
    QVariant data(const QModelIndex &index, int role = Qt::DisplayRole) const;
    QVariant headerData(int section, Qt::Orientation orientation, int role) const;
    //void LoadParticleTable(const QString & filename);
    const std::vector<int>& GetRowToParticle() const {
        return RowToParticle;
    }
    static const int columnNumber;// = 8;
    thermalfist::ThermalParticleSystem* GetTPS() {
        return TPS;
    }
    thermalfist::ThermalParticleSystem* GetTPSorig() {
        return TPScopy;
    }
    void setTPS(thermalfist::ThermalParticleSystem *TPSin) {
      *TPScopy = *TPSin;
      reset();
    }
    void applyChanges() {
      if (TPS != NULL)
        *TPScopy = *TPS;
    }
    bool haveChanges() {
      if (TPS == NULL || TPScopy == NULL)
        return false;
      return !((*TPS) == (*TPScopy));
    }
    void updateParticle(int index) {
        QModelIndex start_ix = createIndex( index, 0 );
        QModelIndex end_ix = createIndex( index, columnNumber-1 );
        emit( dataChanged( start_ix, end_ix ) );
    }
    void reset();
    void addParticle(const thermalfist::ThermalParticle &);
    void removeParticle(int number);
};

class DecayEditorTableModel : public QAbstractTableModel
{
    Q_OBJECT
    thermalfist::ThermalParticle *fParticle;
    thermalfist::ThermalParticleSystem *fTPS;
    int columnNumber;
    double bratioSum;
public:
    DecayEditorTableModel(QObject *parent, thermalfist::ThermalParticle *Particle = NULL, thermalfist::ThermalParticleSystem *TPS = NULL);
    int rowCount(const QModelIndex &parent = QModelIndex()) const ;
    int columnCount(const QModelIndex &parent = QModelIndex()) const;
    QVariant data(const QModelIndex &index, int role = Qt::DisplayRole) const;
    QVariant headerData(int section, Qt::Orientation orientation, int role) const;
    Qt::ItemFlags flags ( const QModelIndex & index ) const {
        if (!index.isValid()) {
            return Qt::ItemIsEnabled;
        }
        if (index.column()!=1) {
            return QAbstractItemModel::flags(index) | Qt::ItemIsEditable;
        }
        return QAbstractItemModel::flags(index);
    }
    bool setData(const QModelIndex & index, const QVariant & value, int role);
    void addDecay() {
        beginInsertRows(QModelIndex(),fParticle->Decays().size(),fParticle->Decays().size());
        fParticle->Decays().push_back(thermalfist::ParticleDecayChannel());
        endInsertRows();
    }
    void removeDecay(int number) {
        if (number>=0 && number<fParticle->Decays().size()) {
            beginRemoveRows(QModelIndex(), number, number);
            fParticle->Decays().erase(fParticle->Decays().begin() + number);
            endRemoveRows();
        }
    }
    void addColumn() {
        beginInsertColumns(QModelIndex(),columnNumber,columnNumber);
        columnNumber++;
        endInsertColumns();
    }
    void removeColumn() {
        int tcolumnNumber = 3;

        bratioSum = 0.;
        for(unsigned int i = 0; i < fParticle->Decays().size(); i++) {
            if (fParticle->Decays()[i].mDaughters.size()+2>tcolumnNumber) tcolumnNumber = 2 + fParticle->Decays()[i].mDaughters.size();
            bratioSum += fParticle->Decays()[i].mBratio;
        }
        if (tcolumnNumber<columnNumber) {
            beginRemoveColumns(QModelIndex(),columnNumber-1,columnNumber-1);
            columnNumber--;
            endRemoveColumns();
        }
    }
};

#endif // LISTTABLEMODEL_H
