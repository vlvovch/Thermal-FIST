#ifndef TABLEMODEL_H
#define TABLEMODEL_H

#include <QAbstractTableModel>

#include "HRGBase/ThermalModelBase.h"
#include "HRGFit/ThermalModelFit.h"
#include "particlespectra.h"


class TableModel : public QAbstractTableModel
{
    Q_OBJECT
    //ThermalParticleSystem TPS;
    thermalfist::ThermalModelBase *model;
    std::vector<int> RowToParticle;
    std::vector<int> RowToParticleAll;
    std::vector<int> RowToParticleStable;
    bool showOnlyStable;
public:
    TableModel(QObject *parent);
    int rowCount(const QModelIndex &parent = QModelIndex()) const ;
    int columnCount(const QModelIndex &parent = QModelIndex()) const;
    QVariant data(const QModelIndex &index, int role = Qt::DisplayRole) const;
    QVariant headerData(int section, Qt::Orientation orientation, int role) const;
    const std::vector<int>& GetRowToParticle() const {
        return RowToParticle;
    }
    static const int columnNumber;// = 8;
    void updateParticle(int index) {
        QModelIndex start_ix = createIndex( index, 0 );
        QModelIndex end_ix = createIndex( index, columnNumber-1 );
        emit( dataChanged( start_ix, end_ix ) );
    }
    void reset();
    void setModel(thermalfist::ThermalModelBase *mod);
    void updateAll() {
        QModelIndex start_ix = createIndex( 0, 0 );
        QModelIndex end_ix = createIndex( RowToParticle.size()-1, columnNumber-1 );
        emit( dataChanged( start_ix, end_ix ) );
    }
    void setOnlyStable(bool onlyStable) {
        beginResetModel();
        showOnlyStable = onlyStable;
        if (showOnlyStable) RowToParticle = RowToParticleStable;
        else RowToParticle = RowToParticleAll;
        endResetModel();
    }
};


class DecayTableModel : public QAbstractTableModel
{
    Q_OBJECT
    thermalfist::ThermalParticle *fParticle;
    thermalfist::ThermalParticleSystem *fTPS;
    int columnNumber;
    double bratioSum;
    //std::vector<int> RowToParticle;
public:
    DecayTableModel(QObject *parent, thermalfist::ThermalParticle *Particle = NULL, thermalfist::ThermalParticleSystem *TPS = NULL);
    int rowCount(const QModelIndex &parent = QModelIndex()) const ;
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
        //beginInsertRows(QModelIndex(),fParticle->Decays().size(),fParticle->Decays().size());
        //fParticle->Decays().push_back(ParticleDecay());
        //endInsertRows();
			beginInsertRows(QModelIndex(), fParticle->Decays().size(), fParticle->Decays().size());
			std::vector<thermalfist::ParticleDecay> decays = fParticle->Decays();
			decays.push_back(thermalfist::ParticleDecay());
			fParticle->SetDecays(decays);
			endInsertRows();
    }
    void removeDecay(int number) {
        if (number>=0 && number<fParticle->Decays().size()) {
            beginRemoveRows(QModelIndex(), number, number);
						std::vector<thermalfist::ParticleDecay> decays = fParticle->Decays();
						decays.erase(decays.begin() + number);
						fParticle->SetDecays(decays);
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

        //qDebug() << Particle->fDecays.size();
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

    /*void LoadParticleTable(const QString & filename);
    const std::vector<int>& GetRowToParticle() const {
        return RowToParticle;
    }
    ThermalParticleSystem* GetTPS() {
        return &TPS;
    }*/
};



class QuantitiesModel : public QAbstractTableModel
{
    Q_OBJECT
    //ThermalParticle *fParticle;
    std::vector<thermalfist::FittedQuantity> *quantities;
    thermalfist::ThermalModelBase *model;
    //ThermalParticleSystem *fTPS;
    //int columnNumber;
    //double bratioSum;
    //std::vector<int> RowToParticle;
public:
    QuantitiesModel(QObject *parent, std::vector<thermalfist::FittedQuantity>  *quants = NULL, thermalfist::ThermalModelBase *modelop = NULL);
    int rowCount(const QModelIndex &parent = QModelIndex()) const ;
    int columnCount(const QModelIndex &parent = QModelIndex()) const;
    QVariant data(const QModelIndex &index, int role = Qt::DisplayRole) const;
    QVariant headerData(int section, Qt::Orientation orientation, int role) const;
		Qt::ItemFlags flags(const QModelIndex & index) const;
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
    void addQuantity() {
        beginInsertRows(QModelIndex(),quantities->size(),quantities->size());
        //fParticle->Decays().push_back(ParticleDecay());
        quantities->push_back(thermalfist::FittedQuantity());
        endInsertRows();
    }
    void removeQuantity(int number) {
        if (number>=0 && number<quantities->size()) {
            beginRemoveRows(QModelIndex(), number, number);
            quantities->erase(quantities->begin() + number);
            endRemoveRows();
        }
    }

    void setModel(thermalfist::ThermalModelBase *mod);
    void setQuantities(std::vector<thermalfist::FittedQuantity>  *quants);

    void updateAll() {
        QModelIndex start_ix = createIndex( 0, 0 );
        QModelIndex end_ix = createIndex( quantities->size()-1, 5-1 );
        emit( dataChanged( start_ix, end_ix ) );
    }

    /*void LoadParticleTable(const QString & filename);
    const std::vector<int>& GetRowToParticle() const {
        return RowToParticle;
    }
    ThermalParticleSystem* GetTPS() {
        return &TPS;
    }*/
};


class SpectraModel : public QAbstractTableModel
{
    Q_OBJECT
    //ThermalParticleSystem TPS;
    ParticlesSpectra *spectra;
    std::vector<int> RowToParticle;
public:
    SpectraModel(QObject *parent);
    int rowCount(const QModelIndex &parent = QModelIndex()) const ;
    int columnCount(const QModelIndex &parent = QModelIndex()) const;
    QVariant data(const QModelIndex &index, int role = Qt::DisplayRole) const;
    QVariant headerData(int section, Qt::Orientation orientation, int role) const;
    const std::vector<int>& GetRowToParticle() const {
        return RowToParticle;
    }
    static const int columnNumber;// = 8;
    void updateParticle(int index) {
        QModelIndex start_ix = createIndex( index, 0 );
        QModelIndex end_ix = createIndex( index, columnNumber-1 );
        emit( dataChanged( start_ix, end_ix ) );
    }
    void reset();
    void setSpectra(ParticlesSpectra *spec);
    void updateAll() {
        QModelIndex start_ix = createIndex( 0, 0 );
        QModelIndex end_ix = createIndex( RowToParticle.size()-1, columnNumber-1 );
        emit( dataChanged( start_ix, end_ix ) );
    }
//    void setOnlyStable(bool onlyStable) {
//        beginResetModel();
//        showOnlyStable = onlyStable;
//        if (showOnlyStable) RowToParticle = RowToParticleStable;
//        else RowToParticle = RowToParticleAll;
//        endResetModel();
//    }
};


#endif // TABLEMODEL_H
