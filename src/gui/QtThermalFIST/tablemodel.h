/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef TABLEMODEL_H
#define TABLEMODEL_H

#include <QAbstractTableModel>

#include "HRGBase/ThermalModelBase.h"
#include "HRGFit/ThermalModelFit.h"
#include "particlespectra.h"


class TableModel : public QAbstractTableModel
{
    Q_OBJECT
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


#endif // TABLEMODEL_H
