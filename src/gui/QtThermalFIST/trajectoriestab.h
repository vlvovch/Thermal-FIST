/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef TRAJECTORIESTAB_H
#define TRAJECTORIESTAB_H

#include <QWidget>
#include <QTabWidget>


class TrajectoriesTab : public QWidget
{
    Q_OBJECT

    QTabWidget *tabWidget;

public:
    TrajectoriesTab(QWidget *parent = 0, const QVector<QWidget*>& tabs = QVector<QWidget*>(),
            const QVector<QString>& tabnames = QVector<QString>());

signals:

public slots:

private:
};


#endif // EQUATIONOFSTATETAB_H
