#ifndef RESULTDIALOG_H
#define RESULTDIALOG_H

#include <QDialog>
#include <QTableView>
#include <QTableWidget>
#include <QPushButton>
#include <QTextEdit>

#include "tablemodel.h"

class ThermalModelBase;

class ResultDialog : public QDialog
{
    Q_OBJECT

    //ThermalParticle *fParticle;
    //ThermalParticleSystem *fTPS;

    ThermalModelBase *model;


    QTextEdit *parameters;
    QTextEdit *results;

    QString GetParameters();
    QString GetResults();
public:
    explicit  ResultDialog(QWidget *parent = 0, ThermalModelBase *mod = NULL);

signals:

public slots:
    void checkFixTableSize();
};

#endif // DECAYSEDITOR_H
