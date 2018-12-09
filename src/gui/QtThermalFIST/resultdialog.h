#ifndef RESULTDIALOG_H
#define RESULTDIALOG_H

#include <QDialog>
#include <QTableView>
#include <QTableWidget>
#include <QPushButton>
#include <QTextEdit>

#include "HRGBase/ThermalModelBase.h"
#include "BaseStructures.h"
#include "tablemodel.h"


class ResultDialog : public QDialog
{
	Q_OBJECT

		//ThermalParticle *fParticle;
		//ThermalParticleSystem *fTPS;

  thermalfist::ThermalModelBase *model;
	ChargesFluctuations *flucts;


	QTextEdit *parameters;
	QTextEdit *results;

	QString GetParameters();
	QString GetResults();
public:
	explicit  ResultDialog(QWidget *parent = 0, thermalfist::ThermalModelBase *mod = NULL, ChargesFluctuations *flucts_in = NULL);

signals:

	public slots :
		void checkFixTableSize();
};

#endif // DECAYSEDITOR_H
