/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "calculationtabledialog.h"

#include <QLayout>
#include <QHeaderView>
#include <QTimer>
#include <QLabel>
#include <QApplication>
#include <QClipboard>
#include <QMenu>
#include <QKeyEvent>
#include <QDebug>
#include <QFileDialog>
#include <QMessageBox>
#include <QBuffer>

#include "HelperRoutines.h"
#include "WasmFileIO.h"

#include "HRGBase/ThermalModelBase.h"

using namespace thermalfist;


void CalculationTable::clear() {
  parameter_name.clear();
  parameter_values.clear();
  temperature_values.clear();
  quantities_names.clear();
  quantities_values.clear();
  densities_names.clear();
  densities_values.clear();
}


CalculationTableDialog::CalculationTableDialog(QWidget* parent, const CalculationTable& table) :
  QDialog(parent), m_table(table)
{
  QVBoxLayout* layout = new QVBoxLayout();

  QHBoxLayout *layLeftTop = new QHBoxLayout();
  layLeftTop->setAlignment(Qt::AlignLeft);

  QLabel *labelSel = new QLabel(tr("Quantity:"));
  comboQuantity = new QComboBox();
  comboQuantity->addItem(tr("Equation of state"));
  comboQuantity->addItem(tr("Hadron densities [fm^-3]"));
  comboQuantity->addItem(tr("Hadron densities [T^3]"));
  comboQuantity->setCurrentIndex(0);
  connect(comboQuantity, SIGNAL(currentIndexChanged(int)), this, SLOT(fillTable()));


  labelFeeddown = new QLabel(tr("Feeddown:"));
  comboFeeddown = new QComboBox();
  comboFeeddown->addItem(tr("Primordial"));
  comboFeeddown->addItem(tr("Stability flags"));
  comboFeeddown->addItem(tr("Strong+EM+weak"));
  comboFeeddown->addItem(tr("Strong+EM"));
  comboFeeddown->addItem(tr("Strong"));
  comboFeeddown->setCurrentIndex(0);
  connect(comboFeeddown, SIGNAL(currentIndexChanged(int)), this, SLOT(fillTable()));

  layLeftTop->addWidget(labelSel);
  layLeftTop->addWidget(comboQuantity);
  layLeftTop->addSpacing(20);
  layLeftTop->addWidget(labelFeeddown);
  layLeftTop->addWidget(comboFeeddown);

  tableValues = new QTableWidget();
  tableValues->installEventFilter(this);
  tableValues->setEditTriggers(QAbstractItemView::NoEditTriggers);
  tableValues->setHorizontalScrollBarPolicy(Qt::ScrollBarAsNeeded);
  configureTableRowHeight(tableValues);

  buttonWriteToFile = new QPushButton(tr("Write to file..."));
  connect(buttonWriteToFile, SIGNAL(clicked()), this, SLOT(writeTableToFile()));

  layout->addLayout(layLeftTop);
  layout->addWidget(tableValues);
  layout->addWidget(buttonWriteToFile, 0, Qt::AlignLeft);


  setLayout(layout);

  this->setWindowTitle(tr("Equation of state properties"));

  fillTable();
  tableValues->setContextMenuPolicy(Qt::CustomContextMenu);
  connect(tableValues, SIGNAL(customContextMenuRequested(QPoint)), this, SLOT(contextMenuRequest(QPoint)));
  tableValues->installEventFilter(this);

  //QTimer::singleShot(0, this, SLOT(checkFixTableSize()));
}

void CalculationTableDialog::checkFixTableSize() {
  int tle = 0;
  for (int i = 0; i < tableValues->horizontalHeader()->count(); ++i)
    tle += tableValues->horizontalHeader()->sectionSize(i);
  tableValues->setMinimumWidth(tle + 23);
  tle = 0;
  for (int i = 0; i < tableValues->horizontalHeader()->count(); ++i)
    tle += tableValues->horizontalHeader()->sectionSize(i);
  tableValues->setMinimumWidth(tle + 23);
}

void CalculationTableDialog::fillTable()
{
  labelFeeddown->setVisible(comboQuantity->currentIndex() != 0);
  comboFeeddown->setVisible(comboQuantity->currentIndex() != 0);


  tableValues->clear();

  tableValues->setColumnCount(m_table.quantities_names.size() + 1);
  if (comboQuantity->currentIndex() != 0) {
    tableValues->setColumnCount(m_table.densities_names.size() + 1);
  }
  tableValues->setRowCount(m_table.parameter_values.size());

  tableValues->setHorizontalHeaderItem(0, new QTableWidgetItem(m_table.parameter_name));
  if (comboQuantity->currentIndex() == 0) {
    for(int ic = 1; ic < m_table.quantities_names.size() + 1; ++ic) {
      tableValues->setHorizontalHeaderItem(ic, new QTableWidgetItem(m_table.quantities_names[ic - 1]));
    }
  } else {
    for(int ic = 1; ic < m_table.densities_names.size() + 1; ++ic) {
      tableValues->setHorizontalHeaderItem(ic, new QTableWidgetItem(m_table.densities_names[ic - 1]));
    }
  }
  tableValues->verticalHeader()->hide();

  for(int ir = 0; ir < m_table.parameter_values.size(); ++ir) {
    tableValues->setItem(ir, 0, new QTableWidgetItem(QString::number(m_table.parameter_values[ir])));
    if (comboQuantity->currentIndex() == 0) {
      for (int ic = 1; ic < m_table.quantities_names.size() + 1; ++ic) {
        tableValues->setItem(ir, ic, new QTableWidgetItem(QString::number(m_table.quantities_values[ic - 1][ir])));
      }
    } else {
      double mult = 1.;
      if (comboQuantity->currentIndex() == 2) {
        mult = 1. / pow(m_table.temperature_values[ir] * xMath::GeVtoifm(), 3);
      }
      int feeddown = comboFeeddown->currentIndex();
      for (int ic = 1; ic < m_table.densities_names.size() + 1; ++ic) {
        tableValues->setItem(ir, ic, new QTableWidgetItem(QString::number(mult * m_table.densities_values[ir][feeddown][ic - 1])));
      }
    }
  }
  tableValues->resizeColumnsToContents();
}

bool CalculationTableDialog::eventFilter(QObject* obj, QEvent* event)
{
  if (obj == tableValues) {
    if (event->type() == QEvent::KeyPress && static_cast<QKeyEvent*>(event)->matches(QKeySequence::Copy)) {
      copyTableViewSelectionToClipBoard(tableValues);
      return true;
    }
    else {
      //return false;
      return QDialog::eventFilter(obj, event);
    }
  }
  else {
    // pass the event on to the parent class
    return QDialog::eventFilter(obj, event);
  }
}

void CalculationTableDialog::contextMenuRequest(QPoint pos)
{
  QMenu* menu = new QMenu(this);
  menu->setAttribute(Qt::WA_DeleteOnClose);

  menu->addAction("Copy table to clipboard", this, SLOT(copyTable()));
  menu->addAction("Write table to file...", this, SLOT(writeTableToFile()));

  menu->popup(tableValues->mapToGlobal(pos));
}

void CalculationTableDialog::copyTable() {
  QString table_text = "";

  table_text.append(m_table.parameter_name);
  if (comboQuantity->currentIndex()==0) {
    for (int ic = 0; ic < m_table.quantities_names.size(); ++ic) {
      table_text.append('\t');
      table_text.append(m_table.quantities_names[ic]);
    }
  } else {
    for (int ic = 0; ic < m_table.densities_names.size(); ++ic) {
      table_text.append('\t');
      table_text.append(m_table.densities_names[ic]);
    }
  }
  table_text.append('\n');

  for (int i = 0; i < tableValues->rowCount(); ++i) {
    for (int j = 0; j < tableValues->columnCount(); ++j) {
      if (tableValues->item(i, j) != NULL)
        table_text.append(tableValues->item(i, j)->text());
      if (j != tableValues->columnCount() - 1)
        table_text.append('\t');
      else
        table_text.append('\n');
    }
  }

  QApplication::clipboard()->setText(table_text);
}

void CalculationTableDialog::writeTableToFile() {
#ifdef Q_OS_WASM
  // WASM: Write to buffer and trigger browser download
  QByteArray data;
  QBuffer buffer(&data);
  buffer.open(QIODevice::WriteOnly | QIODevice::Text);

  const int tabsize = 20;

  QTextStream out(&buffer);

  if (comboQuantity->currentIndex() != 0) {
    out << "# Hadron densities table" << Qt::endl;
    QString units = (comboQuantity->currentIndex() == 1) ? "fm^-3" : "T^3";
    out << "# Each density is in units of " << units << Qt::endl;
    out << "# Feeddown: " << comboFeeddown->currentText() << Qt::endl;
  }

  out.setFieldWidth(tabsize);
  out.setFieldAlignment(QTextStream::AlignLeft);
  out << m_table.parameter_name;
  if (comboQuantity->currentIndex()==0) {
    for(int ic = 0; ic < m_table.quantities_names.size(); ++ic) {
      out << m_table.quantities_names[ic];
    }
  } else {
    for(int ic = 0; ic < m_table.densities_names.size(); ++ic) {
      out << m_table.densities_names[ic];
    }
  }
  out << qSetFieldWidth(0) << Qt::endl << qSetFieldWidth(tabsize);

  int feeddown = comboFeeddown->currentIndex();
  for(int ir = 0; ir < m_table.parameter_values.size(); ++ir) {
    out << m_table.parameter_values[ir];
    if (comboQuantity->currentIndex()==0) {
      for(int ic = 0; ic < m_table.quantities_names.size(); ++ic) {
        out << m_table.quantities_values[ic][ir];
      }
    } else {
      double mult = 1.;
      if (comboQuantity->currentIndex() == 2) {
        mult = 1. / pow(m_table.temperature_values[ir] * xMath::GeVtoifm(), 3);
      }
      for(int ic = 0; ic < m_table.densities_names.size(); ++ic) {
        out << mult * m_table.densities_values[ir][feeddown][ic];
      }
    }
    out << qSetFieldWidth(0) << Qt::endl << qSetFieldWidth(tabsize);
  }

  buffer.close();
  WasmFileIO::saveFile(data, "table.txt");

#else
  // Native: Use standard file dialog
  QString fileName = QFileDialog::getSaveFileName(this, tr("Save table to file"), "", tr("Text files (*.txt)"));
  if (fileName.isEmpty())
    return;

  QFile file(fileName);
  if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
    QMessageBox::critical(this, tr("Error"), tr("Cannot open file %1 for writing").arg(fileName));
    return;
  }

  const int tabsize = 20;

  QTextStream out(&file);

  if (comboQuantity->currentIndex() != 0) {
    out << "# Hadron densities table" << Qt::endl;
    QString units = (comboQuantity->currentIndex() == 1) ? "fm^-3" : "T^3";
    out << "# Each density is in units of " << units << Qt::endl;
    out << "# Feeddown: " << comboFeeddown->currentText() << Qt::endl;
  }

  out.setFieldWidth(tabsize);
  out.setFieldAlignment(QTextStream::AlignLeft);
  out << m_table.parameter_name;
  if (comboQuantity->currentIndex()==0) {
    for(int ic = 0; ic < m_table.quantities_names.size(); ++ic) {
      out << m_table.quantities_names[ic];
    }
  } else {
    for(int ic = 0; ic < m_table.densities_names.size(); ++ic) {
      out << m_table.densities_names[ic];
    }
  }
  out << qSetFieldWidth(0) << Qt::endl << qSetFieldWidth(tabsize);

  int feeddown = comboFeeddown->currentIndex();
  for(int ir = 0; ir < m_table.parameter_values.size(); ++ir) {
    out << m_table.parameter_values[ir];
    if (comboQuantity->currentIndex()==0) {
      for(int ic = 0; ic < m_table.quantities_names.size(); ++ic) {
        out << m_table.quantities_values[ic][ir];
      }
    } else {
      double mult = 1.;
      if (comboQuantity->currentIndex() == 2) {
        mult = 1. / pow(m_table.temperature_values[ir] * xMath::GeVtoifm(), 3);
      }
      for(int ic = 0; ic < m_table.densities_names.size(); ++ic) {
        out << mult * m_table.densities_values[ir][feeddown][ic];
      }
    }
    out << qSetFieldWidth(0) << Qt::endl << qSetFieldWidth(tabsize);
  }

  file.close();
#endif
}
