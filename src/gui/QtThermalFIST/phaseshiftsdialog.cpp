/*
 * Thermal-FIST package
 *
 * Copyright (c) 2026 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "phaseshiftsdialog.h"

#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QLabel>
#include <QCheckBox>
#include <QLineEdit>
#include <QPushButton>
#include <QTableWidget>
#include <QTableWidgetItem>
#include <QHeaderView>
#include <QAbstractItemView>
#include <QDialogButtonBox>
#include <QFileDialog>
#include <QFileInfo>
#include <QMessageBox>

#include <set>
#include <string>
#include <vector>
#include <map>

#include "ThermalFISTConfig.h"
#include "HRGPhaseShifts/PhaseShiftChannel.h"

using namespace thermalfist;

PhaseShiftsDialog::PhaseShiftsDialog(QWidget* parent,
                                     bool enabled,
                                     const QStringList& configs,
                                     const QSet<QString>& disabledChannels)
  : QDialog(parent), m_configs(configs), m_disabled(disabledChannels), m_populating(false)
{
  setWindowTitle(tr("Phase shifts (S-matrix / Beth-Uhlenbeck)"));

  m_chkEnable = new QCheckBox(tr("Enable phase shifts"));
  m_chkEnable->setToolTip(tr("Add S-matrix / phase-shift channels (from the config files below)\n"
                             "as effective Beth-Uhlenbeck degrees of freedom."));
  m_chkEnable->setChecked(enabled);
  connect(m_chkEnable, SIGNAL(toggled(bool)), this, SLOT(onEnableToggled(bool)));

  QLabel* labelConf = new QLabel(tr("Config files:"));
  m_leConfigs = new QLineEdit();
  m_leConfigs->setReadOnly(true);
  m_btnLoad = new QPushButton(tr("Load configs..."));
  m_btnLoad->setToolTip(tr("Choose one or more phase-shift config (.conf) files."));
  connect(m_btnLoad, SIGNAL(clicked()), this, SLOT(loadConfigs()));

  QHBoxLayout* confLay = new QHBoxLayout();
  confLay->addWidget(labelConf);
  confLay->addWidget(m_leConfigs, 1);
  confLay->addWidget(m_btnLoad);

  m_table = new QTableWidget();
  m_table->setColumnCount(3);
  QStringList headers;
  headers << tr("Channel") << tr("Wave") << tr("Source");
  m_table->setHorizontalHeaderLabels(headers);
  m_table->verticalHeader()->setVisible(false);
  m_table->setSelectionMode(QAbstractItemView::NoSelection);
  m_table->setEditTriggers(QAbstractItemView::NoEditTriggers);
  m_table->horizontalHeader()->setStretchLastSection(true);
  m_table->setMinimumSize(520, 280);
  connect(m_table, SIGNAL(itemChanged(QTableWidgetItem*)), this, SLOT(onItemChanged(QTableWidgetItem*)));

  m_btnAll  = new QPushButton(tr("Select all"));
  m_btnNone = new QPushButton(tr("Select none"));
  connect(m_btnAll,  SIGNAL(clicked()), this, SLOT(selectAll()));
  connect(m_btnNone, SIGNAL(clicked()), this, SLOT(selectNone()));
  QHBoxLayout* selLay = new QHBoxLayout();
  selLay->addWidget(m_btnAll);
  selLay->addWidget(m_btnNone);
  selLay->addStretch(1);

  QDialogButtonBox* buttonBox = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
  connect(buttonBox, SIGNAL(accepted()), this, SLOT(accept()));
  connect(buttonBox, SIGNAL(rejected()), this, SLOT(reject()));

  QVBoxLayout* mainLay = new QVBoxLayout();
  mainLay->addWidget(m_chkEnable);
  mainLay->addLayout(confLay);
  mainLay->addWidget(new QLabel(tr("Channels (toggle individually):")));
  mainLay->addWidget(m_table, 1);
  mainLay->addLayout(selLay);
  mainLay->addWidget(buttonBox);
  setLayout(mainLay);

  updateConfigLabel();
  repopulate();
  onEnableToggled(m_chkEnable->isChecked());
}

bool PhaseShiftsDialog::phaseShiftsEnabled() const
{
  return m_chkEnable->isChecked();
}

QSet<QString> PhaseShiftsDialog::disabledChannels() const
{
  return m_disabled;
}

void PhaseShiftsDialog::updateConfigLabel()
{
  QStringList names;
  for (int i = 0; i < m_configs.size(); ++i)
    names << QFileInfo(m_configs[i]).fileName();
  m_leConfigs->setText(names.join(", "));
  m_leConfigs->setToolTip(m_configs.join("\n"));
}

void PhaseShiftsDialog::repopulate()
{
  m_populating = true;
  m_table->setRowCount(0);

  // One row per PARTIAL WAVE, keyed by "<channel>:<wave>" so each wave can be
  // toggled separately (e.g. pipi_I2 -> an S row and a D row). Rows are deduped by
  // that key across configs (the last config wins for the shown source), in
  // first-appearance order.
  struct Row { QString channel, wave, source, key; };
  std::vector<Row> rows;
  std::map<QString, size_t> rowIndex;
  for (int i = 0; i < m_configs.size(); ++i) {
    std::vector<PhaseShifts::PhaseShiftConfigChannel> chans;
    try {
      chans = PhaseShifts::ListPhaseShiftConfigChannels(m_configs[i].toStdString());
    } catch (const std::exception& e) {
      QMessageBox::warning(this, tr("Phase shifts"),
        tr("Could not read config:\n%1\n\n%2").arg(m_configs[i]).arg(e.what()));
      continue;
    }
    for (size_t c = 0; c < chans.size(); ++c) {
      const PhaseShifts::PhaseShiftConfigChannel& ci = chans[c];
      QString source;
      if (ci.reusesResonance) {
        QStringList codes;
        for (size_t k = 0; k < ci.reusedCodes.size(); ++k)
          codes << QString::number(ci.reusedCodes[k]);
        source = tr("replaces resonance(s): %1").arg(codes.join(", "));
      } else {
        source = tr("synthetic cluster");
      }
      for (size_t w = 0; w < ci.waves.size(); ++w) {
        Row r;
        r.channel = QString::fromStdString(ci.name);
        r.wave    = QString::fromStdString(ci.waves[w]);
        r.source  = source;
        r.key     = r.channel + ":" + r.wave;
        std::map<QString, size_t>::iterator f = rowIndex.find(r.key);
        if (f == rowIndex.end()) { rowIndex[r.key] = rows.size(); rows.push_back(r); }
        else { rows[f->second] = r; }   // later config overrides the shown info
      }
    }
  }

  QSet<QString> present;
  for (size_t r = 0; r < rows.size(); ++r) {
    present.insert(rows[r].key);
    const int row = m_table->rowCount();
    m_table->insertRow(row);
    QTableWidgetItem* itName = new QTableWidgetItem(rows[r].channel);
    itName->setFlags(Qt::ItemIsUserCheckable | Qt::ItemIsEnabled);
    itName->setCheckState(m_disabled.contains(rows[r].key) ? Qt::Unchecked : Qt::Checked);
    itName->setData(Qt::UserRole, rows[r].key);
    m_table->setItem(row, 0, itName);
    m_table->setItem(row, 1, new QTableWidgetItem(rows[r].wave));
    m_table->setItem(row, 2, new QTableWidgetItem(rows[r].source));
  }

  // Drop remembered disables for waves no longer present (keeps the set tidy).
  QSet<QString> stillDisabled;
  for (QSet<QString>::const_iterator it = m_disabled.constBegin(); it != m_disabled.constEnd(); ++it)
    if (present.contains(*it))
      stillDisabled.insert(*it);
  m_disabled = stillDisabled;

  m_table->resizeColumnsToContents();
  m_table->horizontalHeader()->setStretchLastSection(true);
  m_populating = false;
}

void PhaseShiftsDialog::onItemChanged(QTableWidgetItem* item)
{
  if (m_populating || !item || item->column() != 0)
    return;
  const QString key = item->data(Qt::UserRole).toString();   // "<channel>:<wave>"
  if (key.isEmpty())
    return;
  if (item->checkState() == Qt::Checked)
    m_disabled.remove(key);
  else
    m_disabled.insert(key);
}

void PhaseShiftsDialog::loadConfigs()
{
  QString prefix = QString(ThermalFIST_INPUT_FOLDER) + "/list/phaseshifts";
  if (!m_configs.isEmpty())
    prefix = QFileInfo(m_configs.first()).absolutePath();
  QStringList paths = QFileDialog::getOpenFileNames(this, tr("Open phase-shift config file(s)"), prefix,
                   tr("Config files (*.conf *.txt);;All files (*)"));
  if (paths.isEmpty())
    return;
  m_configs = paths;
  updateConfigLabel();
  repopulate();
  if (!m_chkEnable->isChecked())
    m_chkEnable->setChecked(true);   // loading configs implies wanting them on
}

void PhaseShiftsDialog::onEnableToggled(bool on)
{
  m_leConfigs->setEnabled(on);
  m_btnLoad->setEnabled(on);
  m_table->setEnabled(on);
  m_btnAll->setEnabled(on);
  m_btnNone->setEnabled(on);
}

void PhaseShiftsDialog::selectAll()
{
  m_populating = true;
  for (int r = 0; r < m_table->rowCount(); ++r)
    m_table->item(r, 0)->setCheckState(Qt::Checked);
  m_populating = false;
  m_disabled.clear();
}

void PhaseShiftsDialog::selectNone()
{
  m_populating = true;
  m_disabled.clear();
  for (int r = 0; r < m_table->rowCount(); ++r) {
    QTableWidgetItem* it = m_table->item(r, 0);
    it->setCheckState(Qt::Unchecked);
    m_disabled.insert(it->data(Qt::UserRole).toString());
  }
  m_populating = false;
}
