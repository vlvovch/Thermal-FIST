/*
 * Thermal-FIST package
 *
 * Copyright (c) 2026 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef PHASESHIFTSDIALOG_H
#define PHASESHIFTSDIALOG_H

#include <QDialog>
#include <QStringList>
#include <QSet>
#include <QString>

class QCheckBox;
class QLineEdit;
class QPushButton;
class QTableWidget;
class QTableWidgetItem;

/**
 * \brief Configure the S-matrix / phase-shift (Beth-Uhlenbeck) degrees of freedom.
 *
 * A master "enable" checkbox, a button to load one or more phase-shift config
 * files, and a table of the channels those configs define - each individually
 * toggleable. The owner (MainWindow) reads the result back via the getters and
 * rebuilds the particle list, applying only the enabled channels.
 */
class PhaseShiftsDialog : public QDialog
{
  Q_OBJECT
public:
  PhaseShiftsDialog(QWidget* parent,
                    bool enabled,
                    const QStringList& configs,
                    const QSet<QString>& disabledChannels);

  bool phaseShiftsEnabled() const;
  QStringList configs() const { return m_configs; }
  /// Partial waves switched OFF in the table, as "<channel>:<wave>" skip keys
  /// (consumed by AddPhaseShiftChannelsFromFile's skip set).
  QSet<QString> disabledChannels() const;

private slots:
  void loadConfigs();
  void onEnableToggled(bool on);
  void onItemChanged(QTableWidgetItem* item);
  void selectAll();
  void selectNone();

private:
  void repopulate();
  void updateConfigLabel();

  QCheckBox*    m_chkEnable;
  QLineEdit*    m_leConfigs;
  QPushButton*  m_btnLoad;
  QPushButton*  m_btnAll;
  QPushButton*  m_btnNone;
  QTableWidget* m_table;

  QStringList   m_configs;     ///< loaded config files (applied in order)
  QSet<QString> m_disabled;    ///< channels currently switched off (authoritative)
  bool          m_populating;  ///< guard so programmatic fills don't fire onItemChanged
};

#endif // PHASESHIFTSDIALOG_H
