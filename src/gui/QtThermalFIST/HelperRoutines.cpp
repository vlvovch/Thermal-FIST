#include "HelperRoutines.h"
#include <QKeyEvent>
#include <QApplication>
#include <QClipboard>


void QTableWidgetCC::keyPressEvent(QKeyEvent *event)
{
  if (event->matches(QKeySequence::Copy))
  {
    copyTableViewSelectionToClipBoard(this);
  }
}

void copyTableViewSelectionToClipBoard(QTableView* view)
{
  QAbstractItemModel* model = view->model();
  QItemSelectionModel* selection = view->selectionModel();
  QModelIndexList indexes = selection->selectedIndexes();

  qSort(indexes);

  if (indexes.size() < 1)
    return;

  QString selected_text;
  for (int i = 0; i < indexes.size(); i++)
  {
    QModelIndex current = indexes.at(i);
    QVariant data = model->data(current);
    QString text = data.toString();
    selected_text.append(text);
    if (i == indexes.size() - 1 || current.row() != indexes.at(i + 1).row())
      selected_text.append('\n');
    else
      selected_text.append('\t');
  }

  QApplication::clipboard()->setText(selected_text);
}
