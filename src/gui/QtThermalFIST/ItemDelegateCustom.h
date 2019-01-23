/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef ITEMDELEGATECUSTOM
#define ITEMDELEGATECUSTOM

/** Custom item delegate to center-align a checkbow in a QTextView. */
/** Based on https://wiki.qt.io/Technical_FAQ#How_can_I_align_the_checkboxes_in_a_view.3F */

#include <QStyledItemDelegate>
#include <QApplication>
#include <QMouseEvent>
#include <QKeyEvent>

class ItemDelegateCustom : public QStyledItemDelegate
{
public:
  ItemDelegateCustom(QObject *parent = 0)
    : QStyledItemDelegate(parent)
  {
  }

  void paint(QPainter * painter, const QStyleOptionViewItem & option, const QModelIndex & index) const
  {
    //QStyleOptionViewItemV4 viewItemOption(option);

    //if (index.column() == 1) {
    //  const int textMargin = QApplication::style()->pixelMetric(QStyle::PM_FocusFrameHMargin) + 1;
    //  QRect newRect = QStyle::alignedRect(option.direction, Qt::AlignCenter,
    //    QSize(option.decorationSize.width() + 5, option.decorationSize.height()),
    //    QRect(option.rect.x() + textMargin, option.rect.y(),
    //      option.rect.width() - (2 * textMargin), option.rect.height()));
    //  viewItemOption.rect = newRect;
    //}
    //QStyledItemDelegate::paint(painter, viewItemOption, index);

    // Based on http://www.qtcentre.org/threads/19157-QTableView-checkbox-center-with-stylesheet?p=181413#post181413
    //QVariant value = index.data();
    //if (index.column() == 1) {
    //  bool boolVal = index.model()->data(index, Qt::CheckStateRole).toBool();

    //  QStyle *style = qApp->style();

    //  QRect checkBoxRect = style->subElementRect(QStyle::SE_CheckBoxIndicator, &option);
    //  int chkWidth = checkBoxRect.width();
    //  int chkHeight = checkBoxRect.height();

    //  int centerX = option.rect.left() + qMax(option.rect.width() / 2 - chkWidth / 2, 0);
    //  int centerY = option.rect.top() + qMax(option.rect.height() / 2 - chkHeight / 2, 0);
    //  QStyleOptionViewItem modifiedOption(option);
    //  modifiedOption.rect.moveTo(centerX, centerY);
    //  modifiedOption.rect.setSize(QSize(chkWidth, chkHeight));
    //  if (boolVal) {
    //    modifiedOption.state |= QStyle::State_On;
    //  }

    //  style->drawPrimitive(QStyle::PE_IndicatorItemViewItemCheck, &modifiedOption, painter);
    //}
    //else {
    //  QStyledItemDelegate::paint(painter, option, index);
    //}

    //return;

    // Based on http://www.qtcentre.org/threads/29668-How-to-centering-the-checkbox-in-qtableview
    if (index.column() == 1)
    {
      bool data = index.model()->data(index, Qt::CheckStateRole).toBool();
      QStyleOptionButton checkboxstyle;
      QRect checkbox_rect = QApplication::style()->subElementRect(QStyle::SE_CheckBoxIndicator, &checkboxstyle);
      checkboxstyle.rect = option.rect;
      checkboxstyle.rect.setLeft(option.rect.x() +
        option.rect.width() / 2 - checkbox_rect.width() / 2);
      if (data)
        checkboxstyle.state = QStyle::State_On | QStyle::State_Enabled;
      else
        checkboxstyle.state = QStyle::State_Off | QStyle::State_Enabled;


      QApplication::style()->drawControl(QStyle::CE_CheckBox, &checkboxstyle, painter);

      //QStyleOptionViewItemV4 viewItemOption(option);
      //viewItemOption.rect.setLeft(option.rect.x() + option.rect.width() / 2 - checkbox_rect.width() / 2);
      //QStyledItemDelegate::paint(painter, viewItemOption, index);
    }
    else
    {
      QStyledItemDelegate::paint(painter, option, index);
    }
  }

  virtual bool editorEvent(QEvent *event, QAbstractItemModel *model, const QStyleOptionViewItem &option,
    const QModelIndex &index)
  {
    Q_ASSERT(event);
    Q_ASSERT(model);

    // make sure that the item is checkable
    Qt::ItemFlags flags = model->flags(index);
    if (!(flags & Qt::ItemIsUserCheckable) || !(flags & Qt::ItemIsEnabled))
      return false;
    // make sure that we have a check state
    QVariant value = index.data(Qt::CheckStateRole);
    if (!value.isValid())
      return false;
    // make sure that we have the right event type
    if (event->type() == QEvent::MouseButtonRelease) {
      const int textMargin = QApplication::style()->pixelMetric(QStyle::PM_FocusFrameHMargin) + 1;
      QRect checkRect = QStyle::alignedRect(option.direction, Qt::AlignCenter,
        option.decorationSize,
        QRect(option.rect.x() + (2 * textMargin), option.rect.y(),
          option.rect.width() - (2 * textMargin),
          option.rect.height()));
      if (!checkRect.contains(static_cast<QMouseEvent*>(event)->pos()))
        return false;
    }
    else if (event->type() == QEvent::KeyPress) {
      if (static_cast<QKeyEvent*>(event)->key() != Qt::Key_Space&& static_cast<QKeyEvent*>(event)->key() != Qt::Key_Select)
        return false;
    }
    else {
      return false;
    }
    Qt::CheckState state = (static_cast<Qt::CheckState>(value.toInt()) == Qt::Checked
      ? Qt::Unchecked : Qt::Checked);
    return model->setData(index, state, Qt::CheckStateRole);
  }
};

#endif