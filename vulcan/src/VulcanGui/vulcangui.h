#ifndef VULCANGUI_H
#define VULCANGUI_H

#include <QDialog>
#include <QSystemTrayIcon>
#include <QDebug>
#include <QMenu>
#include <QAction>
#include <QIcon>
#include <QPixmap>
#include <QLocalSocket>
#include <QFileSystemWatcher>
#include <QDir>
#include <QList>
#include <QFileInfoList>
#include <QSet>
#include <QStringList>
#include <QAbstractListModel>
#include "instancemodel.h"

namespace Ui {
class VulcanGui;
}

class VulcanGui : public QDialog
{
    Q_OBJECT

public:
    explicit VulcanGui(QWidget *parent = 0);
    ~VulcanGui();

private slots:
    void iconActivated(QSystemTrayIcon::ActivationReason reason);
    void quit();
    void message();
    void directoryChanged(QString path);
    void newConnectionSlot();
    void newInstanceAdded(QString path);

private:
    void createTrayIcon();
    void createActions();
    Ui::VulcanGui *ui;
    QSystemTrayIcon *systray;
    QMenu *systrayMenu;
    QAction *quitAction;
    QAction *showAction;
    QIcon *systrayIcon;
    QPixmap *iconPixmap;
    QFileSystemWatcher *fileWatcher;
    QDir *sockDir;
    QList<QLocalSocket*> sockList;
    QFileInfoList infoDir;
    QSet<QString> currentSocks;
    QStringList dirList;
    instanceModel *model;

protected:
    bool eventFilter(QObject *obj, QEvent *event);
};

#endif // VULCANGUI_H
