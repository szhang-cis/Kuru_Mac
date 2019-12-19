#ifndef INSTANCEMODEL_H
#define INSTANCEMODEL_H

#include <QAbstractListModel>
#include <QModelIndex>
#include <QHash>
#include "socketinstance.h"
#include <QDebug>
#include <QIODevice>

class instanceModel : public QAbstractListModel
{
    Q_OBJECT
public:
    explicit instanceModel(QObject *parent = 0);
    int rowCount(const QModelIndex &parent = QModelIndex()) const;
    QVariant data(const QModelIndex &index, int role) const;
    void addInstance(QString path);
    QHash<QString, socketInstance*> socketList;

signals:
    void newInstanceAdded(QString path);

public slots:

private:
    QList<QString> pathList;

};

#endif // INSTANCEMODEL_H
