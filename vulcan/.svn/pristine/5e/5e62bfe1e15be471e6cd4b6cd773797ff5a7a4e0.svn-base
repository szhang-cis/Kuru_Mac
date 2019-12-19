#include "instancemodel.h"

instanceModel::instanceModel(QObject *parent) :
    QAbstractListModel(parent)
{
}

void instanceModel::addInstance(QString path) {
    if (!this->pathList.contains(path)){
        this->beginInsertRows(QModelIndex(), this->pathList.size(),this->pathList.size());
        if (!this->socketList.contains(path)) {
            pathList.append(path);
            socketInstance *in = new socketInstance(this);
            in->connectToServer(path, QIODevice::ReadWrite);
            this->socketList.insert(path, in);
            emit newInstanceAdded(path);
        }
        this->endInsertRows();
    }
}

QVariant instanceModel::data(const QModelIndex &index, int role) const {
    if (index.isValid() and role == Qt::DisplayRole) {
        return QVariant(this->pathList.at(index.row()));
    } else {
        return QVariant();
    }
}

int instanceModel::rowCount(const QModelIndex &parent) const {
    Q_UNUSED(parent);
    return this->pathList.size();
}
