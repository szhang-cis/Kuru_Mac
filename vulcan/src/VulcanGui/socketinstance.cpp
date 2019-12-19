#include "socketinstance.h"

socketInstance::socketInstance(QObject *parent) :
    QLocalSocket(parent)
{
    connect(this, SIGNAL(connected()), this, SLOT(connectedD()));
    connect(this, SIGNAL(error(QLocalSocket::LocalSocketError)), this, SLOT(errorD(QLocalSocket::LocalSocketError)));
    connect(this, SIGNAL(readyRead()), this, SLOT(dataArrived()));
    this->status = 0;
    qDebug() << "sock Created";
}

void socketInstance::connectedD() {
    qDebug() << "Connected";
    this->status = 1;
}

void socketInstance::errorD(QLocalSocket::LocalSocketError socketError) {
    Q_UNUSED(socketError);
    qDebug() << "Error" << this->errorString();
}

void socketInstance::connectToServer(const QString &name, OpenMode openMode) {
    QLocalSocket::connectToServer(name, openMode);
}

void socketInstance::dataArrived() {
    char myData[100] = "";
    this->read(myData, 1);
    int datalen = 0;
    qDebug() << "Command" << myData[0];
    if (myData[0] == 'N') {
        while (datalen < 87) {
            int recvd = this->read(myData+datalen*sizeof(myData), 100);
            if (recvd == -1){
                qDebug() << "Check this";
            }
            datalen = datalen + recvd;
            qDebug() << myData << " - " << datalen;
        }
        QString qtData = myData;
        this->pid = qtData.mid(1,6).toInt();
        this->max_iter = qtData.mid(8,4).toInt();
        this->max_interval = qtData.mid(13,3).toInt();
        this->start_date = qtData.mid(17,19);
        this->case_name = qtData.mid(37,49);
        qDebug() << "PID:" << this->pid;
        qDebug() << "MaxIter:" << this->max_iter;
        qDebug() << "MaxInter:" << this->max_interval;
        qDebug() << "StartDate:" << this->start_date;
        qDebug() << "CaseName:" << this->case_name;
    } else if (myData[0] == 'U') {
        while (datalen < 19) {
            int recvd = this->read(myData+datalen*sizeof(myData), 100);
            if (recvd == -1){
                qDebug() << "Check this";
            }
            datalen = datalen + recvd;
            qDebug() << myData << " - " << datalen;
        }
        qDebug() << myData;
    }
    qDebug() << "func trig" << myData;
}

QStringList socketInstance::baseStatus() {
    QStringList ret;
    ret << QString().setNum(this->pid) << QString().setNum(this->max_iter) << QString().setNum(this->max_interval) << this->start_date << this->case_name;
    return ret;
}
