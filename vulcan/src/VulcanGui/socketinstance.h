#ifndef SOCKETINSTANCE_H
#define SOCKETINSTANCE_H

#include <QLocalSocket>
#include <QString>
#include <QStringList>

class socketInstance : public QLocalSocket
{
    Q_OBJECT
public:
    explicit socketInstance(QObject *parent = 0);

signals:

public slots:
    void connectedD();
    void errorD(QLocalSocket::LocalSocketError socketError);
    void dataArrived();

public:
    void connectToServer(const QString &name, OpenMode openMode);
    QStringList baseStatus();

private:
    int status;
    int pid, max_iter, max_interval, cur_iter, cur_interval;
    QString case_name, start_date;

};

#endif // SOCKETINSTANCE_H
