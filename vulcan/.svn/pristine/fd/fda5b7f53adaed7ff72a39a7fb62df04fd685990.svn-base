#include "vulcangui.h"
#include "ui_vulcangui.h"

VulcanGui::VulcanGui(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::VulcanGui)
{
    ui->setupUi(this);
    iconPixmap = new QPixmap(":/images/test.png");
    systrayIcon = new QIcon(*iconPixmap);
    this->setWindowIcon(*systrayIcon);
    createActions();
    createTrayIcon();
    systray->show();
    this->setWindowFlags(Qt::Dialog | Qt::Tool | Qt::CustomizeWindowHint);
    this->installEventFilter(this);
    this->currentSocks = QSet<QString>();
    this->fileWatcher = new QFileSystemWatcher(this);
    this->sockDir = new QDir(QDir::homePath() + "/.vulcanSock/");
    this->sockDir->setFilter(QDir::System | QDir::NoSymLinks | QDir::NoDotAndDotDot);
    this->sockDir->setSorting(QDir::Name);
    if (!this->sockDir->exists()) {
        this->sockDir->mkpath(this->sockDir->path());
    }
    this->model = new instanceModel();
    ui->casesCombo->setModel(model);
    this->fileWatcher->addPath(sockDir->path());
    connect(systray, SIGNAL(activated(QSystemTrayIcon::ActivationReason)), this, SLOT(iconActivated(QSystemTrayIcon::ActivationReason)));
    connect(this->ui->killCaseButton, SIGNAL(clicked()), this, SLOT(message()));
    this->ui->killCaseButton->setEnabled(true);
    connect(this->fileWatcher, SIGNAL(directoryChanged(QString)), this, SLOT(directoryChanged(QString)));
    connect(model, SIGNAL(newInstanceAdded(QString)), this, SLOT(newInstanceAdded(QString)));
    this->infoDir = this->sockDir->entryInfoList(QDir::System | QDir::NoSymLinks | QDir::NoDotAndDotDot);
    for (int i=0;i<this->infoDir.size();i++) {
        this->currentSocks.insert(this->infoDir.at(i).absoluteFilePath());
        this->model->addInstance(this->infoDir.at(i).absoluteFilePath());
        //sockList.last()->connectToServer(this->infoDir.at(i).absoluteFilePath(), QIODevice::ReadWrite);
        //connect(sockList.last(), SIGNAL(connected()), this, SLOT(newConnectionSlot()));
    }

}

void VulcanGui::newInstanceAdded(QString path) {
    systray->showMessage("Nueva Instancia", path);
}

void VulcanGui::directoryChanged(QString path) {
    Q_UNUSED(path);
    this->infoDir = QDir(sockDir->path()).entryInfoList(QDir::System | QDir::NoDotAndDotDot | QDir::NoSymLinks);
    for (int i=0;i<this->infoDir.size();i++) {
        if (!this->currentSocks.contains(this->infoDir.at(i).absoluteFilePath())){
            qDebug() << this->infoDir.at(i).absoluteFilePath();
            this->currentSocks.insert(this->infoDir.at(i).absoluteFilePath());
            this->model->addInstance(this->infoDir.at(i).absoluteFilePath());
        }
    }

}

void VulcanGui::newConnectionSlot() {
    qDebug() << "Connected";
}

void VulcanGui::message() {
    systray->showMessage("Hola", "Hola que tal este es un mensaje");
    foreach (const QString &value, this->currentSocks) {
        qDebug() << value;
    }
    QStringList l = this->model->socketList[this->ui->casesCombo->currentText()]->baseStatus();
    this->ui->pidLabel->setText(l[0]);
    this->ui->currentStepLabel->setText(l[1]);
    this->ui->currentIntervalLabel->setText(l[2]);
    this->ui->launchDateLabel->setText(l[3]);
    this->ui->caseNameLabel->setText(l[4]);
}

void VulcanGui::createActions() {
    showAction = new QAction(tr("&Show"), this);
    connect(showAction, SIGNAL(triggered()), this, SLOT(show()));
    quitAction = new QAction(tr("&Quit"), this);
    connect(quitAction, SIGNAL(triggered()), this, SLOT(quit()));
}

void VulcanGui::quit() {
    QApplication::quit();
}

void VulcanGui::iconActivated(QSystemTrayIcon::ActivationReason reason) {
    switch (reason) {
        case QSystemTrayIcon::Trigger:
        case QSystemTrayIcon::DoubleClick:
            if (this->isActiveWindow()){
                this->hide();
            } else {
                this->show();
                this->activateWindow();
                this->raise();
                this->setFocus();
            }
        break;
        default:
        break;
    }
}

void VulcanGui::createTrayIcon() {
    systrayMenu = new QMenu(this);
    systrayMenu->addAction(showAction);
    systrayMenu->addSeparator();
    systrayMenu->addAction(quitAction);

    systray = new QSystemTrayIcon(this);
    systray->setContextMenu(systrayMenu);
    systray->setIcon(*systrayIcon);
}

bool VulcanGui::eventFilter(QObject *obj, QEvent *event) {
    if (event->type() == QEvent::ActivationChange) {
        if (!this->isActiveWindow()){
            this->hide();
        }
    }
    return QObject::eventFilter(obj, event);
}

VulcanGui::~VulcanGui()
{
    delete ui;
}
