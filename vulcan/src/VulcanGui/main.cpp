#include "vulcangui.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    VulcanGui w;

    return a.exec();
}
