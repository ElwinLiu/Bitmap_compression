#include "bmp_compression.h"

#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    BMP_compression w;
    w.show();
    return a.exec();
}
