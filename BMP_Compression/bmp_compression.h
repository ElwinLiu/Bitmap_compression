#ifndef BMP_COMPRESSION_H
#define BMP_COMPRESSION_H

#include <QMainWindow>
#include "Write_data.h"

QT_BEGIN_NAMESPACE
namespace Ui { class BMP_compression; }
QT_END_NAMESPACE

int recover(int m, unsigned int* b, unsigned int* q, unsigned char* datas);
void myFunc1(unsigned int& tmp, int need, int give, unsigned char data);
void getPixel(unsigned char *pBmpBuf, unsigned char * pDatas, int m,
              unsigned int* b, unsigned int* q, unsigned char last, unsigned int rest_bits);
unsigned char get_hex(int start, int bits, unsigned char ch);

class BMP_compression : public QMainWindow
{
    Q_OBJECT

public:
    BMP_compression(QWidget *parent = nullptr);
    ~BMP_compression();

private slots:
    void on_btn_import_click();

    void on_btn_unzip_clicked();

private:
    Ui::BMP_compression *ui;
};
#endif // BMP_COMPRESSION_H
