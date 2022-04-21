#include <iostream>
#include "bmp_compression.h"
#include "ui_bmp_compression.h"
#include "qfile.h"
#include "qfiledialog.h"
#include "qdebug.h"
#include "bitmap.h"
#include "qmessagebox.h"
#include "fstream"
using namespace std;

BITMAPFILEHEADER fileHeader;
BITMAPINFOHEADER infoHeader;
PixelInfo* pColorTable;
/*存储原图的像素宽度高度和位图深度*/
unsigned int Height = 0;
unsigned int Width = 0;
unsigned int bitCount = 0;
float process_value = 0;
QProgressBar* bar;
QTextBrowser* browser;

unsigned int length(unsigned int p);
void dp(int n, unsigned int p[], unsigned int s[], unsigned int q[], unsigned int b[]);
void Traceback(int n, int& m, unsigned int s[], unsigned int l[]);
int Output(unsigned int s[], unsigned int l[], unsigned int b[], int n, int& m);

BMP_compression::BMP_compression(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::BMP_compression)
{
    ui->setupUi(this);
    bar = ui->progressBar;
    browser = ui->textBrowser;

    ui->progressBar->setRange(0,50000 - 1);

    ui->progressBar->setValue(0);
    connect(ui->btn_import, SIGNAL(clicked(bool)), this, SLOT(on_btn_import_click()));
}

BMP_compression::~BMP_compression()
{
    delete ui;
}

void BMP_compression::on_btn_import_click()
{
    QString path = QFileDialog::getOpenFileName(
                this, "选择要压缩的bmp文件",
                "/",
                "*.bmp ");
    qDebug() << "path=" << path;

    //QString 转 Char*
    char*  ch;
    QByteArray ba = path.toLatin1();
    ch=ba.data();

    FILE* fp;
    fp = fopen(ch, "rb");//读取同目录下的image.bmp文件。
    if(fp == NULL) return;

    //如果不先读取bifType，根据C语言结构体Sizeof运算规则——整体大于部分之和，从而导致读文件错位
    unsigned short  fileType;
    fread(&fileType,1,sizeof (unsigned short), fp);

    if (fileType == 0x4d42)
    {
        //rewind(fp);
        //获取 位图信息 中的图像宽度 高度 和 位深度
        fread(&fileHeader, sizeof(BITMAPFILEHEADER), 1, fp);

        fread(&infoHeader, sizeof(BITMAPINFOHEADER), 1, fp);

        Width = infoHeader.biWidth;
        Height = infoHeader.biHeight;
        bitCount = infoHeader.biBitCount;

        if(bitCount == 8)
        {
            pColorTable=new PixelInfo[256];
            fread(pColorTable, sizeof(PixelInfo), 256, fp);
        }

        /*获取 bmp 像素数据*/
        unsigned char *pBmpBuf;
        int lineByte = (Width*bitCount / 8 + 3) / 4 * 4;  // 确保每行 像素个数为 4 的倍数
        int n = lineByte*Height;
        pBmpBuf = new unsigned char[n];
        fread(pBmpBuf, n, 1, fp);

        fclose(fp);

        unsigned int* p = new unsigned int[n];  //第i个灰度像素的二进制数
        unsigned int* s = new unsigned int[n+1];  //存储i个像素所需要的最小位数
        unsigned int* b = new unsigned int[n+1];  //第i段像素所占位数
        unsigned int* q = new unsigned int[n+1];  //第i段有几个像素
        int m = 0; //最终段数
        ui->textBrowser->insertPlainText("正在读取图像像素信息...\n");
        for (int i = 0; i < n; i++)
        {
            process_value += 10000/(float)n;
            ui->progressBar->setValue(process_value);
            p[i] = pBmpBuf[i];
        }

        delete[] pBmpBuf; //释放内存

        dp(n, p, s, q, b);
        Output(s, q, b, n, m);//重写q, b 其中q为每段存储的像素数，b为该段像素数的最大位数

        //准备目录、文件名
        path.chop(0);
        path += ".666";
        ba = path.toLatin1();
        ch=ba.data();

        FILE *fout;
        fout = fopen(ch, "wb");

        fwrite(&fileType, sizeof(unsigned short), 1, fout);
        fwrite(&fileHeader, sizeof(BITMAPFILEHEADER), 1, fout);
        fwrite(&infoHeader, sizeof(BITMAPINFOHEADER), 1, fout);

        if(bitCount == 8)
        {
            fwrite(pColorTable, sizeof(PixelInfo), 256, fout);
        }

        int total_bits = 32 + m * 11; //32位是描述有多少段（即m）;m*11是有m段，每段的8位用来描述段的长度，3位用来描述像素位数
        browser->insertPlainText("正在写入文件...\n");
        for (int i = 1; i <= m; i++)
        {
            process_value += 10000 / (float)m;
            bar->setValue(process_value);
            total_bits += q[i] * b[i]; //算上所有像素有多少位
        }

        int total_bytes = total_bits % 8 == 0 ? total_bits / 8 : total_bits / 8 + 1; //向上取整计算字节数

        Writer w(fout, total_bytes);
        w.Write(m, n, q, b, p);
        fclose(fout);

        bar->setValue(50000 - 1);
        browser->insertPlainText("压缩完成！");
        QMessageBox::information(this, tr("提示"), "压缩已完成！");
        bar->setValue(0);
        process_value = 0;

        delete[] p; delete[] s; delete[] q; delete[] b; delete[] pColorTable;
        pColorTable = nullptr;
    }
}

unsigned int length(unsigned int p)
{
    int count = 0;
    while (p != 0)
    {
        p = p >> 1;
        count++;
    }
    if (count == 0) count = 1; //避免产生0位，出现bug
    return count;
}

void dp(int n, unsigned int p[], unsigned int s[], unsigned int q[], unsigned int b[])
{
    int Lmax = 256, header = 11;
    s[0] = 0;
    browser->insertPlainText("正在对像素进行压缩...\n");
    for (int i = 1; i <= n; i++)
    {
        process_value += 30000 / (float)n;
        bar->setValue(process_value);
        b[i] = length(p[i-1]);//计算像素点p需要的存储位数
        unsigned int bmax = b[i];
        s[i] = s[i - 1] + bmax;
        q[i] = 1;
        for (int j = 2; j <= i && j <= Lmax; j++) {
            if (bmax < b[i - j + 1])
                bmax = b[i - j + 1];
            if (s[i] > s[i - j] + j * bmax) {
                s[i] = s[i - j] + j * bmax;
                q[i] = j;
            }
        }
        s[i] += header;
    }
}

int Maxb(unsigned int s[], unsigned int l[], unsigned int b[], int j)
{
    int bmax = 0;
            if (j == 1) {
                bmax = b[1];
                for (int i = 2; i <= s[j]; i++) {
                    if (bmax < b[i])
                        bmax = b[i];
                }
            } else {
                bmax = b[s[j - 1] + 1];
                for (int i = s[j - 1] + 2; i <= s[j]; i++) {
                    if (bmax < b[i]) {
                        bmax = b[i];
                    }
                }
            }
            return bmax;
}

void Traceback(int n, int& m, unsigned int s[], unsigned int l[])
{
    if (n == 0)
        return;
    Traceback(n - l[n], m, s, l);
    s[m++] = n - l[n]; //重新为s[]数组赋值，用来存储分段位置
}

int Output(unsigned int s[], unsigned int l[], unsigned int b[], int n, int& m)
{
    Traceback(n, m, s, l);
    s[m] = n;
    for (int j = 1; j <= m; j++)
    {
        l[j] = l[s[j]];
        //b[j] = b[s[j]];
        b[j] = Maxb(s, l, b, j);
    }
    return m;
}

void BMP_compression::on_btn_unzip_clicked()
{
    QString path = QFileDialog::getOpenFileName(
                this, "选择要解压的bmp文件",
                "/",
                "*.bmp.666 ");
    qDebug() << "path=" << path;

    //QString 转 Char*
    char*  ch;
    QByteArray ba = path.toLatin1();
    ch=ba.data();
    //获取解压后文件名
    int index = path.lastIndexOf('.');//查找最后一次出现'.'在字符串中的索引
    path = path.mid(0, index);
    index = path.lastIndexOf('.');//查找最后一次出现'.'在字符串中的索引
    path.insert(index, "_recovered");

    FILE* fp;
    fp = fopen(ch, "rb");//读取.bmp.666文件。
    if(fp == NULL) return;


    ba = path.toLatin1();
    ch=ba.data();
    FILE* fout;
    fout = fopen(ch, "wb");

    //如果不先读取bifType，根据C语言结构体Sizeof运算规则——整体大于部分之和，从而导致读文件错位
    unsigned short  fileType;
    fread(&fileType,1,sizeof (unsigned short), fp);

    //获取 位图信息 中的图像宽度 高度 和 位深度
    fread(&fileHeader, sizeof(BITMAPFILEHEADER), 1, fp);
    fread(&infoHeader, sizeof(BITMAPINFOHEADER), 1, fp);
    //读取颜色表信息
    if(bitCount == 8)
    {
        pColorTable=new PixelInfo[256];
        fread(pColorTable, sizeof(PixelInfo), 256, fp);
    }

    fwrite(&fileType, sizeof(unsigned short), 1, fout);
    fwrite(&fileHeader, sizeof(BITMAPFILEHEADER), 1, fout);
    fwrite(&infoHeader, sizeof(BITMAPINFOHEADER), 1, fout);
    if(bitCount == 8) fwrite(pColorTable, sizeof(PixelInfo), 256, fout);

    int m = 0; //记录有多少段

    fread(&m, sizeof(int), 1, fp); //读取有多少段
    unsigned int* b = new unsigned int[m + 1];  //第i段像素所占位数
    unsigned int* q = new unsigned int[m + 1];  //第i段有几个像素

    //还原b、q信息
    int rest_bits = 0;//记录剩下的位数
    int section_num = m * 11 % 8 == 0 ? m * 11 / 8 : m * 11 / 8 + 1;//向上取整
    if (m * 11 % 8 != 0) {rest_bits = 8 - (m * 11 - m * 11 / 8 * 8);}
    unsigned char* datas = new unsigned char[section_num];
    fread(datas, sizeof(char), section_num, fp);
    recover(m, b, q, datas);

    //计算像素总字节数,压缩后的像素字节数
    browser->insertPlainText("\n正在读取压缩文件信息... \n");
    int p_sum = 0;
    for (int i = 1; i <= m; i++)
    {
        process_value += 10000 / (float)m;
        bar->setValue(process_value);
        p_sum += q[i];
    }
    int pBit_sum = 0, qByte_sum = 0;
    for (int i = 1; i <= m; i++)
    {
        process_value += 10000 / (float)m;
        bar->setValue(process_value);
        pBit_sum += q[i]*b[i];
    }
    qByte_sum = pBit_sum % 8 == 0 ? pBit_sum / 8 : pBit_sum / 8 + 1;

    unsigned char *pBmpBuf = new unsigned char[p_sum]; //存储最终像素数据
    unsigned char *pBmpDatas = new unsigned char[qByte_sum];//存储未解析的像素数据

    //解析像素信息
    browser->insertPlainText("正在解析压缩文件信息... \n");
    fread(pBmpDatas, sizeof(unsigned char), qByte_sum, fp); //读取数据
    browser->insertPlainText("正在写入还原后文件信息... \n");
    getPixel(pBmpBuf, pBmpDatas, m, b, q, datas[section_num - 1], rest_bits);

    fwrite(pBmpBuf, p_sum, 1, fout); //写入像素信息
    delete[] b; delete[] q; delete[] pBmpBuf; delete[] datas; delete[] pBmpDatas; delete[] pColorTable;
    pColorTable = nullptr;

    fclose(fout);
    fclose(fp);

    bar->setValue(50000 - 1);
    QMessageBox::information(this, tr("提示"), "解压已完成！");
    browser->insertPlainText("解压完成！");


    bar->setValue(0);
    process_value = 0;
    return;
}

void myFunc1(unsigned int& tmp, int need, int give, unsigned char data)
{
    data = data << (8 - give);
    data = data >> (8 - need);
    tmp = data;
    return;
}

int recover(int m, unsigned int* b, unsigned int* q, unsigned char* datas)
{
    int section_num = m * 11 % 8 == 0 ? m * 11 / 8 : m * 11 / 8 + 1;//向上取整
    unsigned int tmp_q = 0, tmp_b = 0; //存放临时q、b的元素
    bool is_q = true;
    int rest_bit_q = 8, rest_bit_b = 3, rest_bit_datas = 8, cnt = 1;
    for (int i = 0; i < section_num; i++)
    {
        process_value += 10000/(float)section_num;
        bar->setValue(process_value);
        rest_bit_datas = 8;
        while (rest_bit_datas != 0)
        {
            if (is_q)   //如果现在准备读q
            {
                if (rest_bit_datas < rest_bit_q)//如果这一轮不够
                {
                    myFunc1(tmp_q, rest_bit_q, rest_bit_datas, datas[i]);
                    rest_bit_q = rest_bit_q - rest_bit_datas; //更新需求的位数
                    rest_bit_datas = 0; //当前8位被读完
                }
                else//如果这一轮够了
                {
                    unsigned char datas_i =  datas[i] << (8 - rest_bit_datas);
                    datas_i = datas_i >> (rest_bit_datas - rest_bit_q);
                    tmp_q = tmp_q + datas_i;
                    rest_bit_datas -= rest_bit_q; //更新剩余位数
                    rest_bit_q = 8; //开启新的一轮
                    is_q = false; //准备读3位的像素位
                    q[cnt] = tmp_q + 1; //完成第cnt个元素
                    tmp_q = 0;
                }
            }
            else    //如果现在准备读b
            {
                if (rest_bit_datas < rest_bit_b)//如果这一轮不够
                {
                    myFunc1(tmp_b, rest_bit_b, rest_bit_datas, datas[i]);
                    rest_bit_b = rest_bit_b - rest_bit_datas; //更新需求的位数
                    rest_bit_datas = 0; //当前8位被读完
                }
                else//如果这一轮够了
                {
                    //unsigned char c = datas[i] << (8 - rest_bit_datas);
                    //c = c >>  (rest_bit_datas - rest_bit_b);
                    //char datas_i =  (datas[i] << (8 - rest_bit_datas) >> (rest_bit_datas - rest_bit_b));
                    unsigned char datas_i = datas[i] << (8 - rest_bit_datas);
                    datas_i = datas_i >> (8 - rest_bit_b);
                    tmp_b = tmp_b + datas_i;
                    rest_bit_datas -= rest_bit_b; //更新剩余位数
                    rest_bit_b = 3; //开启新的一轮
                    is_q = true; //准备读8位的段长
                    b[cnt] = tmp_b + 1; //完成第cnt个元素
                    cnt++;
                    tmp_b = 0;
                }
            }
        }
    }
    return rest_bit_datas; //返回空出来多少位
}

void getPixel(unsigned char *pBmpBuf, unsigned char * pDatas, int m,
              unsigned int* b, unsigned int* q, unsigned char last, unsigned int rest_bits)
{
    int cnt = 1, p_cnt = 0, tmp_need = 8, data_cnt = 0;
    int over_ = 0; //记录一组有多少像素未读取
    unsigned char tmp = 0;

    //处理最后一个字节中剩余的位
    while (rest_bits > 0)
    {
        unsigned int group = q[cnt]*b[cnt];  //计算一段像素的总位数
        if (group <= rest_bits)
        {//如果这一段总位数不超过rest_bits
            for (unsigned int i = 0; i < q[cnt]; i++)
            {
                pBmpBuf[p_cnt++] = get_hex(rest_bits - i * b[cnt], b[cnt], last);
            }
            rest_bits -= group; //rest_bits减少group位
            cnt++; //进入下一组
        }
        else
        {//如果这一组总位数跨越rest_bits
            while (rest_bits >= b[cnt])
            {//把没跨越的部分读完
                pBmpBuf[p_cnt++] = get_hex(rest_bits, b[cnt], last);
                rest_bits -= b[cnt];
                over_++;
            }
            over_ = q[cnt] - over_;
            if (rest_bits != 0)
            {//如果像素位跨越rest_bits
                tmp = get_hex(rest_bits, rest_bits, last);
                tmp = tmp << (b[cnt] - rest_bits);
                tmp_need = b[cnt];
                tmp_need -= rest_bits; //该像素位还需要tmp_need位来填充一个像素值
                rest_bits = 0; //剩余位数置0
            }
        }
    }
    if (over_ == 0)
    {
        cnt++; //切换下一组
        over_ = q[cnt]; //更新未读的像素数量
    }
    //开始对pData进行处理
    rest_bits = 8; //记录待解析的8位还剩多少
    while (cnt <= m)
    {
        process_value += 20000/(float)(m-cnt+1);
        bar->setValue(process_value);
        while (rest_bits >= tmp_need)
        {//剩余位数满足当前位数需要
            tmp += get_hex(rest_bits, tmp_need, pDatas[data_cnt]);
            pBmpBuf[p_cnt++] = tmp;
            if (p_cnt == 262143)
            {
                int jidlas = 1;
            }
            tmp = 0; //将tmp置0
            rest_bits -= tmp_need;
            over_--; //该组未读的像素量-1
            if (over_ == 0)
            {
                if (cnt == m)
                {
                    cnt++;
                    break;
                }
                cnt++; //切换下一组
                over_ = q[cnt]; //更新未读的像素数量
            }
            tmp_need = b[cnt];  //更新需要的像素量
        }
        if (rest_bits != 0)
        {//如果留了一部分没有读完，同时也不够
            tmp = get_hex(rest_bits, rest_bits, pDatas[data_cnt]);
            tmp = tmp << (b[cnt] - rest_bits);
            tmp_need -= rest_bits; //该像素位还需要tmp_need位来填充一个像素值
            rest_bits = 8; //剩余位数置8
        }
        rest_bits = 8; //进入下一个8位
        data_cnt++; //读取下一个字节
    }
    //    int cnt = 0, rest_pix = 0;
    //    for (int i = 0; i < )
}

//获得从第start位开始，共bits位的十进制
unsigned char get_hex(int start, int bits, unsigned char ch)
{
    ch = ch << (8 - start);
    ch = ch >> (8 - bits);
    return ch;
}
