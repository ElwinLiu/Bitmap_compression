#include "Write_data.h"

Writer::Writer(FILE* f, int n)
{
        datas = new char[n];
        datas[0] = 0;
        index = 0;
        pos = 0;
        size = n;
        fout = f;
}

void Writer::store(unsigned char data, int bits){
    int rest_bits = 8 - pos;
    if (rest_bits > bits){
        //不需要写入下一个字节
        data = data << (rest_bits - bits);
        datas[index] = datas[index] | data;
        pos += bits;
    }
    else{
        //需要下一个字节
        datas[index] = datas[index] | (data >> (bits - rest_bits));
        index++;
        datas[index] = data << (8 - bits + rest_bits);
        pos = (bits - rest_bits);
    }
}

void Writer::Write(int m, int n, unsigned int q[], unsigned int b[], unsigned int p[])
{
    fwrite(&m, sizeof(int), 1, fout);
    for (int i = 1; i <= m; i++)
    {
        store((unsigned char)q[i] - 1, 8); //八位段长
        store((unsigned char)b[i] - 1, 3);//三位像素位数
    }
    int pBit = b[1], cnt = 1, pLen = q[1];
    for (int i = 0; i < n; i++)
    {
        if (pLen == 0) {cnt++; pBit = b[cnt]; pLen = q[cnt];}
        store((unsigned char)p[i], pBit);
        pLen--;
    }
    fwrite(datas, sizeof(char), size, fout);
}

