#include <qfile.h>
#include <fstream>
using namespace std;

class Writer
{
private:
    FILE* fout;
    int size;//写入的字节总数
    char* datas;//数据
    int index;//当前写入到的字节位置
    int pos;//当前写入到的字节的位位置
public:
    Writer(FILE* f, int n);
    void store(unsigned char data, int bits);
    void Write(int m, int n, unsigned int q[], unsigned int b[], unsigned int p[]);
};

