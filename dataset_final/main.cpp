#include <iostream>
#include <iomanip>
#include <math.h>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <stdlib.h>


#include <string>
#include <cstdlib>
#include <sstream>
#include <cstring> 

#define SIZE COL*ROW
#define COL 1000
#define INI 61
#define ROW 3
#define randint(a, b) (rand() % (b-a))+ a + 1

using namespace std;

void randpoint(double p[COL][ROW])
{
    double a[SIZE];
    srand((int)time(0));
    for (int i = 0; i < SIZE; i++)
    {
        a[i] = rand() / (double)(RAND_MAX / 200) - 100; // -100~100
        a[i] = floor(a[i] * 100) / 100;
    }
}

void randedge(double p[COL][ROW],int index)
{
    double b[30000][3];
    int edgenum = 0;
    int k = 0;
    for (int i = 0; i < COL - 1; i++) // 0~29
    {
        int v[COL - 1];
        int count = 0;
        for (int j = i; j < randint(i, COL - 1); j++)
        {
            int temp = randint(i, COL - 1);
            for (int l = 0; l < count; l++)
            {
                if (temp == v[l])
                {
                    temp = randint(i, COL - 1);
                    l = 0;
                }
            }
            if (temp != v[0])
            {
                v[count++] = temp;
                b[k][0] = i;
                b[k][1] = temp;
                double err = rand() / (double)(RAND_MAX / 20) - 10;
                b[k][2] = sqrt(pow((p[i][0] - p[temp][0]), 2) + pow((p[i][1] - p[temp][1]), 2) + pow((p[i][2] - p[temp][2]), 2)) + err;
                if (b[k][2] < 0) {
                    b[k][2] = sqrt(pow((p[i][0] - p[temp][0]), 2) + pow((p[i][1] - p[temp][1]), 2) + pow((p[i][2] - p[temp][2]), 2));
                }
                b[k][2] = floor(b[k][2] * 100) / 100;
                k++;
                edgenum++;
            }
        }
    }
    cout << edgenum << endl;

    ofstream Edge;





    

    string str;
    ostringstream oss;
    oss << "P" << index << "-" << COL << ".txt";
    str = oss.str();


    Edge.open(str, std::ios::out | std::ios::trunc); // 只写；打开文件，清空内容
    if (!Edge.is_open())
    {
        return;
    }
    for (int i = 0; i < edgenum; i++)
    {
        Edge << b[i][0] << " " << b[i][1] << " " << b[i][2] << "\n";
    }
    Edge.close();
}

int main(int argc, char* argv[])
{

    for (int i = INI; i <= (INI + 9); i++) {
        double k[COL][ROW];
        randpoint(k);
        randedge(k,i);
        }
        


}