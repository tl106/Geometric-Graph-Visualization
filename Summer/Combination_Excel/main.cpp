#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cmath>
#include <string.h>
#include <iomanip>

#define PI 3.1415926535
#define INFINITE 65535
#define VERTEX_NUM 1000               // node个数
#define AOTO_READ_FILE_START 61     // 文件起始编号
#define EDGE_LENGTH 15              // 初始化球的半径

using namespace std;

// read file
void read_file(float** distance, const char* file_name) {
    // initialize distance
    for (int i = 0; i < VERTEX_NUM; i++) {
        for (int j = 0; j < VERTEX_NUM; j++) {
            distance[i][j] = INFINITE;
        }
    }

    ifstream read;
    read.open(file_name, ios::in);
    if (!read.is_open()) {
        cout << "Input file open failure." << endl;
        exit(EXIT_FAILURE);
    }

    int x, y;
    float d;
    while (!read.eof()) {
        read >> x >> y >> d;
        distance[y][x] = d;
        distance[x][y] = d;
    }

    cout << "Input file read successfully." << endl;
    read.close();
}

//////////////////////////////////////////
//////////////////////////////////////////
// 三维球面上的Marsaglia方法
class cRandom
{
public:
    int seed;
    double random;

    cRandom(int x, double y) :seed(x), random(y) {};
    cRandom() :seed(0), random(0) {};
};

cRandom my_random(int z)
// 16807 way to create random numbers
// z is the seed number, num is the total random number to create
{
    //z(n+1)=(a*z(n)+b) mod m
    //describe m=a*q+r to avoid that the number is large than the computer can bear
    const int m = (int)(pow(2, 31) - 1);
    const int a = 16807;
    const int q = 127773;
    const int r = 2836;

    int temp = a * (z % q) - r * (z / q);

    if (temp < 0)
    {
        temp = m + temp;
    }
    //z is the seed number
    z = temp;
    double t = z * 1.0 / m;

    cRandom cr;
    cr.random = t;
    cr.seed = z;

    return cr;
}

// 初始化每个顶点的坐标，使顶点在一个直径为L_0的球面上较为均匀的分布
void initialize_positions(float** position) {
    int z1 = 20;    // 种子
    int z2 = 112;   // 种子
    cRandom sita(z1, 0.0);
    cRandom pesi(z2, 0.0);

    int count = 0;
    while (count < VERTEX_NUM) {
        sita = my_random(pesi.seed);
        pesi = my_random(sita.seed);

        double u = 2.0 * sita.random - 1.0; // -1 ~ 1
        double v = 2.0 * pesi.random - 1.0;

        double r2 = pow(u, 2) + pow(v, 2);
        if (r2 < 1) {
            // 点(x,y,z)在半径为1的球面上
            double x = 2 * u * sqrt(1 - r2);
            double y = 2 * v * sqrt(1 - r2);
            double z = 1 - 2 * r2;

            position[count][0] = EDGE_LENGTH * x; // EDGE_LENGTH: 半径
            position[count][1] = EDGE_LENGTH * y;
            position[count][2] = EDGE_LENGTH * z;

            count++;
        }
    }
}
//////////////////////////////////////////
//////////////////////////////////////////

// swap for bubble sort
void swap(int& a, int& b) {
    a = a + b;
    b = a - b;
    a = a - b;
}

// bubble sort
// descending order
// 将点的id根据所连线的数量降序排列
void bubble_sort(float* count_edge, int* vertex_desc) {
    for (int i = 0; i < VERTEX_NUM; i++) {
        for (int j = 0; j < VERTEX_NUM - i - 1; j++) {
            if (count_edge[j] < count_edge[j + 1])
            {
                swap(count_edge[j], count_edge[j + 1]);
                swap(vertex_desc[j], vertex_desc[j + 1]);
            }
        }
    }
}

//找到与当前点相连且state值为1或2的点，分两次根据向量计算新的坐标
void adjust_position_error(float** distance, float** position, int index, int* index_error) {
    int relation[VERTEX_NUM];
    for (int i = 0; i < VERTEX_NUM; i++) {
        relation[i] = -1;
    }
    int relation_num = 0;

    for (int i = 0; i < VERTEX_NUM; i++) {
        if (distance[index][i] != INFINITE && index != i) {
            relation[relation_num] = i;
            relation_num++;
        }
    }

    float x0 = position[index][0];
    float y0 = position[index][1];
    float z0 = position[index][2];

    float x, y, z;

    for (int j = 0; j < 100; j++) {             // 循环次数
        for (int i = 0; i < relation_num; i++) {

            x = position[(relation[i])][0];
            y = position[(relation[i])][1];
            z = position[(relation[i])][2];

            float rate = 0.1;   // 每次改变幅度

            x0 = x0 + rate * (sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2)) - distance[index][(relation[i])]) * (x - x0) / sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2));
            y0 = y0 + rate * (sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2)) - distance[index][(relation[i])]) * (y - y0) / sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2));
            z0 = z0 + rate * (sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2)) - distance[index][(relation[i])]) * (z - z0) / sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2));
        }
    }

    position[index][0] = x0;
    position[index][1] = y0;
    position[index][2] = z0;

}

void re_adjust(float** position, float** distance) {
    //error count
    float total_error[VERTEX_NUM];
    for (int i = 0; i < VERTEX_NUM; i++) {
        total_error[i] = 0;
        for (int j = 0; j < VERTEX_NUM; j++) {
            if (distance[i][j] != INFINITE) {
                float x1 = position[i][0];
                float y1 = position[i][1];
                float z1 = position[i][2];
                float x2 = position[j][0];
                float y2 = position[j][1];
                float z2 = position[j][2];
                total_error[i] += abs(sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2)) - distance[i][j]);
            }
        }
    }

    //sort
    int error_index[VERTEX_NUM];
    for (int i = 0; i < VERTEX_NUM; i++) {
        error_index[i] = i;
    }
    bubble_sort(total_error, error_index);

    //re_adjust
    for (int i = 0; i < VERTEX_NUM; i++) {
        int index = error_index[i];
        //adjust_position_error
        adjust_position_error(distance, position, index, error_index);
    }
}

//输出
void output(float** position) {
    ofstream white;
    white.open("output.txt");    // output文件
    for (int i = 0; i < VERTEX_NUM; i++) {
        white << std::fixed << position[i][0] << " " << position[i][1] << " " << position[i][2] << endl;
    }
}

//计算error并返回
float evaluate(float **position, float **distance, float *output_error, float *output_edge, float *output_distance, float *output_error_by_distance)
{
    float total_error = 0;
    float total_distance = 0;
    int total_edge = 0;
    for (int i = 0; i < VERTEX_NUM; i++)
    {
        for (int j = i; j < VERTEX_NUM; j++)
        {
            if (distance[i][j] != INFINITE)
            {
                float x1 = position[i][0];
                float y1 = position[i][1];
                float z1 = position[i][2];
                float x2 = position[j][0];
                float y2 = position[j][1];
                float z2 = position[j][2];

                total_distance += distance[i][j];
                total_error += abs(sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2)) - distance[i][j]);
                total_edge++;
            }
        }
    }
    // printf("Total error: %f\n", total_error); // output文件中所有连线的误差总和
    // printf("Total edge: %d\n", total_edge);
    // printf("Total distance: %f\n", total_distance); // input文件中所有连线的距离总和
    // printf("total_error / total_edge: %f\n", total_error / total_edge);
    // printf("total_error / total_distance: %f\n", total_error / total_distance);

    *output_error = total_error;
    *output_edge = total_edge;
    *output_distance = total_distance;
    *output_error_by_distance = total_error / total_distance;

    return total_error / total_edge;
}

int main(int argc, char* argv[]) {
    // auto-test
    ofstream file;
    file.open("output_combination.csv", ios::out | ios::trunc);

    if (file)
    {
        file << "file"
             << ","
             << "radius"
             << ","
             << "total error"
             << ","
             << "total edge"
             << ","
             << "total distance"
             << ","
             << "error / edge"
             << ","
             << "error / distance"
             << ","
             << "time_span"
             << endl;
    }

    for (int count = AOTO_READ_FILE_START; count <= 100; count++) // 改
    {
        clock_t start, end, time_span;
        start = clock();
        //get file name
        string str;
        ostringstream oss;
        oss << "P" << count << "-" << VERTEX_NUM << ".txt"; // file path
        str = oss.str();

        char file_name[100];

        strcpy(file_name, str.c_str());

        cout << "file name is : " << file_name << endl;

        //start per file

        // distance
        float distance[VERTEX_NUM][VERTEX_NUM];
        float* p2distance[VERTEX_NUM];
        for (int i = 0; i < VERTEX_NUM; i++)
        {
            p2distance[i] = distance[i];
        }
        read_file(p2distance, file_name);
        
        // position
        float position[VERTEX_NUM][3];
        float* p2position[VERTEX_NUM];
        for (int i = 0; i < VERTEX_NUM; i++)
        {
            p2position[i] = position[i];
        }

        initialize_positions(p2position);

        // initial error
        float output_error = 0, output_edge = 0, output_distance = 0, output_error_by_distance = 0;
        float error_rate = evaluate(p2position, p2distance, &output_error, &output_edge, &output_distance, &output_error_by_distance); // edge

        // adjust positions
        float err_edge = error_rate;
        int count_iteration = 0;
        while (true) {
            re_adjust(p2position, p2distance);
            error_rate = evaluate(p2position, p2distance, &output_error, &output_edge, &output_distance, &output_error_by_distance); // edge
            if (err_edge - error_rate < 0.01) {
                ///////////////////////////////
                end = clock();
                time_span = (float)(end - start) * 1000 / CLOCKS_PER_SEC;
                file << fixed << setprecision(6) << file_name << "," << EDGE_LENGTH << "," << output_error << "," << output_edge << "," << output_distance << "," << error_rate << "," << output_error_by_distance << "," << time_span << endl;
                cout << "Output file write successfully." << endl;
                ///////////////////////////////
                break;
            }
            err_edge = error_rate;
            count_iteration ++;
        }

        // output
        // output(p2position);
    }
    // cout << "Number of iterations is " << count_iteration << endl;

    return 0;
}