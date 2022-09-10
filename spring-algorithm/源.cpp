#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <iostream>
using namespace std;

#define INFINITE 65535
#define MAXSIZE 100
#define edge_length 50
#define pi 3.1415926535897932384626433
#define K 1

typedef struct graph
{
    // int *vexs;                       // 顶点集合
    float* vexs;
    int vexnum;                         // 顶点数
    // double matrix[MAXSIZE][MAXSIZE]; // 邻接矩阵，在二维数组中，INFINITE代表两个顶点没有连接，数字代表两顶点间的权值
    float matrix[MAXSIZE][MAXSIZE];
}Graph;

//从txt文件读入顶点和边,并创建一个无向图
void CreateGraph(Graph* graph) {

    float node[2 * MAXSIZE];            // 连续两个数字为一组顶点
    graph->vexs = node;                 // graph中顶点集合指向node数组
    graph->vexnum = 0;                  // 初始化顶点数

    float x, y;                         // 顶点坐标
    int first, second;                  // 两个顶点编号
    float weight;                       // 两个顶点之间距离

    FILE* fp1, * fp2 = NULL;

    fp1 = fopen("N3-10.txt", "r");    // node.txt

    while (fscanf(fp1, "%f %f", &x, &y) != EOF) {
        node[2 * graph->vexnum] = x;
        node[2 * graph->vexnum + 1] = y;
        //        printf("%.2f %.2f\n", x, y);
        graph->vexnum++;
    }

    fclose(fp1);

    fp2 = fopen("P3-.txt", "r");  // path.txt

    for (int i = 0; i < graph->vexnum; i++) {
        for (int j = 0; j < graph->vexnum; j++) {
            graph->matrix[i][j] = INFINITE;
        }
    }

    while (fscanf(fp2, "%d %d %f", &first, &second, &weight) != EOF) {
        graph->matrix[first - 1][second - 1] = weight;
        graph->matrix[second - 1][first - 1] = weight;
    }

    fclose(fp2);
}

//弗洛伊德算法求任意两点的最短距离（多源最短路径）d_ij
void Floyd(Graph* graph, double** shortestPath) {
    //void Floyd(Graph* graph){

    double shortestPath_copy[MAXSIZE][MAXSIZE];

    for (int i = 0; i < graph->vexnum; i++) {
        for (int j = 0; j < graph->vexnum; j++) {
            shortestPath_copy[i][j] = graph->matrix[i][j];
        }
    }

    for (int k = 0; k < graph->vexnum; k++) {         //k为中间节点
        for (int i = 0; i < graph->vexnum; i++) {
            for (int j = 0; j < graph->vexnum; j++) {
                if (shortestPath_copy[i][j] > (shortestPath_copy[i][k] + shortestPath_copy[k][j])) {
                    shortestPath_copy[i][j] = shortestPath_copy[i][k] + shortestPath_copy[k][j];      //更新最短路径
                }
            }
        }
    }

    for (int i = 0; i < graph->vexnum; i++) {
        for (int j = 0; j < graph->vexnum; j++) {
            shortestPath[i][j] = shortestPath_copy[i][j];
        }
    }
}


//根据paper中的公式（2）算出l_ij, l_ij = L * d_ij
//两点i和j之间的l_ij == the desirable length between them in the drawing
//L is the desirable length of a single edge in the display plane, L = L_0 / max(d_ij)
//L_0 is the length of a side of display square area
void compute_desirableLength(double** shortestPath, double** l, Graph graph) {
    double L_0 = 100;              // -50 < x < 50, -50 < y < 50
    double max_d_ij = 0;
    for (int i = 0; i < graph.vexnum; i++) {
        for (int j = 0; j < graph.vexnum; j++) {
            if (shortestPath[i][j] > max_d_ij) {
                max_d_ij = shortestPath[i][j];
            }
        }
    }

    double L;
    L = L_0 / max_d_ij;

    for (int i = 0; i < graph.vexnum; i++) {
        for (int j = 0; j < graph.vexnum; j++) {
            l[i][j] = L * shortestPath[i][j];
        }
    }
}

//k_ij is the strength of the spring between P_i and P_j, k_ij = K / d_ij^2,  K is a constant
void compute_springStrength(double** shortestPath, double** k, Graph graph) {
    for (int i = 0; i < graph.vexnum; i++) {
        for (int j = 0; j < graph.vexnum; j++) {
            k[i][j] = K / shortestPath[i][j] / shortestPath[i][j];
        }
    }
}

//初始化每个顶点的坐标，把顶点放置在一个直径为L_0的圆内接正n多边形上
void initialize_positions(double* xCoordinate, double* yCoordinate, int num) {
    double angle = 0;
    for (int i = 0; i < num; i++) {
        xCoordinate[i] = edge_length * cos(angle);
        yCoordinate[i] = edge_length * sin(angle);
        angle += 2 * pi / num;
    }
}

//根据公式（5）计算出当前的总弹簧能量
double calculateEnergy(double** k, double** l, double* xCoordinate, double* yCoordinate, int num) {
    double e = 0;
    for (int i = 0; i < num - 1; i++) { // num
        for (int j = i + 1; j < num; j++) {
            double xDif = pow(xCoordinate[i] - xCoordinate[j], 2);
            double yDif = pow(yCoordinate[i] - yCoordinate[j], 2);
            e += 0.5 * k[i][j] * (xDif + yDif + pow(l[i][j], 2) - 2 * l[i][j] * sqrt(xDif + yDif)); // eq 5
        }
    }
    return e;
}

//根据公式（9）计算出 Delta_m, 不断迭代的过程中需要不停的调用这个函数进行计算，来判断是否达到local minimum
double calculateDeltaM(double** k, double** l, double* xCoordinate, double* yCoordinate, int m, int num) {
    // partial derivative
    double xPd = 0, yPd = 0;
    for (int i = 0; i < num; i++) {
        if (m != i) {
            double xDif = xCoordinate[m] - xCoordinate[i], yDif = yCoordinate[m] - yCoordinate[i];
            double denominator = pow(pow(xDif, 2) + pow(yDif, 2), 0.5);
            xPd += k[m][i] * (xDif - (l[m][i] * xDif) / denominator); // eq 7
            yPd += k[m][i] * (yDif - (l[m][i] * yDif) / denominator); // eq 8
        }
    }
    double deltam = sqrt(pow(xPd, 2) + pow(yPd, 2)); // eq 9
    return deltam;
}

//初始化之后开始迭代减小E，直到达到最小的能量E
void minimizeEnergy(double** k, double** l, double* xCoordinate, double* yCoordinate, int num) {
    double e = 0.1;
    while (true) {
        // energy
        double energy = calculateEnergy(k, l, xCoordinate, yCoordinate, num);

        printf("energy %f", energy);

        printf("\n\n\n");

        // maxDeltai
        int m = 0;
        double maxDeltai = calculateDeltaM(k, l, xCoordinate, yCoordinate, 0, num);
        for (int i = 1; i < num; i++) {
            double temp = calculateDeltaM(k, l, xCoordinate, yCoordinate, i, num);
            if (maxDeltai < temp) {
                maxDeltai = temp;
                m = i;
            }
        }

        if (maxDeltai <= e) {
            break;
        }

        while (true) {
            //deltam of pm
            double deltam = calculateDeltaM(k, l, xCoordinate, yCoordinate, m, num);
            if (deltam <= e) {
                break;
            }

            // eq 13-16
            double param13 = 0, param14 = 0, param15 = 0, param16 = 0;
            double x_energy = 0.0, y_energy = 0.0;
            for (int i = 0; i < num; i++) {
                if (m != i) {
                    double xDif = xCoordinate[m] - xCoordinate[i], yDif = yCoordinate[m] - yCoordinate[i];
                    double denominator = pow(pow(xDif, 2) + pow(yDif, 2), 1.5);
                    param13 += k[m][i] * (1 - l[m][i] * pow(yDif, 2) / denominator);
                    param14 += k[m][i] * l[m][i] * xDif * yDif / denominator;
                    param15 += k[m][i] * l[m][i] * xDif * yDif / denominator;
                    param16 += k[m][i] * (1 - l[m][i] * pow(xDif, 2) / denominator);

                    double distance = pow(pow(xDif, 2) + pow(yDif, 2), 0.5);
                    x_energy += xDif * k[m][i] * (1.0 - l[m][i] / distance); // eq7
                    y_energy += yDif * k[m][i] * (1.0 - l[m][i] / distance); // eq8
                }
            }

            // eq 11-12
            double xUncertainty, yUncertainty;
            // xUncertainty = (param14 / param13 - 1) / (1 - param14 * param15 / (param13 * param16));
            // yUncertainty = -(1 + param15 / param16 * xUncertainty);
            double denom = param13 * param16 - param14 * param15;
            xUncertainty = (param14 * y_energy - param16 * x_energy) / denom; // eq11
            yUncertainty = (param14 * x_energy - param13 * y_energy) / denom; // eq12

            // eq 10
            xCoordinate[m] += xUncertainty;
            yCoordinate[m] += yUncertainty;
        }
    }
}

//输出横纵坐标到txt文件中, 按照顶点的ID排序，从1开始
void outputCoordinate(double* xCoordinate, double* yCoordinate, int num) {
    ofstream white;
    white.open("output.txt");
    for (int i = 0; i < num; i++) {
        white << std::fixed << xCoordinate[i] << "  " << yCoordinate[i] << endl;
    }
}

int main(int argc, const char* argv[]) {
    Graph graph;
    CreateGraph(&graph);                                //根据文件输入创建一个无向图，记录顶点和边的信息

    int vnum = graph.vexnum; // the number of vertices in graph

    double shortestPath[MAXSIZE][MAXSIZE];
    double* p2shortestPath[MAXSIZE];
    for (int i = 0; i < graph.vexnum; i++) {
        p2shortestPath[i] = shortestPath[i];
    }

    double l[MAXSIZE][MAXSIZE];              //desirable length between P_i and P_j in the drawing
    double* p2l[MAXSIZE];
    for (int i = 0; i < graph.vexnum; i++) {
        p2l[i] = l[i];
    }

    double k[MAXSIZE][MAXSIZE];              //strength of the spring between P_i and P_j
    double* p2k[MAXSIZE];
    for (int i = 0; i < graph.vexnum; i++) {
        p2k[i] = k[i];
    }

    double xCoordinate[MAXSIZE];            //顶点X坐标
    double yCoordinate[MAXSIZE];            //顶点y坐标

    Floyd(&graph, p2shortestPath);                       //Floyd算法求出每对点间的最短距离 d_ij

    compute_desirableLength(p2shortestPath, p2l, graph);         //计算所有顶点间的l_ij,

    compute_springStrength(p2shortestPath, p2k, graph);          //计算所有顶点间弹簧的strength
    initialize_positions(xCoordinate, yCoordinate, vnum);   //初始化每个顶点的坐标
    minimizeEnergy(p2k, p2l, xCoordinate, yCoordinate, vnum); //移动每一个顶点的位置以减小总弹簧能量，最终达到最小能量
    outputCoordinate(xCoordinate, yCoordinate, vnum);       //输出横纵坐标到文件中

    return 0;
}