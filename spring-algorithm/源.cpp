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
    // int *vexs;                       // ���㼯��
    float* vexs;
    int vexnum;                         // ������
    // double matrix[MAXSIZE][MAXSIZE]; // �ڽӾ����ڶ�ά�����У�INFINITE������������û�����ӣ����ִ�����������Ȩֵ
    float matrix[MAXSIZE][MAXSIZE];
}Graph;

//��txt�ļ����붥��ͱ�,������һ������ͼ
void CreateGraph(Graph* graph) {

    float node[2 * MAXSIZE];            // ������������Ϊһ�鶥��
    graph->vexs = node;                 // graph�ж��㼯��ָ��node����
    graph->vexnum = 0;                  // ��ʼ��������

    float x, y;                         // ��������
    int first, second;                  // ����������
    float weight;                       // ��������֮�����

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

//���������㷨�������������̾��루��Դ���·����d_ij
void Floyd(Graph* graph, double** shortestPath) {
    //void Floyd(Graph* graph){

    double shortestPath_copy[MAXSIZE][MAXSIZE];

    for (int i = 0; i < graph->vexnum; i++) {
        for (int j = 0; j < graph->vexnum; j++) {
            shortestPath_copy[i][j] = graph->matrix[i][j];
        }
    }

    for (int k = 0; k < graph->vexnum; k++) {         //kΪ�м�ڵ�
        for (int i = 0; i < graph->vexnum; i++) {
            for (int j = 0; j < graph->vexnum; j++) {
                if (shortestPath_copy[i][j] > (shortestPath_copy[i][k] + shortestPath_copy[k][j])) {
                    shortestPath_copy[i][j] = shortestPath_copy[i][k] + shortestPath_copy[k][j];      //�������·��
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


//����paper�еĹ�ʽ��2�����l_ij, l_ij = L * d_ij
//����i��j֮���l_ij == the desirable length between them in the drawing
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

//��ʼ��ÿ����������꣬�Ѷ��������һ��ֱ��ΪL_0��Բ�ڽ���n�������
void initialize_positions(double* xCoordinate, double* yCoordinate, int num) {
    double angle = 0;
    for (int i = 0; i < num; i++) {
        xCoordinate[i] = edge_length * cos(angle);
        yCoordinate[i] = edge_length * sin(angle);
        angle += 2 * pi / num;
    }
}

//���ݹ�ʽ��5���������ǰ���ܵ�������
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

//���ݹ�ʽ��9������� Delta_m, ���ϵ����Ĺ�������Ҫ��ͣ�ĵ�������������м��㣬���ж��Ƿ�ﵽlocal minimum
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

//��ʼ��֮��ʼ������СE��ֱ���ﵽ��С������E
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

//����������굽txt�ļ���, ���ն����ID���򣬴�1��ʼ
void outputCoordinate(double* xCoordinate, double* yCoordinate, int num) {
    ofstream white;
    white.open("output.txt");
    for (int i = 0; i < num; i++) {
        white << std::fixed << xCoordinate[i] << "  " << yCoordinate[i] << endl;
    }
}

int main(int argc, const char* argv[]) {
    Graph graph;
    CreateGraph(&graph);                                //�����ļ����봴��һ������ͼ����¼����ͱߵ���Ϣ

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

    double xCoordinate[MAXSIZE];            //����X����
    double yCoordinate[MAXSIZE];            //����y����

    Floyd(&graph, p2shortestPath);                       //Floyd�㷨���ÿ�Ե�����̾��� d_ij

    compute_desirableLength(p2shortestPath, p2l, graph);         //�������ж�����l_ij,

    compute_springStrength(p2shortestPath, p2k, graph);          //�������ж���䵯�ɵ�strength
    initialize_positions(xCoordinate, yCoordinate, vnum);   //��ʼ��ÿ�����������
    minimizeEnergy(p2k, p2l, xCoordinate, yCoordinate, vnum); //�ƶ�ÿһ�������λ���Լ�С�ܵ������������մﵽ��С����
    outputCoordinate(xCoordinate, yCoordinate, vnum);       //����������굽�ļ���

    return 0;
}