#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <cmath>

using namespace std;

#define INFINITE 65535
#define MAXSIZE 100
#define EDGE_LENGTH 50
#define PI 3.1415926535897932384626433
#define K 1

typedef struct graph {
    int vexnum;                                 // ������
    float matrix[MAXSIZE][MAXSIZE];             // �ڽӾ����ڶ�ά�����У�INFINITE������������û�����ӣ����ִ�����������Ȩֵ
}Graph;

// ��txt�ļ����붥��ͱ�,������һ������ͼ
void CreateGraph(Graph* graph) {
    graph->vexnum = 0;                          // ��ʼ��������

    int first, second;                          // ����������
    float weight;                               // ��������֮�����

    FILE* fp1 = NULL;
    if ((fp1 = fopen("P3-20.txt", "r")) == NULL) { // ��path.txt�ļ�
        cout << "Fail to open the file." << endl;
        exit(-1);
    }

    // ��matrix������ֵ��ʼ��ΪINFINITE
    for (int i = 0; i < MAXSIZE; i++) {
        for (int j = 0; j < MAXSIZE; j++) {
            graph->matrix[i][j] = INFINITE;
        }
    }

    // ��ȡ�ļ��������ӵ�������λ������matrix
    // ȷ��vexnum�����������
    while (fscanf(fp1, "%d %d %f", &first, &second, &weight) != EOF) {
        graph->matrix[first][second] = weight;
        graph->matrix[second][first] = weight;

        if (first > second && first > graph->vexnum) {
            graph->vexnum = first;
        }
        else if (second > first && second > graph->vexnum) {
            graph->vexnum = second;
        }
    }
    graph->vexnum++;

    fclose(fp1);

    // debug
    /*cout << "vexnum: " << graph->vexnum << endl;
    cout << "matrix:\n";
    for (int i = 0; i < graph->vexnum; i++) {
        for (int j = 0; j < graph->vexnum; j++) {
            cout << graph->matrix[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;*/
}

// ���������㷨�������������̾��루��Դ���·����d_ij
void Floyd(Graph* graph, double** d) {
    // ��ʼ��
    for (int i = 0; i < graph->vexnum; i++) {
        for (int j = 0; j < graph->vexnum; j++) {
            d[i][j] = graph->matrix[i][j];
        }
    }

    for (int k = 0; k < graph->vexnum; k++) {       // kΪ�м�ڵ�
        for (int i = 0; i < graph->vexnum; i++) {
            for (int j = 0; j < graph->vexnum; j++) {
                if (d[i][j] > (d[i][k] + d[k][j])) {
                    d[i][j] = d[i][k] + d[k][j];    // �������·��
                }
            }
        }
    }

    // debug
    /*cout << "d:\n";
    for (int i = 0; i < graph->vexnum; i++) {
        for (int j = 0; j < graph->vexnum; j++) {
            cout << d[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;*/
}

// �����������ĳ���l_ij
// l_ij is the desirable length between them in the drawing
// L is the desirable length of a single edge in the display plane, L = L_0 / max(d_ij)
// L_0 is the length of a side of display square area
void compute_desirableLength(double** d, double** l, Graph graph) {
    // -50 ~ 50
    double L_0 = (double)EDGE_LENGTH * 2;

    // ��d�����ֵ������max_d_ij
    double max_d_ij = 0;
    for (int i = 0; i < graph.vexnum; i++) {
        for (int j = 0; j < graph.vexnum; j++) {
            if (d[i][j] > max_d_ij) {
                max_d_ij = d[i][j];
            }
        }
    }

    // eq3
    double L = L_0 / max_d_ij;

    // eq2
    for (int i = 0; i < graph.vexnum; i++) {
        for (int j = 0; j < graph.vexnum; j++) {
            l[i][j] = L * d[i][j];
        }
    }

    // debug
    /*cout << "l:\n";
    for (int i = 0; i < graph.vexnum; i++) {
        for (int j = 0; j < graph.vexnum; j++) {
            cout << l[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;*/
}

//k_ij is the strength of the spring between P_i and P_j, k_ij = K / d_ij^2,  K is a constant
void compute_springStrength(double** d, double** k, Graph graph) {
    for (int i = 0; i < graph.vexnum; i++) {
        for (int j = 0; j < graph.vexnum; j++) {
            k[i][j] = K / pow(d[i][j], 2); // eq4
        }
    }

    // debug
    /*cout << "k:\n";
    for (int i = 0; i < graph.vexnum; i++) {
        for (int j = 0; j < graph.vexnum; j++) {
            cout << k[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;*/
}

// ��ά�����ϵ�Marsaglia����
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

// ��ʼ��ÿ����������꣬ʹ������һ��ֱ��ΪL_0�������Ͻ�Ϊ���ȵķֲ�
void initialize_positions(double* xCordinate, double* yCordinate, double* zCordinate, int vexnum) {
    int z1 = 20;    // ����
    int z2 = 112;   // ����
    cRandom sita(z1, 0.0);
    cRandom pesi(z2, 0.0);

    int count = 0;
    while (count < vexnum) {
        sita = my_random(pesi.seed);
        pesi = my_random(sita.seed);

        double u = 2.0 * sita.random - 1.0; // -1 ~ 1
        double v = 2.0 * pesi.random - 1.0;

        double r2 = pow(u, 2) + pow(v, 2);
        if (r2 < 1) {
            // ��(x,y,z)�ڰ뾶Ϊ1��������
            double x = 2 * u * sqrt(1 - r2);
            double y = 2 * v * sqrt(1 - r2);
            double z = 1 - 2 * r2;

            xCordinate[count] = EDGE_LENGTH * x; // EDGE_LENGTH: �뾶
            yCordinate[count] = EDGE_LENGTH * y;
            zCordinate[count] = EDGE_LENGTH * z;

            // debug
            // cout << x * 50 << " " << y * 50 << " " << z * 50 << endl;

            count++;
        }
    }

    // debug
    /*for (int i = 0; i < vexnum; i++) {
        cout << xCordinate[i] << " " << yCordinate[i] << " " << zCordinate[i] << " " << endl;
    }*/
}

// ���㵱ǰ���ܵ�������
double calculateEnergy(double** k, double** l, double* xCoordinate, double* yCoordinate, double* zCoordinate, int num) {
    double energy = 0;
    for (int i = 0; i < num - 1; i++) { // num
        for (int j = i + 1; j < num; j++) {
            double square_dif_x, square_dif_y, square_dif_z;
            square_dif_x = pow(xCoordinate[i] - xCoordinate[j], 2);
            square_dif_y = pow(yCoordinate[i] - yCoordinate[j], 2);
            square_dif_z = pow(zCoordinate[i] - zCoordinate[j], 2);
            energy += 0.5 * k[i][j] * (square_dif_x + square_dif_y + square_dif_z + pow(l[i][j], 2) - 2 * l[i][j] * sqrt(square_dif_x + square_dif_y + square_dif_z)); // eq 5
        }
    }
    return energy;
}

// ����Delta_m, ���ϵ����Ĺ�������Ҫ��ͣ�ĵ�������������м��㣬���ж��Ƿ�ﵽlocal minimum
double calculateDeltaM(double** k, double** l, double* xCoordinate, double* yCoordinate, double* zCoordinate, int m, int num) {
    // partial derivative
    double energy_x = 0, energy_y = 0, energy_z = 0;
    for (int i = 0; i < num; i++) {
        if (m != i) {
            double dif_x = xCoordinate[m] - xCoordinate[i], dif_y = yCoordinate[m] - yCoordinate[i], dif_z = zCoordinate[m] - zCoordinate[i];
            double distance = sqrt(pow(dif_x, 2) + pow(dif_y, 2) + pow(dif_z, 2)); // distance between vertex m & i
            energy_x += k[m][i] * dif_x * (1 - l[m][i] / distance); // eq 7
            energy_y += k[m][i] * dif_y * (1 - l[m][i] / distance); // eq 8
            energy_z += k[m][i] * dif_z * (1 - l[m][i] / distance); // eq 9
        }
    }
    double deltam = sqrt(pow(energy_x, 2) + pow(energy_y, 2) + pow(energy_z, 2)); // eq 10
    return deltam;
}

// ��ʼ��֮��ʼ������СE��ֱ���ﵽ��С������E
void minimizeEnergy(double** k, double** l, double* xCoordinate, double* yCoordinate, double* zCoordinate, int num) {
    double e = 0.1;

    int count1 = 0;
    const int MAX_ITERATE = 50;
    while (count1 < MAX_ITERATE) {
        // energy
        double energy = calculateEnergy(k, l, xCoordinate, yCoordinate, zCoordinate, num);
        cout << "energy: " << energy << "\n\n\n";

        // maxDeltai
        int m = 0;
        double maxDeltai = calculateDeltaM(k, l, xCoordinate, yCoordinate, zCoordinate, 0, num);
        for (int i = 1; i < num; i++) {
            double temp = calculateDeltaM(k, l, xCoordinate, yCoordinate, zCoordinate, i, num);
            if (maxDeltai < temp) {
                maxDeltai = temp;
                m = i;
            }
        }

        if (maxDeltai <= e) {
            break;
        }

        int count2 = 0;
        while (count2 < MAX_ITERATE) {
            //deltam of pm
            double deltam = calculateDeltaM(k, l, xCoordinate, yCoordinate, zCoordinate, m, num);
            if (deltam <= e) {
                break;
            }

            //  eq 15 - 23
            double energy_xx = 0, energy_xy = 0, energy_xz = 0, energy_yz = 0, energy_zz = 0, energy_yy = 0;
            double energy_x = 0, energy_y = 0, energy_z = 0;
            for (int i = 0; i < num; i++) {
                if (m != i) {
                    double dif_x = xCoordinate[m] - xCoordinate[i], dif_y = yCoordinate[m] - yCoordinate[i], dif_z = zCoordinate[m] - zCoordinate[i];
                    double square_dif_x = pow(dif_x, 2), square_dif_y = pow(dif_y, 2), square_dif_z = pow(dif_z, 2);
                    double denominator = pow(square_dif_x + square_dif_y + square_dif_z, 1.5);

                    energy_xx += k[m][i] * (1 - l[m][i] * (square_dif_y + square_dif_z) / denominator); // eq 15
                    energy_zz += k[m][i] * (1 - l[m][i] * (square_dif_x + square_dif_y) / denominator); // eq 22
                    energy_yy += k[m][i] * (1 - l[m][i] * (square_dif_x + square_dif_z) / denominator); // eq 23
                    energy_xy += k[m][i] * l[m][i] * (square_dif_x + square_dif_y) / denominator; // eq 17
                    energy_xz += k[m][i] * l[m][i] * (square_dif_x + square_dif_z) / denominator; // eq 18
                    energy_yz += k[m][i] * l[m][i] * (square_dif_y + square_dif_z) / denominator; // eq 21

                    double distance = pow(pow(dif_x, 2) + pow(dif_y, 2) + pow(dif_z, 2), 0.5); // distance between vertex m & i
                    energy_x += k[m][i] * dif_x * (1 - l[m][i] / distance); // eq 7
                    energy_y += k[m][i] * dif_y * (1 - l[m][i] / distance); // eq 8
                    energy_z += k[m][i] * dif_z * (1 - l[m][i] / distance); // eq 9
                }
            }
            energy_x = -energy_x;
            energy_y = -energy_y;
            energy_z = -energy_z;

            // eq 12 - 14 (eq 24 - 26)
            double xUncertainty, yUncertainty, zUncertainty;

            // D =  afk + bgi + cej - cfi - agj - bek
            // X = dfk + bgl + chj - cfl - dgj - bhk
            // Y = ahk + dgi + cel - chi - dek - agl
            // Z = afl + bhi + dje - dfi - bel - ahj
            double D = (energy_xx * energy_yy * energy_zz) + (energy_xy * energy_yz * energy_xz) + (energy_xz * energy_xy * energy_yz) - (energy_xz * energy_yy * energy_xz) - (energy_xx * energy_yz * energy_yz) - (energy_xy * energy_xy * energy_zz);
            double X = (energy_x * energy_yy * energy_zz) + (energy_xy * energy_yz * energy_z) + (energy_xz * energy_y * energy_yz) - (energy_xz * energy_yy * energy_z) - (energy_x * energy_yz * energy_yz) - (energy_xy * energy_y * energy_zz);
            double Y = (energy_xx * energy_y * energy_zz) + (energy_x * energy_yz * energy_xz) + (energy_xz * energy_xy * energy_z) - (energy_xz * energy_y * energy_xz) - (energy_x * energy_xy * energy_zz) - (energy_xx * energy_yz * energy_z);
            xUncertainty = X / D;
            yUncertainty = Y / D;
            zUncertainty = (energy_x - energy_xx * xUncertainty - energy_xy * yUncertainty) / energy_xz;

            // eq 11
            xCoordinate[m] += xUncertainty;
            yCoordinate[m] += yUncertainty;
            zCoordinate[m] += zUncertainty;
            count2++;
        }
        count1++;
    }
}

//����������굽txt�ļ���, ���ն����ID���򣬴�1��ʼ
void outputCoordinate(double* xCoordinate, double* yCoordinate, double* zCoordinate, int num) {
    ofstream white;
    white.open("output.txt");
    for (int i = 0; i < num; i++) {
        white << std::fixed << xCoordinate[i] << " " << yCoordinate[i] << " " << zCoordinate[i] << endl;
    }
}

void evaluate(double* xCoordinate, double* yCoordinate, double* zCoordinate, Graph* graph) {
    double total_err = 0;
    int count_edge = 0;
    for (int i = 0; i < graph->vexnum; i++) {
        for (int j = i; j < graph->vexnum; j++) {
            if (graph->matrix[i][j] != INFINITE) {
                double x_dif = xCoordinate[i] - xCoordinate[j], y_dif = yCoordinate[i] - yCoordinate[j], z_dif = zCoordinate[i] - zCoordinate[j];
                double output_distance = sqrt(pow(x_dif, 2) + pow(y_dif, 2) + pow(z_dif, 2));
                double temp_err = output_distance - graph->matrix[i][j];
                if (temp_err < 0) {
                    temp_err = -temp_err;
                }
                total_err += temp_err;
                count_edge++;
            }
        }
    }
    cout << "number of edge: " << count_edge << endl;
    cout << "total error: " << total_err << endl;

    double avg_err1 = total_err / count_edge;
    cout << "average error by edge: " << avg_err1 << endl;

    double avg_err2 = total_err / graph->vexnum;
    cout << "average error by vertices: " << avg_err2 << endl;
}

int main(int argc, const char* argv[]) {
    Graph graph;
    CreateGraph(&graph);            //�����ļ����봴��һ������ͼ����¼����ͱߵ���Ϣ

    int vnum = graph.vexnum;        // ����ͼ�еĶ�������

    double d[MAXSIZE][MAXSIZE];     // shortest distance dij between vertices
    double* p2d[MAXSIZE];
    for (int i = 0; i < graph.vexnum; i++) {
        p2d[i] = d[i];
    }

    double l[MAXSIZE][MAXSIZE];     // desirable length between P_i and P_j in the drawing
    double* p2l[MAXSIZE];
    for (int i = 0; i < graph.vexnum; i++) {
        p2l[i] = l[i];
    }

    double k[MAXSIZE][MAXSIZE];     // strength of the spring between P_i and P_j
    double* p2k[MAXSIZE];
    for (int i = 0; i < graph.vexnum; i++) {
        p2k[i] = k[i];
    }

    Floyd(&graph, p2d);                         // Floyd�㷨���ÿ�Ե�����̾��� d_ij
    compute_desirableLength(p2d, p2l, graph);   // �������ж�����l_ij,
    compute_springStrength(p2d, p2k, graph);    // �������ж���䵯�ɵ�strength

    double xCoordinate[MAXSIZE];                // ����x������
    double yCoordinate[MAXSIZE];                // ����y������
    double zCoordinate[MAXSIZE];                // ����z������

    initialize_positions(xCoordinate, yCoordinate, zCoordinate, vnum);       //��ʼ��ÿ�����������
    minimizeEnergy(p2k, p2l, xCoordinate, yCoordinate, zCoordinate, vnum);   //�ƶ�ÿһ�������λ���Լ�С�ܵ������������մﵽ��С����
    outputCoordinate(xCoordinate, yCoordinate, zCoordinate, vnum);           //����������굽�ļ���
    evaluate(xCoordinate, yCoordinate, zCoordinate, &graph);

    return 0;
}