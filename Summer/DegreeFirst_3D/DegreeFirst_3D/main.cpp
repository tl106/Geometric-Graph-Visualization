#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cmath>

#define INFINITE 65535
#define VERTEX_NUM 20       // node个数
#define PI 3.1415926535

using namespace std;

// read file
void read_file(float** distance) {
    FILE* fp1 = NULL;
    fp1 = fopen("/Users/hanzhekaiscomputer/Desktop/P20-20.txt", "r"); // 修改文件名

    int x, y;
    float d;
    while (fscanf(fp1, "%d %d %f", &x, &y, &d) != EOF) {
        distance[y][x] = d;
        distance[x][y] = d;

    }
    fclose(fp1);
}

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


//读入文件，创建n*n的数组distance，所有元素值初始化为INFINITE，若两点之间有连线，则相关元素值设为他们的距离
//(确定vertex_relevance)
//根据每个点与其他点的连线个数降序排列，在vertex_desc数组中记录排名
//找到连线最多的点，将其在vertex_state中的值改为1
//在position中将该点的坐/.,,标定为（0，0）
void init(float** distance, int* vertex_desc, int* vertex_state, float** position) {

    for (int i = 0; i < VERTEX_NUM; i++) {
        for (int j = 0; j < VERTEX_NUM; j++) {
            distance[i][j] = INFINITE;
        }
    }

    for (int i = 0; i < VERTEX_NUM; i++) {
        vertex_state[i] = 0;
    }

    read_file(distance);

    float count_edge[VERTEX_NUM];
    for (int m = 0; m < VERTEX_NUM; m++) {
        count_edge[m] = 0;
    }

    for (int i = 0; i < VERTEX_NUM; i++) {
        for (int j = 0; j < VERTEX_NUM; j++) {
            if (distance[i][j] != INFINITE) {
                count_edge[i] ++;
            }
        }
    }

    //    // test
    //    for (int i = 0; i < VERTEX_NUM; i++) {
    //        printf("%d\n", count_edge[i]);
    //    }

    for (int k = 0; k < VERTEX_NUM; k++) {
        vertex_desc[k] = k;
    }

    bubble_sort(count_edge, vertex_desc);

    for (int l = 0; l < VERTEX_NUM; l++) {
        if (l == 0) {
            vertex_state[0] = 1;
        }
        else {
            vertex_state[l] = 0;
        }
    }

    //    // test
    //    printf("\n\n\n");
    //    for (int i = 0; i < VERTEX_NUM; i++) {
    //        printf("%d\n", vertex_desc[i]);
    //    }
    //    printf("\n\n\n");

    position[(vertex_desc[0])][0] = 0;
    position[(vertex_desc[0])][1] = 0;
    position[(vertex_desc[0])][2] = 0;
}

int find_state_by_index(int index, int* vertex_state, int* vertex_desc) {
    for (int i = 0; i < VERTEX_NUM; i++) {
        if (vertex_desc[i] == index) {
            return vertex_state[i];
        }
    }
    return -1;
}

//找到与当前点相连且state值为1或2的点，分两次根据向量计算新的坐标
void adjust_position(float** distance, float** position, int index, int* vertex_state, int* vertex_desc) {
    int relation[VERTEX_NUM];
    for (int i = 0; i < VERTEX_NUM; i++) {
        relation[i] = -1;
    }
    int relation_num = 0;

    for (int i = 0; i < VERTEX_NUM; i++) {
        if (distance[index][i] != INFINITE && index != i) {
            if (find_state_by_index(i, vertex_state, vertex_desc) > 0) {
                relation[relation_num] = i;
                relation_num++;
            }
        }
    }

    float x0 = position[index][0];
    float y0 = position[index][1];
    float z0 = position[index][2];

    //    // test
    //    printf("%f, %f\n", x0, y0);

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

    //    // test
    //    printf("%f, %f\n", x0, y0);

}


//按照当前点state为0的子节点，随机分配在当前点的周围，距离以实际距离为准。
//将这些点的state改为1，并将position中这些点的坐标写好
//将当前点的state改为2
void expand(int* vertex_state, int index, float** position, float** distance, int* vertex_desc) {
    int constant = 123456;
    float var1 = constant;
    float var2 = constant;
    //    float var = 654321;
    for (int i = 0; i < VERTEX_NUM; i++) {
        float idistance = distance[index][i];
        if (idistance != INFINITE && find_state_by_index(i, vertex_state, vertex_desc) == 0) {
            var1 += constant / 321 * 123456;
            //            var1 += 654321 / 321 * 123456;
            srand(var1);
            float alpha = (float)((rand() % 100000) * 0.00001);

            var2 += constant / 321 * 654321;
            //            var2 += 654321 / 321 * 654321;
            srand(var2);
            float beta = (float)((rand() % 100000) * 0.00001);

            alpha = 2 * PI * alpha - PI;
            beta = 2 * PI * beta - PI;

            float dx = position[index][0];
            float dy = position[index][1];
            float dz = position[index][2];

            position[i][0] = dx + idistance * sin(alpha) * cos(beta);
            position[i][1] = dy + idistance * sin(alpha) * sin(beta);
            position[i][2] = dz + idistance * cos(alpha);


            for (int j = 0; j < VERTEX_NUM; j++) {
                if (vertex_desc[j] == i) {
                    vertex_state[j] = 1;
                }
            }
        }
    }
    for (int j = 0; j < VERTEX_NUM; j++) {
        if (vertex_desc[j] == index) {
            vertex_state[j] = 2;
        }
    }
}


//找到第一个state为1的点，并返回它的index
int getCurrent(int* vertex_state, int* vertex_desc) {
    for (int i = 0; i < VERTEX_NUM; i++) {
        if (vertex_state[i] == 1) {
            return vertex_desc[i];
        }
    }
    return -1;
}

//计算error并返回
float evaluate(float** position, float** distance) {
    float total_error = 0;
    float total_distance = 0;
    int total_edge = 0;
    for (int i = 0; i < VERTEX_NUM; i++) {
        for (int j = i; j < VERTEX_NUM; j++) {
            if (distance[i][j] != INFINITE) {
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
    //printf("Total distance: %f\n", total_distance);     // input文件中所有连线的距离总和
    printf("Total error: %f\n", total_error);           // output文件中所有连线的误差总和
    printf("Total edge: %d\n", total_edge);
    printf("Average error: %f\n", total_error / total_edge);

    return total_error / total_edge;
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

    //    // test
    //    printf("%f, %f\n", x0, y0);

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

    //    // test
    //    printf("%f, %f\n", x0, y0);

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
    white.open("/Users/hanzhekaiscomputer/Desktop/P20-20-output.txt");    // output文件
    for (int i = 0; i < VERTEX_NUM; i++) {
        white << std::fixed << position[i][0] << " " << position[i][1] << " " << position[i][2] << endl;
    }
}


int main(int argc, char* argv[]) {

    float distance[VERTEX_NUM][VERTEX_NUM];
    float* p2distance[VERTEX_NUM];
    for (int i = 0; i < VERTEX_NUM; i++)
    {
        p2distance[i] = distance[i];
    }

    //0 未处理的点; 1 已出现在图上但未确定位置的点; 2 已确定位置的点
    //首先需要全部初始化为0
    int vertex_state[VERTEX_NUM];

    //记录每个点与其他点连线的条数
//    int vertex_relevance[VERTEX_NUM];

    //根据每个点与其他点的连线个数降序排列，在vertex_desc数组中记录排名
    int vertex_desc[VERTEX_NUM];

    float position[VERTEX_NUM][3];
    float* p2position[VERTEX_NUM];
    for (int i = 0; i < VERTEX_NUM; i++)
    {
        p2position[i] = position[i];
    }

    //    int current_vertex;

    init(p2distance, vertex_desc, vertex_state, p2position);

    for (int i = 0; i < VERTEX_NUM; i++) {
        int index = getCurrent(vertex_state, vertex_desc);
        if (index == -1) {
            break;
        }
        adjust_position(p2distance, p2position, index, vertex_state, vertex_desc);

        expand(vertex_state, index, p2position, p2distance, vertex_desc);
    }

    float error_rate = evaluate(p2position, p2distance); // edge

    float err_edge = error_rate;
    int count_iteration = 0;
    while (true) {
        re_adjust(p2position, p2distance);
        cout << "------------" << endl;
        error_rate = evaluate(p2position, p2distance); // edge
        if (err_edge - error_rate < 0.01) {
            break;
        }
        err_edge = error_rate;
        count_iteration ++;
    }

    cout << "Number of iterations is " << count_iteration << endl;
    
    output(p2position);

    return 0;
}

