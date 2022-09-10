#include <stdio.h>
#include <stdlib.h>

#define INFINITE 65535
#define MAXSIZE 100


typedef struct graph        
{
    int *vexs              // 顶点集合
    int vexnum;            // 顶点数
    double matrix[MAXSIZE][MAXSIZE]; // 邻接矩阵，在二维数组中，INFINITE代表两个顶点没有连接，数字代表两顶点间的权值
}Graph

//从txt文件读入顶点和边,并创建一个无向图
void CreateGraph(Graph* graph){  
}

//弗洛伊德算法求任意两点的最短距离（多源最短路径）d_ij
void Floyd(Graph* graph, double* shortestPath){       
    
    for(int i = 0; i < graph->vexnum; i++){
        for(int j = 0; j < graph->vexnum; j++){
            (*shortestPath)[i][j] = graph->matrix[i][j];
        }     
    }
        
    for(int k; k < graph->vexnum; k++){         //k为中间节点
        for(int i = 0; i < graph->vexnum; i++){
            for(int j = 0; j < graph->vexnum; j++){
                if((*shortestPath)[i][j] > (*shortestPath)[i][k] + (*shortestPath)[k][j]){
                    (*shortestPath)[i][j] = (*shortestPath)[i][k] + (*shortestPath)[k][j];      //更新最短路径
                }                 
            }                
        }                         
    }       
}


//根据paper中的公式（2）算出l_ij, l_ij = L * d_ij 
//两点i和j之间的l_ij == the desirable length between them in the drawing
//L is the desirable length of a single edge in the display plane, L = L_0 / max(d_ij)
//L_0 is the length of a side of display square area
void compute_desirableLength(double* shortestPath, double* l){                                                     
}

//k_ij is the strength of the spring between P_i and P_j, k_ij = K / d_ij^2,  K is a constant
void compute_springStrength(double* k){
}

//初始化每个顶点的坐标，把顶点放置在一个直径为L_0的圆内接正n多边形上
void initialize_positions(double* xCoordinate, double* yCoordinate){

}

//根据公式（5）计算出当前的总弹簧能量
void calculateEnergy(){

}

//根据公式（9）计算出 Delta_m, 不断迭代的过程中需要不停的调用这个函数进行计算，来判断是否达到local minimum
void calculateDeltaM(){

}

//初始化之后开始迭代减小E，直到达到最小的能量E
void minimizeEnergy(double* k, double* l, double* xCoordinate, double* yCoordinate){
    while(){
        while(){

        }
    }
    //calculateEnergy();
    //calculateDeltaM();
} 

//输出横纵坐标到txt文件中, 按照顶点的ID排序，从1开始
void outputCoordinate(double* xCoordinate, double* yCoordinate){

}

int void main(int argc,char argc[])
{
    Graph graph;
    double shortestPath[graph.vexnum][graph.vexnum];
    double l[graph.vexnum][graph.vexnum];              //desirable length between P_i and P_j in the drawing
    double k[graph.vexnum][graph.vexnum];              //strength of the spring between P_i and P_j
    double xCoordinate[graph.graph.vexnum];            //顶点X坐标
    double yCoordinate[graph.graph.vexnum];            //顶点y坐标

    CreateGraph(&graph);                                //根据文件输入创建一个无向图，记录顶点和边的信息
    Floyd(&graph, &shortestPath);                       //Floyd算法求出每对点间的最短距离 d_ij
    compute_desirableLength(&shortestPath, &l);         //计算所有顶点间的l_ij, 
    compute_springStrength(&shortestPath, &k);          //计算所有顶点间弹簧的strength
    initialize_positions(&xCoordinate, &yCoordinate);   //初始化每个顶点的坐标
    minimizeEnergy(&k, &l, &xCoordinate, &yCoordinate); //移动每一个顶点的位置以减小总弹簧能量，最终达到最小能量
    outputCoordinate(&xCoordinate, &yCoordinate);       //输出横纵坐标到文件中
}
