#define GL_SILENCE_DEPRECATION

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <GL/freeglut.h>
#include <iostream>
#include <cmath>
#include <iomanip>

#define WINDOW_WIDTH 800
#define WINDOW_HEIGHT 800
#define INFINITE 65535
#define MAXSIZE 1000
#define PI 3.14159265358979323846

using namespace std;

typedef struct graph
{
    float *vexs;
    int vexnum;
    float matrix[MAXSIZE][MAXSIZE];
}Graph;

GLfloat G_fDistance = 0.0f;
GLfloat G_fAngle_horizon = 0.0f;
GLfloat G_fAngle_vertical = 0.0f;
GLfloat G_fLength = 75.0f;
GLfloat G_fRate = 1.0f;
GLfloat G_fRadius_sphere = 1.0f;
GLfloat G_fRadius_cylinder = 0.25f;

void myinit(void);
void myReshape(GLsizei w, GLsizei h);
void display(void);
void processSpecialKeys(int key, int x, int y);
void processNormalKeys(unsigned char key,int x,int y);
void cyLinder(float x0, float y0, float z0, float x1, float y1, float z1, float red, float blue, float g);
float cal_radius(float length, float desirable);
static void torus(float x0, float y0, float z0, float x1, float y1, float z1, float desirable);
static void torus2(float x0, float y0, float z0, float x1, float y1, float z1, float desirable);

GLuint theTorus;

float G_vLit0Position[4] = { 100.0f, 100.0f, 100.0f, 1.0f };
float G_vLit0Ambient[4] = { 0.2f, 0.2f, 0.2f, 1.0f };
float G_vLit0Diffuse[4] = { 0.8f, 0.8f, 0.8f, 1.0f };
float G_vLit0Specular[4] = { 0.5f, 0.0f, 0.5f, 1.0f };
//float G_vMaterialSpecu[4] = { 0.0f, 0.0f, 0.0f, 1.0f };
//
float G_vLit1Position[4] = { -100.0f, -100.0f, -100.0f, 1.0f };
float G_vLit1Ambient[4] = { 0.2f, 0.2f, 0.2f, 1.0f };
float G_vLit1Diffuse[4] = { 0.8f, 0.8f, 0.8f, 1.0f };
float G_vLit1Specular[4] = { 0.5f, 0.0f, 0.5f, 1.0f };

GLfloat mat_specular[] = { 1.0, 0.0, 0.0, 1.0 };
GLfloat mat_shininess[] = { 50.0 };
GLfloat lmodel_ambient[] = { 1.0, 0.0 ,0.0 ,1.0 }; //太阳颜色为红色
GLfloat white_light[] = { 1.0, 1.0, 1.0, 1.0 };



int main(int argc, char* argv[])
{
    glutInit(&argc, argv);
    
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);

    glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);
    glutInitWindowPosition(100, 100);
    
    glutCreateWindow("OpenGL");
    
    myinit();

    glutReshapeFunc(myReshape);

    glutSpecialFunc(processSpecialKeys);
    glutKeyboardFunc(processNormalKeys);
    
    glutDisplayFunc(display);

    glutIdleFunc(display);

    glutMainLoop();

    return 0;
    
}

void myinit(void)
{
    glEnable(GL_DEPTH_TEST);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHT1);
    glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, 1.0);

    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
}

void myReshape(GLsizei w, GLsizei h)
{
//    glViewport(0, 0, w, h);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
//    gluPerspective(120, 1.0*(GLfloat)w/(GLfloat)h, 1.0, 300.0);
    
    glOrtho(-G_fLength, G_fLength, -G_fLength, G_fLength, 500.0, -500.0);
    glutSwapBuffers();
}

void display(void)
{
    Graph graph;
    float node[3*MAXSIZE];             // 连续两个数字为一组顶点
    graph.vexs = node;                 // graph中顶点集合指向node数组
    graph.vexnum = 0;                  // 初始化顶点数
    
    float x, y, z;                      // 顶点坐标
    int first, second;                  // 两个顶点编号
    float weight;                       // 两个顶点之间距离
    
    FILE *fp1, *fp2 = NULL;
    
    fp1 = fopen("/Users/hanzhekaiscomputer/Desktop/P20-20-output.txt", "r");    // node.txt
    
    while (fscanf(fp1, "%f %f %f", &x, &y, &z) != EOF) {
        node[3*graph.vexnum] = x;
        node[3*graph.vexnum+1] = y;
        node[3*graph.vexnum+2] = z;
        graph.vexnum++;
    }
//    for (int i = 0; i<3*graph.vexnum; i = i+3) {
//        printf("%.2f %.2f %.2f\n", node[i], node[i+1], node[i+2]);
//    }
//    printf("%d\n", graph.vexnum);
    
    fclose(fp1);
    
    fp2 = fopen("/Users/hanzhekaiscomputer/Desktop/P20-20.txt", "r");  // path.txt
    
    for (int i =0; i < graph.vexnum; i++) {
        for (int j = 0; j < graph.vexnum; j++) {
            graph.matrix[i][j] = INFINITE;
        }
    }
    
    while (fscanf(fp2, "%d %d %f", &first, &second, &weight) != EOF) {
        graph.matrix[first][second] = weight;
        graph.matrix[second][first] = weight;
    }

//    for (int i = 0; i < graph.vexnum; i++) {
//        for (int j = 0; j < graph.vexnum; j++) {
//            printf("%.2f ", graph.matrix[i][j]);
//        }
//        printf("\n");
//    }
    
    fclose(fp2);
    
    
    
    glClearColor(0.0f,0.0f,0.0f,0.0f);
//    glClearColor(1.0f,1.0f,1.0f,0.0f);
    glClearDepth(1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    glDisable(GL_LIGHTING);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glEnable(GL_LIGHTING);
    
//    glMatrixMode(GL_MODELVIEW);
//    glLoadIdentity();
//    glLightfv(GL_LIGHT0, GL_POSITION, G_vLit0Position);
//    glLightfv(GL_LIGHT0, GL_AMBIENT, G_vLit0Ambient);
//    glLightfv(GL_LIGHT0, GL_DIFFUSE, G_vLit0Diffuse);
//    glLightfv(GL_LIGHT0, GL_SPECULAR, G_vLit0Specular);
//
//    glMatrixMode(GL_MODELVIEW);
//    glLoadIdentity();
//    glLightfv(GL_LIGHT1, GL_POSITION, G_vLit1Position);
//    glLightfv(GL_LIGHT1, GL_AMBIENT, G_vLit1Ambient);
//    glLightfv(GL_LIGHT1, GL_DIFFUSE, G_vLit1Diffuse);
//    glLightfv(GL_LIGHT1, GL_SPECULAR, G_vLit1Specular);
    
    glMatrixMode(GL_MODELVIEW);
    
    glLoadIdentity();
    glTranslatef(0.0, 0.0, -G_fDistance);
    glRotatef(G_fAngle_horizon, 0.0, 1.0, 0.0);
    glRotatef(G_fAngle_vertical, 1.0, 0.0, 0.0);
    
    
    
//    glPushMatrix();
//
//    glTranslatef(0.0, 0.0, -G_fDistance);
//    glRotatef(G_fAngle_horizon, 0.0, 1.0, 0.0);
//    glRotatef(G_fAngle_vertical, 1.0, 0.0, 0.0);
//
//    glColor3f(0.0f, 1.0f, 1.0f);
//    glutSolidSphere(0.5, 50, 50);
//
//    glPopMatrix();
    
    
    float max = 0;
    float min = 0;
    float d1, d2;
    float error;
    
    for (int i = 0; i < graph.vexnum; i++) {
        for (int j = i; j < graph.vexnum; j++) {
            if (graph.matrix[i][j] != INFINITE) {
                
                d1 = graph.matrix[i][j];
                d2 = sqrt(pow(node[3*i]-node[3*j], 2) + pow(node[3*i+1]-node[3*j+1], 2) + pow(node[3*i+2]-node[3*j+2], 2));
                error = d1 - d2;
                if (error > max) {
                    max = error;
                }
                if (error < min) {
                    min = error;
                }

            }
        }
    }
    
//    printf("\n\n\n%f\n\n\n%f\n\n\n", max, min);
    
//    if ((max + min) >= 0) {
//        min = -max;
//    } else {
//        max = -min;
//    }

    
    float red, blue;
    error = 0;
    
    for (int i = 0; i < graph.vexnum; i++) {
        for (int j = i; j < graph.vexnum; j++) {
            if (graph.matrix[i][j] != INFINITE) {
                glPushMatrix();

//                glLineWidth(5.0);
//                glBegin(GL_LINES);
//                    glColor3f(1.0, 0.0, 0.0);
//                    glVertex3f(node[3*i], node[3*i+1], node[3*i+2]);
//                    glVertex3f(node[3*j], node[3*j+1], node[3*j+2]);
//                glEnd();
                
//                glColor3f(0.5, 0.0, 0.5);

                error = graph.matrix[i][j] - sqrt(pow(node[3*i]-node[3*j], 2) + pow(node[3*i+1]-node[3*j+1], 2) + pow(node[3*i+2]-node[3*j+2], 2));
                
                if (error >= 0) {
                    red = error / max;
                    blue = 1 - red;
//                    glColor3f(red, 0.5, 0.0);
                } else {
                    red = error / (-min);
                    blue = 1 - red;
//                    glColor3f(0.0, 0.5, blue);
                }
                
//                red = 0.5 + error * 0.5 / max;
//                blue = 1 - red;
                
                glColor3f(red, 0.0, blue);
                
                
                float length = sqrt(pow(node[3*i]-node[3*j], 2) + pow(node[3*i+1]-node[3*j+1], 2) + pow(node[3*i+2]-node[3*j+2], 2));
                
                printf("%f %f\n", length, graph.matrix[i][j]);
                
                if (length < (graph.matrix[i][j]*0.95)) {
                    theTorus = glGenLists(1);  //建立显示列表空间
                    glNewList(theTorus, GL_COMPILE);  //创建一个显示列表
//                    torus(node[3*i]*G_fRate, node[3*i+1]*G_fRate, node[3*i+2]*G_fRate, node[3*j]*G_fRate, node[3*j+1]*G_fRate, node[3*j+2]*G_fRate, graph.matrix[i][j]*G_fRate);
                    torus(node[3*i]*G_fRate, node[3*i+1]*G_fRate, node[3*i+2]*G_fRate, (node[3*i]*G_fRate+node[3*j]*G_fRate)/2, (node[3*i+1]*G_fRate+node[3*j+1]*G_fRate)/2, (node[3*i+2]*G_fRate+node[3*j+2]*G_fRate)/2, graph.matrix[i][j]*G_fRate/2);
                    torus2((node[3*i]*G_fRate+node[3*j]*G_fRate)/2, (node[3*i+1]*G_fRate+node[3*j+1]*G_fRate)/2, (node[3*i+2]*G_fRate+node[3*j+2]*G_fRate)/2, node[3*j]*G_fRate, node[3*j+1]*G_fRate, node[3*j+2]*G_fRate, graph.matrix[i][j]*G_fRate/2);

                    glEndList();  //现实列表完成
                    glCallList(theTorus);

                } else {
                    cyLinder(node[3*i]*G_fRate, node[3*i+1]*G_fRate, node[3*i+2]*G_fRate, node[3*j]*G_fRate, node[3*j+1]*G_fRate, node[3*j+2]*G_fRate, red, blue, graph.matrix[i][j]*G_fRate);
                }
                
//                cyLinder(node[3*i]*G_fRate, node[3*i+1]*G_fRate, node[3*i+2]*G_fRate, node[3*j]*G_fRate, node[3*j+1]*G_fRate, node[3*j+2]*G_fRate, red, blue, graph.matrix[i][j]*G_fRate);
                
                glPopMatrix();
            }
        }
    }
    
    for (int i = 0; i<3*graph.vexnum; i = i+3) {
        
        glPushMatrix();
        
        glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
        glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
        glMaterialfv(GL_FRONT, GL_DIFFUSE, lmodel_ambient);
        glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);

        glTranslatef(node[i]*G_fRate, node[i+1]*G_fRate, node[i+2]*G_fRate);
//        glRotatef(G_fAngle_horizon, 0.0, 1.0, 0.0);
//        glRotatef(G_fAngle_vertical, 1.0, 0.0, 0.0);
        
        glColor3f(0.0f, 1.0f, 0.0f);
        glutSolidSphere(G_fRadius_sphere, 50, 50);
        
        glPopMatrix();
    }
    
    glutSwapBuffers();
    
}

void processSpecialKeys(int key, int x, int y)
{
    switch(key) {
        case GLUT_KEY_LEFT:
            G_fAngle_horizon -= 5.0f;
            break;
        case GLUT_KEY_RIGHT:
            G_fAngle_horizon += 5.0f;
            break;
        case GLUT_KEY_UP:
            G_fAngle_vertical -= 5.0f;
            break;
        case GLUT_KEY_DOWN:
            G_fAngle_vertical += 5.0f;
            break;
    }
    glutPostRedisplay();
}



void processNormalKeys(unsigned char key,int x,int y)
{
    switch(key) {
        case 97:
            G_fDistance -= 1.0f;
            break;
        case 65:
            G_fDistance += 1.0f;
            break;
        case 98:
            G_fRate *= 1.1;
            G_fRadius_sphere *= 1.1;
            G_fRadius_cylinder *= 1.1;
            break;
        case 66:
            G_fRate /= 1.1;
            G_fRadius_sphere /= 1.1;
            G_fRadius_cylinder /= 1.1;
            break;
        case 99:
            
            break;
        case 67:
            
            break;
        case 48:
            G_fDistance = 0.0f;
            break;
        case 27:
            exit(0);
    }
    glutPostRedisplay();
}



void cyLinder(float x0, float y0, float z0, float x1, float y1, float z1, float red, float blue, float g) {
    //如果要在AB两点间画一个圆柱体，其可以
    //先在y轴上画一个同长度的圆柱，然后
    //求出旋转矩阵，将其移至AB
    GLdouble
        dir_x = x1 - x0,
        dir_y = y1 - y0,
        dir_z = z1 - z0;

    GLdouble cy_length = sqrt(dir_x * dir_x + dir_y * dir_y + dir_z * dir_z);   //获得圆柱的长度

    static GLUquadricObj* quad_obj = NULL;

    if (quad_obj == NULL)

        quad_obj = gluNewQuadric();

    gluQuadricDrawStyle(quad_obj, GLU_FILL);

    gluQuadricNormals(quad_obj, GLU_SMOOTH);

    glPushMatrix();

    glTranslated(x0, y0, z0);

    //获得AB的长度
    double length;

    length = sqrt(dir_x * dir_x + dir_y * dir_y + dir_z * dir_z);

    dir_x /= length;
    dir_y /= length;
    dir_z /= length;

    GLdouble
        up_x = 0.0,
        up_y = 1.0,
        up_z = 0.0;

    GLdouble side_x, side_y, side_z;

    //实现向量的叉乘
    side_x = up_y * dir_z - up_z * dir_y;
    side_y = up_z * dir_x - up_x * dir_z;
    side_z = up_x * dir_y - up_y * dir_x;

    length = sqrt(dir_x * dir_x + dir_y * dir_y + dir_z * dir_z);

    side_x /= length;
    side_y /= length;
    side_z /= length;

    up_x = dir_y * side_z - dir_z * side_y;
    up_y = dir_z * side_x - dir_x * side_z;
    up_z = dir_x * side_y - dir_y * side_x;


    //得到变换矩阵
    GLdouble m[] = {
        side_x,side_y,side_z,0.0,
        up_x,up_y,up_z,0.0,
        dir_x,dir_y,dir_z,0.0,
        0.0,0.0,0.0,1.0
    };

    glMultMatrixd(m);           //用m矩阵乘以当前矩阵
    
//    for (int i = 0; i < 16; i++) {
//        cout << m[i] << endl;
//    }
    
//    glColor3f(red, 0.0, blue);
//    glColor3f(0.5, 0.0, 0.5);
    
    GLdouble radius = G_fRadius_cylinder;      //这些参数可以自己设置
    GLdouble slices = 50.0;
    GLdouble stack = 50.0;
    GLfloat a[4];
    
    glColor3f(0.5, 0.0, 0.5);
    a[0] = 0.5;
    a[1] = 0.0;
    a[2] = 0.5;
    a[3] = 1.0;
    glMaterialfv(GL_FRONT, GL_DIFFUSE, a);
    gluCylinder(quad_obj, radius, radius, cy_length, slices, stack);
    
//    glTranslated((g/2)/cy_length*(x1-x0), (g/2)/cy_length*(y1-y0), (g/2)/cy_length*(z1-z0));
    
    glColor3f(red, 0.0, blue);
    a[0] = red;
    a[1] = 0.0;
    a[2] = blue;
    a[3] = 1.0;
    glMaterialfv(GL_FRONT, GL_DIFFUSE, a);
    gluCylinder(quad_obj, radius, radius, cy_length-g/2, slices, stack);
    
//    glTranslated((cy_length-g)/cy_length*(x1-x0), (cy_length-g)/cy_length*(y1-y0), (cy_length-g)/cy_length*(z1-z0));
    
    glColor3f(0.5, 0.0, 0.5);
    a[0] = 0.5;
    a[1] = 0.0;
    a[2] = 0.5;
    a[3] = 1.0;
    glMaterialfv(GL_FRONT, GL_DIFFUSE, a);
    gluCylinder(quad_obj, radius, radius, g/2, slices, stack);
    
    glPopMatrix();
}


static void torus(float x0, float y0, float z0, float x1, float y1, float z1, float desirable)
{
    int numc = 100;
    int numt = 100;
    int i, j, k;
    double s, t, x, y, z, twopi, twopi2, radius;

    //////////////////////////////
    // 计算目标向量
    GLfloat  dx = x1 - x0;
    GLfloat  dy = y1 - y0;
    GLfloat  dz = z1 - z0;

    //获得AB的长度
    double length;

    length = sqrt(dx * dx + dy * dy + dz * dz);
    radius = cal_radius(length, desirable);
//    cout << "length " << length << endl;
//    cout << "radius " << radius << endl;

    float alpha = desirable / radius;

    // 计算平移量
    GLfloat  px = x0;
    GLfloat  py = y0;
    GLfloat  pz = z0;
    // 起始线段的末端点
    GLfloat  bx = px;
    GLfloat  by = -length + py; //
    GLfloat  bz = pz;
    // 计算起始向量
    GLfloat  sx = bx - x0;
    GLfloat  sy = by - y0;
    GLfloat  sz = bz - z0;
    // 计算向量(sx,sy,sz)与向量(dx,dy,dz)的法向量(sy*dz - dy*sz,sz*dx - sx*dz,sx*dy - dx*sy)
    GLfloat fx = sy * dz - dy * sz;
    GLfloat fy = sz * dx - sx * dz;
    GLfloat fz = sx * dy - dx * sy;
    // 求两向量间的夹角
    // 计算第三条边的长度
    GLfloat ax = fabs(x1 - bx);
    GLfloat ay = fabs(y1 - by);
    GLfloat az = fabs(z1 - bz);
    GLfloat len3 = sqrt(ax * ax + ay * ay + az * az);
    // 根据余弦定理计算夹角
    GLfloat angle = acos((length * length * 2 - len3 * len3) / (2 * length * length)) * 180.0f / PI;

    glPushMatrix();
    glTranslatef(x0, y0, z0);
    glRotatef(angle, fx, fy, fz);
//    glRotated(180, 0, 0, 1);
    glTranslatef(0, -length, 0);
    glRotatef(alpha / 2 / PI * 180, 0.0f, 0.0f, -1.0f);
    glTranslatef(-radius, 0, 0);

    twopi = alpha;
    twopi2 = 2 * PI;

    for (i = 0; i < numc; i++) {
        glBegin(GL_QUAD_STRIP);
        for (j = 0; j <= numt-1; j++) {
            for (k = 1; k >= 0; k--) {
                s = (i + k) % numc + 0.5;
                t = j % numt;

                x = (1 + 0.2/radius * cos(s * twopi2 / numc)) * cos(t * twopi / numt);
                y = (1 + 0.2/radius * cos(s * twopi2 / numc)) * sin(t * twopi / numt);//圆环的公式
                z = 0.2/radius * sin(s * twopi2 / numc);
                
//                double sanjiao=sin(s * twopi2 / numc);;
//                x = (1 + 0.01*sanjiao) * cos(t * twopi / numt);
//                y = (1 + 0.01*sanjiao) * sin(t * twopi / numt);//圆环的公式
//                z = 0.01 * sqrt(1-sanjiao*sanjiao);
                x = radius * x;
                y = radius * y;
                z = radius * z;
                glVertex3f(x, y, z);
            }
        }
        glEnd();
        
//        glBegin(GL_QUADS);
//        j = numt;
//            for (k = 1; k >= 0; k--) {
//                s = (i + k) % numc + 0.5;
//                t = j % numt;
//
//                x = (1 + 0.01 * cos(s * twopi2 / numc)) * cos(t * twopi / numt);
//                y = (1 + 0.01 * cos(s * twopi2 / numc)) * sin(t * twopi / numt);//圆环的公式
//                z = 0.01 * sin(s * twopi2 / numc);
//                x = radius * x;
//                y = radius * y;
//                z = radius * z;
//                glVertex3f(x, y, z);
//            }
//        glEnd();
    }

    glPopMatrix();
}

static void torus2(float x0, float y0, float z0, float x1, float y1, float z1, float desirable)
{
    int numc = 100;
    int numt = 100;
    int i, j, k;
    double s, t, x, y, z, twopi, twopi2, radius;

    //////////////////////////////
    // 计算目标向量
    GLfloat  dx = x1 - x0;
    GLfloat  dy = y1 - y0;
    GLfloat  dz = z1 - z0;

    //获得AB的长度
    double length;

    length = sqrt(dx * dx + dy * dy + dz * dz);
    radius = cal_radius(length, desirable);
//    cout << "length " << length << endl;
//    cout << "radius " << radius << endl;

    float alpha = desirable / radius;

    // 计算平移量
    GLfloat  px = x0;
    GLfloat  py = y0;
    GLfloat  pz = z0;
    // 起始线段的末端点
    GLfloat  bx = px;
    GLfloat  by = -length + py; //
    GLfloat  bz = pz;
    // 计算起始向量
    GLfloat  sx = bx - x0;
    GLfloat  sy = by - y0;
    GLfloat  sz = bz - z0;
    // 计算向量(sx,sy,sz)与向量(dx,dy,dz)的法向量(sy*dz - dy*sz,sz*dx - sx*dz,sx*dy - dx*sy)
    GLfloat fx = sy * dz - dy * sz;
    GLfloat fy = sz * dx - sx * dz;
    GLfloat fz = sx * dy - dx * sy;
    // 求两向量间的夹角
    // 计算第三条边的长度
    GLfloat ax = fabs(x1 - bx);
    GLfloat ay = fabs(y1 - by);
    GLfloat az = fabs(z1 - bz);
    GLfloat len3 = sqrt(ax * ax + ay * ay + az * az);
    // 根据余弦定理计算夹角
    GLfloat angle = acos((length * length * 2 - len3 * len3) / (2 * length * length)) * 180.0f / PI;

    glPushMatrix();
    glTranslatef(x0, y0, z0);
    glRotatef(angle, fx, fy, fz);
    glRotated(180, 0, 0, 1);
//    glTranslatef(0, -length, 0);
    glRotatef(alpha / 2 / PI * 180, 0.0f, 0.0f, -1.0f);
    glTranslatef(-radius, 0, 0);

    twopi = alpha;
    twopi2 = 2 * PI;

    for (i = 0; i < numc; i++) {
        glBegin(GL_QUAD_STRIP);
        for (j = 0; j <= numt-1; j++) {
            for (k = 1; k >= 0; k--) {
                s = (i + k) % numc + 0.5;
                t = j % numt;

                x = (1 + 0.2/radius * cos(s * twopi2 / numc)) * cos(t * twopi / numt);
                y = (1 + 0.2/radius * cos(s * twopi2 / numc)) * sin(t * twopi / numt);//圆环的公式
                z = 0.2/radius * sin(s * twopi2 / numc);
                x = radius * x;
                y = radius * y;
                z = radius * z;
                glVertex3f(x, y, z);
            }
        }
        glEnd();
        
//        glBegin(GL_QUADS);
//        j = numt;
//            for (k = 1; k >= 0; k--) {
//                s = (i + k) % numc + 0.5;
//                t = j % numt;
//
//                x = (1 + 0.01 * cos(s * twopi2 / numc)) * cos(t * twopi / numt);
//                y = (1 + 0.01 * cos(s * twopi2 / numc)) * sin(t * twopi / numt);//圆环的公式
//                z = 0.01 * sin(s * twopi2 / numc);
//                x = radius * x;
//                y = radius * y;
//                z = radius * z;
//                glVertex3f(x, y, z);
//            }
//        glEnd();
    }

    glPopMatrix();
}


float cal_radius(float length, float desirable)
{
    //输入弧长desirable，弦长length;
    float low = 0, high = 2 * PI, mid, temp, precision = 0.0000001, radius;

    do {
        mid = (low + high) / 2;
        temp = mid * length - 2 * desirable * sin(mid / 2);//二分法的返回的差值；
        if (temp == 0)
            break;
        else if (temp > 0)
            high = mid;
        else
            low = mid;
    } while (fabs(temp) > precision);
//    cout.setf(ios::fixed);
//    cout.setf(ios::showpoint);
    radius = desirable / mid;
    return radius;
}
