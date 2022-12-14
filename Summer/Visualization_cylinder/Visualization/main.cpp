#define GL_SILENCE_DEPRECATION

#include <stdio.h>
#include<iostream>
#include <GL/freeglut.h>
#include <cmath>

#define WINDOW_WIDTH 800
#define WINDOW_HEIGHT 800
#define INFINITE 65535
#define MAXSIZE 1000

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
void cyLinder(float x0, float y0, float z0, float x1, float y1, float z1);

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

//    glutIdleFunc(display);

    glutMainLoop();

    return 0;
    
}

void myinit(void)
{
    glEnable(GL_DEPTH_TEST);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    glEnable(GL_LIGHT0);
//    glEnable(GL_LIGHT1);
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
    float node[3*MAXSIZE];             // ?????????????????????????????????
    graph.vexs = node;                 // graph?????????????????????node??????
    graph.vexnum = 0;                  // ??????????????????
    
    float x, y, z;                      // ????????????
    int first, second;                  // ??????????????????
    float weight;                       // ????????????????????????
    
    FILE *fp1, *fp2 = NULL;
    
    fp1 = fopen("/Users/hanzhekaiscomputer/Desktop/P30-50-output.txt", "r");    // node.txt
    
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
    
    fp2 = fopen("/Users/hanzhekaiscomputer/Desktop/P30-50.txt", "r");  // path.txt
    
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
                
                cyLinder(node[3*i]*G_fRate, node[3*i+1]*G_fRate, node[3*i+2]*G_fRate, node[3*j]*G_fRate, node[3*j+1]*G_fRate, node[3*j+2]*G_fRate);
                
                glPopMatrix();
            }
        }
    }
    
    for (int i = 0; i<3*graph.vexnum; i = i+3) {
        
        glPushMatrix();
        
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



void cyLinder(float x0, float y0, float z0, float x1, float y1, float z1) {
    //????????????AB???????????????????????????????????????
    //??????y??????????????????????????????????????????
    //?????????????????????????????????AB
    GLdouble
        dir_x = x1 - x0,
        dir_y = y1 - y0,
        dir_z = z1 - z0;

    GLdouble cy_length = sqrt(dir_x * dir_x + dir_y * dir_y + dir_z * dir_z);   //?????????????????????

    static GLUquadricObj* quad_obj = NULL;

    if (quad_obj == NULL)

        quad_obj = gluNewQuadric();

    gluQuadricDrawStyle(quad_obj, GLU_FILL);

    gluQuadricNormals(quad_obj, GLU_SMOOTH);

    glPushMatrix();

    glTranslated(x0, y0, z0);

    //??????AB?????????
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

    //?????????????????????
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


    //??????????????????
    GLdouble m[] = {
        side_x,side_y,side_z,0.0,
        up_x,up_y,up_z,0.0,
        dir_x,dir_y,dir_z,0.0,
        0.0,0.0,0.0,1.0
    };

    glMultMatrixd(m);           //???m????????????????????????
    
//    for (int i = 0; i < 16; i++) {
//        cout << m[i] << endl;
//    }
    
//    glColor3f(red, 0.0, blue);
//    glColor3f(0.5, 0.0, 0.5);
    
    GLdouble radius = G_fRadius_cylinder;      //??????????????????????????????
    GLdouble slices = 50.0;
    GLdouble stack = 50.0;

    gluCylinder(quad_obj, radius, radius, cy_length, slices, stack);

    glPopMatrix();
}
