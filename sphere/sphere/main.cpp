#define GL_SILENCE_DEPRECATION

#include <stdio.h>
#include <GL/freeglut.h>

#define WINDOW_HEIGHT 300
#define WINDOW_WIDTH 300
#define INFINITE 65535
#define MAXSIZE 100

typedef struct graph
{
    float *vexs;
    int vexnum;
    float matrix[MAXSIZE][MAXSIZE];
}Graph;

GLfloat G_fDistance = 0.0f;
GLfloat G_fAngle_horizon = 0.0f;
GLfloat G_fAngle_vertical = 0.0f;

void myinit(void);
void myReshape(GLsizei w, GLsizei h);
void display(void);
void processSpecialKeys(int key, int x, int y);
void processNormalKeys(unsigned char key,int x,int y);

float G_vLit0Position[4] = { 5.0f, 5.0f, 5.0f, 1.0f };
float G_vLit0Ambient[4] = { 0.2f, 0.2f, 0.2f, 1.0f };
float G_vLit0Diffuse[4] = { 0.8f, 0.8f, 0.8f, 1.0f };
float G_vLit0Specular[4] = { 0.5f, 0.5f, 0.5f, 1.0f };
//float G_vMaterialSpecu[4] = { 0.0f, 0.0f, 0.0f, 1.0f };
//
//float G_vLit1Position[4] = { 0.0f, 0.0f, 0.0f, 1.0f };

int main(int argc, char* argv[])
{
    glutInit(&argc, argv);
    
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);

    glutInitWindowSize(500, 500);
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
    glViewport(0, 0, w, h);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(120, 1.0*(GLfloat)w/(GLfloat)h, 1.0, 300.0);
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
    
    fp1 = fopen("/Users/hanzhekaiscomputer/Desktop/N1-100.txt", "r");    // node.txt
    
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
    
    fp2 = fopen("/Users/hanzhekaiscomputer/Desktop/P1-.txt", "r");  // path.txt
    
    for (int i =0; i < graph.vexnum; i++) {
        for (int j = 0; j < graph.vexnum; j++) {
            graph.matrix[i][j] = INFINITE;
        }
    }
    
    while (fscanf(fp2, "%d %d %f", &first, &second, &weight) != EOF) {
        graph.matrix[first-1][second-1] = weight;
        graph.matrix[second-1][first-1] = weight;
    }

//    for (int i = 0; i < graph.vexnum; i++) {
//        for (int j = 0; j < graph.vexnum; j++) {
//            printf("%.2f ", graph.matrix[i][j]);
//        }
//        printf("\n");
//    }
    
    fclose(fp2);
    
    
    
    
    
    glClearColor(0.0f,0.0f,0.0f,0.0f);
    glClearDepth(1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    glDisable(GL_LIGHTING);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glEnable(GL_LIGHTING);
    
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glLightfv(GL_LIGHT0, GL_POSITION, G_vLit0Position);
    glLightfv(GL_LIGHT0, GL_AMBIENT, G_vLit0Ambient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, G_vLit0Diffuse);
    glLightfv(GL_LIGHT0, GL_SPECULAR, G_vLit0Specular);
    
    glMatrixMode(GL_MODELVIEW);
    
    glLoadIdentity();
    glTranslatef(0.0, 0.0, -G_fDistance);
    glRotatef(G_fAngle_horizon, 0.0, 1.0, 0.0);
    glRotatef(G_fAngle_vertical, 1.0, 0.0, 0.0);
    
//    gluLookAt(0, 0, 50, 0, 0, -1, 0, 1, 0);
//    glScalef(0.5, 0.5, 0.5);
    
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
    
    

    
    for (int i = 0; i < graph.vexnum; i++) {
        for (int j = i; j < graph.vexnum; j++) {
            if (graph.matrix[i][j] != INFINITE) {
                glPushMatrix();

                glColor3f(1.0, 0.0, 0.0);
                glLineWidth(5.0);
                glBegin(GL_LINES);
                    glVertex3f(node[3*i], node[3*i+1], node[3*i+2]);
                    glVertex3f(node[3*j], node[3*j+1], node[3*j+2]);
                glEnd();

                glPopMatrix();
            }
        }
    }
    
    for (int i = 0; i<3*graph.vexnum; i = i+3) {
        
        glPushMatrix();
        
        glTranslatef(node[i], node[i+1], node[i+2]);
//        glRotatef(G_fAngle_horizon, 0.0, 1.0, 0.0);
//        glRotatef(G_fAngle_vertical, 1.0, 0.0, 0.0);
        
        glColor3f(0.0f, 1.0f, 0.0f);
        glutSolidSphere(1, 50, 50);
        
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
        case 48:
            G_fDistance = 0.0f;
            break;
        case 27:
            exit(0);
    }
    glutPostRedisplay();
}
