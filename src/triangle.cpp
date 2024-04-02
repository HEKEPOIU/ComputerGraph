
#include "freeglut_std.h"
#include <GL/gl.h>
#include <array>
#include <cmath>
#include <iostream>
#include <ostream>
#include <stdio.h>
/*** freeglut***/
#include <freeglut.h>

void ChangeSize(int, int);
void RenderScene(void);
void MenuCallback(int);
void OnKeyBoardPress(unsigned char, int, int);
void TranslateMatrix(GLfloat, GLfloat, GLfloat);
void RotateMatrix(float angle, float x, float y, float z);
void ScaleMatrix(float x, float y, float z);
void MousePress(int button, int state, int x, int y);

std::array<float, 3> mouseWorldPos1{10, 10, 0};
std::array<float, 3> mouseWorldPos2{10, 10, 0};

bool mousepoint01Status = false;

GLenum glShadeType = GL_SMOOTH;
float r = 0.0;
std::array<float, 3> currentPos{0, 0, 0};
std::array<float, 3> currentRotate{0, 0, 0};
std::array<float, 3> currentScale{1, 1, 1};

int main(int argc, char **argv) {
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
  glutInitWindowSize(800, 800);
  glutInitWindowPosition(600, 80);
  glutCreateWindow("Simple Triangle");
  glutCreateMenu(MenuCallback);
  glutAddMenuEntry("GL_SMOOTH", 1);
  glutAddMenuEntry("GL_FLAT", 2);
  glutAttachMenu(GLUT_RIGHT_BUTTON);
  glutReshapeFunc(ChangeSize);
  glutKeyboardFunc(OnKeyBoardPress);
  glutMouseFunc(MousePress);
  glutDisplayFunc(RenderScene);
  glutMainLoop(); // http://www.programmer-club.com.tw/ShowSameTitleN/opengl/2288.html
  return 0;
}
void ChangeSize(int w, int h) {
  printf("Window Size= %d X %d\n", w, h);
  glViewport(0, 0, w, h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(-10, 10, -10, 10, -50, 50);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
}
void RenderScene(void) {
  glClearColor(0, 0, 0, 1.0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(0, 0, 10.0f, 0, 0, 0, 0, 1, 0);
  glShadeModel(glShadeType);
  glEnable(GL_DEPTH_TEST);

  glBegin(GL_LINES);
  glColor3f(1, 0, 0);
  glVertex3f(mouseWorldPos1[0], mouseWorldPos1[1], mouseWorldPos1[2]);
  glVertex3f(mouseWorldPos2[0], mouseWorldPos2[1], mouseWorldPos2[2]);
  glEnd();

  TranslateMatrix(currentPos[0], currentPos[1], currentPos[2]);
  RotateMatrix(currentRotate[0], mouseWorldPos2[0] - mouseWorldPos1[0],
               mouseWorldPos2[1] - mouseWorldPos1[1],
               0); // x
  RotateMatrix(currentRotate[1], 0, 1, 0);
  RotateMatrix(currentRotate[2], 0, 0, 1);
  ScaleMatrix(currentScale[0], currentScale[1], currentScale[2]);

  glBegin(GL_QUADS);
  glColor3f(1, 0, 0);
  glVertex3f(10, .1, 0);   // RB
  glVertex3f(10, -.1, 0);  // top
  glVertex3f(-10, -.1, 0); // LB
  glVertex3f(-10, .1, 0);  // top
  glEnd();

  glBegin(GL_QUADS);
  glColor3f(0, 1, 0);
  glVertex3f(.1, 10, 0);   // RB
  glVertex3f(-.1, 10, 0);  // top
  glVertex3f(-.1, -10, 0); // LB
  glVertex3f(.1, -10, 0);  // top
  glEnd();

  glBegin(GL_QUADS);
  glColor3f(0, 0, 1);
  glVertex3f(.1, 0, 10);   // RB
  glVertex3f(-.1, 0, 10);  // top
  glVertex3f(-.1, 0, -10); // LB
  glVertex3f(.1, 0, -10);  // top
  glEnd();

  // front
  glBegin(GL_TRIANGLES);
  glColor3f(1, 1, 1);
  glVertex3f(5, -5, 0); // RB

  glColor3f(1, 1, 1);
  glVertex3f(-5, -5, 0); // LB

  glColor3f(1, 1, 1);
  glVertex3f(0, 5, 0); // top
  glEnd();
  // //-------------------------- 1

  glBegin(GL_TRIANGLES);
  glColor3f(0, 1, 0);
  glVertex3f(0, 5, 0);

  glColor3f(0, 1, 0);
  glVertex3f(-5, -5, 0);

  glColor3f(0, 1, 0);
  glVertex3f(-5, -5, -5);
  glEnd();
  // //---------------- 2

  glBegin(GL_TRIANGLES);
  glColor3f(0, 1, 0);
  glVertex3f(0, 5, 0);
  glColor3f(0, 1, 0);
  glVertex3f(0, 5, -5);
  glColor3f(0, 1, 0);
  glVertex3f(-5, -5, -5);
  glEnd();

  // //------------------------ 3
  glBegin(GL_TRIANGLES);
  glColor3f(0, 0, 1);

  glVertex3f(5, -5, 0); // RB

  glColor3f(0, 0, 1);

  glVertex3f(-5, -5, 0); // LB

  glColor3f(0, 0, 1);

  glVertex3f(5, -5, -5);
  glEnd();

  //------------------------ 4
  glBegin(GL_TRIANGLES);
  glColor3f(0, 0, 1);
  glVertex3f(-5, -5, 0); // RB

  glColor3f(0, 0, 1);
  glVertex3f(-5, -5, -5);

  glColor3f(0, 0, 1);
  glVertex3f(5, -5, -5);
  glEnd();
  //------------------------ 5
  glBegin(GL_TRIANGLES);
  glColor3f(1, 0, 1);
  glVertex3f(5, -5, 0); // RB

  glColor3f(1, 0, 1);
  glVertex3f(0, 5, 0); // top

  glColor3f(1, 0, 1);
  glVertex3f(5, -5, -5);
  glEnd();

  //------------------------ 6
  glBegin(GL_TRIANGLES);

  glColor3f(1, 0, 1);
  glVertex3f(0, 5, 0); // top

  glColor3f(1, 0, 1);
  glVertex3f(5, -5, -5);
  glColor3f(1, 0, 1);
  glVertex3f(0, 5, -5);
  glEnd();

  glBegin(GL_TRIANGLES);

  glColor3f(0, 1, 0);
  glVertex3f(5, -5, -5);

  glColor3f(0, 1, 0);
  glVertex3f(-5, -5, -5);

  glColor3f(1, 0, 0);
  glVertex3f(0, 5, -5);

  glEnd();
  glutSwapBuffers();
}

void MenuCallback(int value) {
  switch (value) {
  case 1:
    glShadeType = GL_SMOOTH;
    break;
  case 2:
    glShadeType = GL_FLAT;
    break;
  }
  glutPostRedisplay();
}

void print4mat(GLfloat array[16]) {
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      printf("%-3f ", array[4 * j + i]);
    }
    printf("\n");
  }
}

void TranslateMatrix(GLfloat x, GLfloat y, GLfloat z) {
  printf("translate %f %f %f\n", x, y, z);
  GLfloat rotMatrix[16];
  glGetFloatv(GL_MODELVIEW_MATRIX, rotMatrix);

  rotMatrix[12] = x;
  rotMatrix[13] = y;
  rotMatrix[14] = z;
  print4mat(rotMatrix);
  glMultMatrixf(rotMatrix);
}

void RotateMatrix(float angle, float x, float y, float z) {
  printf("rotate %f %f %f %f\n", angle, x, y, z);
  GLfloat rotMatrix[16];
  float angleDeg2Rad = angle * (M_PI / 180.0f);

  float cosAngle = cos(angleDeg2Rad);
  float sinAngle = sin(angleDeg2Rad);

  float len = sqrt(x * x + y * y + z * z);

  if (len != 0) {
    x /= len;
    y /= len;
    z /= len;
  }

  float omc = 1.0f - cosAngle;

  rotMatrix[0] = x * x * omc + cosAngle;
  rotMatrix[1] = y * x * omc + z * sinAngle;
  rotMatrix[2] = x * z * omc - y * sinAngle;
  rotMatrix[3] = 0.0f;

  rotMatrix[4] = x * y * omc - z * sinAngle;
  rotMatrix[5] = y * y * omc + cosAngle;
  rotMatrix[6] = y * z * omc + x * sinAngle;
  rotMatrix[7] = 0.0f;

  rotMatrix[8] = x * z * omc + y * sinAngle;
  rotMatrix[9] = y * z * omc - x * sinAngle;
  rotMatrix[10] = z * z * omc + cosAngle;
  rotMatrix[11] = 0.0f;

  rotMatrix[12] = 0.0f;
  rotMatrix[13] = 0.0f;
  rotMatrix[14] = 0.0f;
  rotMatrix[15] = 1.0f;
  print4mat(rotMatrix);
  glMultMatrixf(rotMatrix);
}

void ScaleMatrix(float x, float y, float z) {
  printf("scale %f %f %f\n", x, y, z);
  GLfloat scaleMatrix[16]{};
  scaleMatrix[0] = x;
  scaleMatrix[5] = y;
  scaleMatrix[10] = z;
  scaleMatrix[15] = 1.0f;
  print4mat(scaleMatrix);
  glMultMatrixf(scaleMatrix);
}

void OnKeyBoardPress(unsigned char key, int x, int y) {
  std::cout << key << std::endl;

  switch (key) {
  case 'r':
    currentPos[0] = 0;
    currentPos[1] = 0;
    currentPos[2] = 0;
    currentRotate[0] = 0;
    currentRotate[1] = 0;
    currentRotate[2] = 0;
    break;
  case 'd':
    currentPos[0] += 1;
    break;
  case 'a':
    currentPos[0] -= 1;
    break;
  case 'j':
    currentPos[1] += 1;
    break;
  case 'l':
    currentPos[1] -= 1;
    break;
  case 'i':
    currentPos[2] += 1;
    break;
  case 'p':
    currentPos[2] -= 1;
    break;
  case 'x':
    currentRotate[0] += 1;
    break;
  case 'y':
    currentRotate[1] += 1;
    break;
  case 'z':
    currentRotate[2] += 1;
    break;
  case '1':
    currentScale[0] += 0.1;
    break;
  case '2':
    currentScale[1] += 0.1;
    break;
  case '3':
    currentScale[2] += 0.1;
    break;
  case '4':
    currentScale[0] -= 0.1;
    break;
  case '5':
    currentScale[1] -= 0.1;
    break;
  case '6':
    currentScale[2] -= 0.1;
    break;
  }
  glutPostRedisplay();
}
void MousePress(int button, int state, int x, int y) {
  if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
    if (mousepoint01Status == false) {
      mouseWorldPos1[0] = (2 * ((float)x / 800) - 1) * (10);
      mouseWorldPos1[1] = (-2 * ((float)y / 800) + 1) * (10);
      mouseWorldPos1[2] = 0;
      mousepoint01Status = true;
    } else {
      mouseWorldPos2[0] = (2 * ((float)x / 800) - 1) * (10);
      mouseWorldPos2[1] = (-2 * ((float)y / 800) + 1) * (10);
      mouseWorldPos2[2] = 0;
      mousepoint01Status = false;
    }

    // gluInvertMatrix();
    glutPostRedisplay();
  }
}
