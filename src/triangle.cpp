
#include "freeglut_std.h"
#include <GL/gl.h>
#include <array>
#include <cmath>
#include <iostream>
#include <ostream>
#include <stdio.h>
/*** freeglut***/
#include <freeglut.h>
#include <vector>

void ChangeSize(int, int);
void RenderScene(void);
void MenuCallback(int);
void OnKeyBoardPress(unsigned char, int, int);
void TranslateMatrix(GLfloat, GLfloat, GLfloat);
void RotateMatrix(float angle, float x, float y, float z);
void ScaleMatrix(float x, float y, float z);
void MousePress(int button, int state, int x, int y);
void ChangeGridPointSize(std::vector<std::vector<int>> &grid_point, int x,
                         int y);
void DrawOutLine(int x, int y);
void DrawGrid(std::vector<std::vector<int>> &grid_point);

std::array<int, 2> windowSize{800, 800};
std::vector<std::vector<int>> _gridPoint{};
std::vector<bool> _gridDraw{};
int size_x = 21; // -10 ~ 10
int size_y = 21;

int main(int argc, char **argv) {
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
  glutInitWindowSize(800, 800);
  glutInitWindowPosition(600, 80);
  glutCreateWindow("Simple Triangle");
  glutCreateMenu(MenuCallback);
  glutAddMenuEntry("10", 1);
  glutAddMenuEntry("15", 2);
  glutAddMenuEntry("20", 3);
  ChangeGridPointSize(_gridPoint, 21, 21);
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
  glOrtho(-400, 400, -400, 400, -50, 50);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
}
void RenderScene(void) {
  glClearColor(0, 0, 0, 1.0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(0, 0, 10.0f, 0, 0, 0, 0, 1, 0);
  glEnable(GL_DEPTH_TEST);
  glViewport(0, 0, windowSize[0], windowSize[1]);

  DrawGrid(_gridPoint);

  glutSwapBuffers();
}

void MenuCallback(int value) {
  // TODO: change grid size.
  switch (value) {
  case 1:
    ChangeGridPointSize(_gridPoint, 21, 21);
    break;
  case 2:
    ChangeGridPointSize(_gridPoint, 31, 31);
    break;
  case 3:
    ChangeGridPointSize(_gridPoint, 41, 41);
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
  GLfloat rotMatrix[16];
  glGetFloatv(GL_MODELVIEW_MATRIX, rotMatrix);

  rotMatrix[12] = x;
  rotMatrix[13] = y;
  rotMatrix[14] = z;
  glMultMatrixf(rotMatrix);
}

void RotateMatrix(float angle, float x, float y, float z) {
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
  glMultMatrixf(rotMatrix);
}

void ScaleMatrix(float x, float y, float z) {
  GLfloat scaleMatrix[16]{};
  scaleMatrix[0] = x;
  scaleMatrix[5] = y;
  scaleMatrix[10] = z;
  scaleMatrix[15] = 1.0f;
  glMultMatrixf(scaleMatrix);
}

void ChangeGridPointSize(std::vector<std::vector<int>> &grid_point, int totalX,
                         int totalY) {
  size_x = totalX;
  size_y = totalY;
  grid_point.clear();
  _gridDraw.clear();
  int halfX = totalX / 2;
  int halfY = totalY / 2;
  for (int x = -halfX; x <= halfX; x++) {
    for (int y = -halfY; y <= halfY; y++) {
      grid_point.push_back({x, y});
      _gridDraw.push_back(false);
    }
  }
}

void DrawOutLine(int x, int y) {
  int w = windowSize[0];
  int h = windowSize[1];
  int sizeX = w / size_x;
  int sizeY = h / size_y;
  glColor3f(1.0f, 1.0f, 1.0f);

  glBegin(GL_LINES);

  glVertex3f(x - sizeX / 2, y + sizeY / 2, 0);
  glVertex3f(x - sizeX / 2, y - sizeY / 2, 0);

  glVertex3f(x - sizeX / 2, y - sizeY / 2, 0);
  glVertex3f(x + sizeX / 2, y - sizeY / 2, 0);

  glVertex3f(x + sizeX / 2, y - sizeY / 2, 0);
  glVertex3f(x + sizeX / 2, y + sizeY / 2, 0);

  glVertex3f(x + sizeX / 2, y + sizeY / 2, 0);
  glVertex3f(x - sizeX / 2, y + sizeY / 2, 0);

  glEnd();
}

void DrawFillCell(int x, int y) {
  int w = windowSize[0];
  int h = windowSize[1];
  int sizeX = w / size_x;
  int sizeY = h / size_y;

  glColor3f(1.0f, 1.0f, 1.0f);
  glBegin(GL_QUADS);

  glVertex3f(x - sizeX / 2, y + sizeY / 2, 0);
  glVertex3f(x - sizeX / 2, y - sizeY / 2, 0);

  glVertex3f(x + sizeX / 2, y - sizeY / 2, 0);

  glVertex3f(x + sizeX / 2, y + sizeY / 2, 0);

  glEnd();
}

void DrawGrid(std::vector<std::vector<int>> &grid_point) {
  int spaceingX = windowSize[0] / (size_x);
  int spaceingY = windowSize[1] / (size_y);
  for (int i = 0; i < grid_point.size(); i++) {
    if (_gridDraw[i] == true) {
      DrawFillCell(grid_point[i][0] * spaceingX, grid_point[i][1] * spaceingY);
    } else {
      DrawOutLine(grid_point[i][0] * spaceingX, grid_point[i][1] * spaceingY);
    }
  }
}

void OnKeyBoardPress(unsigned char key, int x, int y) { glutPostRedisplay(); }
void ScrenPosToGridPoint(float x, float y, int &gridX, int &gridY) {
  int w = windowSize[0];
  int h = windowSize[1];
  int sizeX = w / (size_x + 1);
  int sizeY = h / (size_y + 1);
  gridX = x / sizeX;
  gridY = y / sizeY;
  std::cout << gridX << " " << gridY << std::endl;
}
void MousePress(int button, int state, int x, int y) {
  if (state == GLUT_DOWN) {
    float mousePosX = (2 * ((float)x / 800) - 1) * (400);
    float mousePosY = (-2 * ((float)y / 800) + 1) * (400);
    int gridX = 0;
    int gridY = 0;
    ScrenPosToGridPoint(mousePosX, mousePosY, gridX, gridY);
    gridX += size_x / 2;
    gridY += size_y / 2;
    _gridDraw[gridY + gridX * size_y] = !_gridDraw[gridY + gridX * size_y];

    glutPostRedisplay();
  }
  // gluInvertMatrix();
}
