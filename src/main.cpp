
#include "freeglut_std.h"
#include <GL/gl.h>
#include <array>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <memory>
#include <ostream>
#include <stdio.h>

/*** freeglut***/
#include <freeglut.h>
#include <vector>

#include "OBJLoader.hpp"

void ChangeSize(int, int);
void RenderScene(void);
void MenuCallback(int);
void OnKeyBoardPress(unsigned char, int, int);
void TranslateMatrix(GLfloat, GLfloat, GLfloat);
void RotateMatrix(float angle, float x, float y, float z);
void ScaleMatrix(float x, float y, float z);
void MousePress(int button, int state, int x, int y);
void LoadObjFileAndRecreateMenu(
    std::vector<std::shared_ptr<DrawableObject>> &obj_list, int &menu_id);

std::unique_ptr<OBJLoader> obj_loader = std::make_unique<OBJLoader>();
std::vector<std::shared_ptr<DrawableObject>> objs{};
std::shared_ptr<DrawableObject> current_display_obj{};
int menu_id = -1;

std::array<int, 2> windowSize{800, 800};
int main(int argc, char **argv) {

  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
  glutInitWindowSize(800, 800);
  glutInitWindowPosition(600, 80);
  glutCreateWindow("Simple Triangle");

  LoadObjFileAndRecreateMenu(objs, menu_id);

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
  glEnable(GL_DEPTH_TEST);

  RotateMatrix(45, 1, 1, 0);

  if (current_display_obj != nullptr) {
    current_display_obj->draw(GL_QUADS, {1, 1, 1});
  }

  glutSwapBuffers();
}

void MenuCallback(int value) {
  if (value > objs.size()) {
    LoadObjFileAndRecreateMenu(objs, menu_id);
    return;
  }
  current_display_obj = objs[value - 1];

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

void OnKeyBoardPress(unsigned char key, int x, int y) {
  std::cout << key << std::endl;
  glutPostRedisplay();
}
void MousePress(int button, int state, int x, int y) {

  // gluInvertMatrix();
  glutPostRedisplay();
}

void LoadObjFileAndRecreateMenu(
    std::vector<std::shared_ptr<DrawableObject>> &obj_list, int &menu_id) {

  if (menu_id != -1) {
    glutDestroyMenu(menu_id);
  }

  menu_id = glutCreateMenu(MenuCallback); // GLUT Not Support remove menu
                                          // item, so we need reCreate one.

  obj_list.clear(); // clear old date, make sure no repeated data.
  int count = 1;
  for (const auto &entry : std::filesystem::directory_iterator(RESOURCE_DIR)) {
    if (std::filesystem::is_regular_file(entry)) {
      std::cout << "file_name: " << entry.path().filename() << std::endl;
      obj_list.push_back(OBJLoader().get_OBJ(entry.path().string()));
      glutAddMenuEntry(entry.path().filename().c_str(), count);
      count++;
    }
  }
  glutAddMenuEntry("Reload Path", count);
  glutAttachMenu(GLUT_RIGHT_BUTTON);
}
