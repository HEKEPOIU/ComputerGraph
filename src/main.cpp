#include "DrawableObject.hpp"
#include <array>
#include <memory>
#include <stdio.h>

/*** freeglut***/
#include <freeglut.h>
#include <freeglut_std.h>
#include <vector>

void ChangeSize(int, int);
void RenderScene(void);
void MenuCallback(int);
void ColorModeMenuCallback(int);
void OnKeyBoardPress(unsigned char, int, int);
void MousePress(int button, int state, int x, int y);
void RenderModeMenuCallback(int);
void RotationModeCallback(int);

std::array<int, 6> ortho_settings = {-10, 10, -10, 10, -50, 50};
std::array<double, 3> camera_pos = {0, 0, 10};
std::array<double, 3> camera_look_at = {0, 0, 0};

std::array<int, 2> windowSize{800, 800};
int main(int argc, char **argv) {

  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
  glutInitWindowSize(800, 800);
  glutInitWindowPosition(600, 80);
  glutCreateWindow("Simple Triangle");

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
  glOrtho(ortho_settings[0], ortho_settings[1], ortho_settings[2],
          ortho_settings[3], ortho_settings[4], ortho_settings[5]);
  glMatrixMode(GL_MODELVIEW);
  gluLookAt(camera_pos[0], camera_pos[1], camera_pos[2], camera_look_at[0],
            camera_look_at[1], camera_look_at[2], 0, 1, 0);
  glLoadIdentity();
}
void RenderScene(void) {
  glClearColor(0, 0, 0, 1.0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  // TODO: we should calculate the up vector, but now I am lazy to do it.
  gluLookAt(camera_pos[0], camera_pos[1], camera_pos[2], camera_look_at[0],
            camera_look_at[1], camera_look_at[2], 0, 1, 0);
  glEnable(GL_DEPTH_TEST);

  glutSwapBuffers();
}

void MenuCallback(int value) { glutPostRedisplay(); }

void ColorModeMenuCallback(int value) { glutPostRedisplay(); }

void RenderModeMenuCallback(int value) { glutPostRedisplay(); }

void RotationModeCallback(int value) { glutPostRedisplay(); }

void OnKeyBoardPress(unsigned char key, int x, int y) { glutPostRedisplay(); }

void MousePress(int button, int state, int x, int y) { glutPostRedisplay(); }

void LoadObjFileAndRecreateMenu(
    std::vector<std::shared_ptr<DrawableObject>> &obj_list, int &menu_id) {}
