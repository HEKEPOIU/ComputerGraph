#include "DrawableObject.hpp"
#include "Transform.hpp"
#include <array>
#include <filesystem>
#include <iostream>
#include <memory>
#include <ostream>
#include <stdio.h>

/*** freeglut***/
#include <freeglut.h>
#include <freeglut_std.h>
#include <vector>

#include "OBJLoader.hpp"

enum class ControlMode {
  CAMERAMOVE,
  LOOKATPOINTMOVE,
  TRANSLATE,
  ROTATE,
  SCALE,
};

enum class CurrentControlAxis {
  X,
  Y,
  Z,
};
ControlMode current_control_mode = ControlMode::CAMERAMOVE;

void ChangeSize(int, int);
void RenderScene(void);
void MenuCallback(int);
void ColorModeMenuCallback(int);
void OnKeyBoardPress(unsigned char, int, int);
void MousePress(int button, int state, int x, int y);
void RenderModeMenuCallback(int);
void RotationModeCallback(int);

std::array<float, 3> mouseViewport1WorldPos1{10, 10, 0};
std::array<float, 3> mouseViewport1WorldPos2{10, 10, 0};
std::array<int, 6> ortho_settings = {-10, 10, -10, 10, -50, 50};
std::array<double, 3> camera_pos = {0, 0, 10};
std::array<double, 3> camera_look_at = {0, 0, 0};

bool mousepoint01Status = false;

std::array<float, 3> generate_random_color();

void LoadObjFileAndRecreateMenu(
    std::vector<std::shared_ptr<DrawableObject>> &obj_list, int &menu_id);

std::unique_ptr<OBJLoader> obj_loader = std::make_unique<OBJLoader>();
std::vector<std::shared_ptr<DrawableObject>> objs{};
std::shared_ptr<DrawableObject> current_display_obj{};
int current_menu_id = -1;
DrawableObject::DrawType current_draw_mode = DrawableObject::DrawType::FACES;
CurrentControlAxis current_control_axis = CurrentControlAxis::X;

std::array<int, 2> windowSize{800, 800};
int main(int argc, char **argv) {

  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
  glutInitWindowSize(800, 800);
  glutInitWindowPosition(600, 80);
  glutCreateWindow("Simple Triangle");

  LoadObjFileAndRecreateMenu(objs, current_menu_id);

  for (auto &obj : objs) {
    obj->set_transform_to_target({0, 0, 0}, ortho_settings);
  }

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

  glBegin(GL_LINES);
  glColor3f(1.0f, 0.0f, 0.0f);
  glVertex3f(100.0f, 0.0f, 0.0f);
  glVertex3f(-100.0f, 0.0f, 0.0f);
  glEnd();

  glBegin(GL_LINES);
  glColor3f(0.0f, 1.0f, 0.0f);
  glVertex3f(0.0f, 100.0f, 0.0f);
  glVertex3f(0.0f, -100.0f, 0.0f);
  glEnd();

  glBegin(GL_LINES);
  glColor3f(0.0f, 0.0f, 1.0f);
  glVertex3f(0.0f, 0.0f, 100.0f);
  glVertex3f(0.0f, 0.0f, -100.0f);
  glEnd();

  if (current_display_obj != nullptr) {
    if (current_display_obj->get_transform()->get_rotation_type() ==
        Transform::RotationType::ANYAXIS) {
      glBegin(GL_LINES);
      glColor3f(1.0f, 1.0f, 1.0f);
      glVertex3fv(mouseViewport1WorldPos1.data());
      glVertex3fv(mouseViewport1WorldPos2.data());
      glEnd();
    }

    current_display_obj->draw(current_draw_mode);
  }

  glutSwapBuffers();
}

void MenuCallback(int value) {
  if (value > objs.size()) {
    LoadObjFileAndRecreateMenu(objs, current_menu_id);
    for (auto &obj : objs) {
      obj->set_transform_to_target({0, 0, 0}, ortho_settings);
    }
    return;
  }
  current_display_obj = objs[value - 1];

  glutPostRedisplay();
}

void ColorModeMenuCallback(int value) {

  for (int i = 0; i < objs.size(); i++) {
    std::vector<std::array<float, 3>> face_colors;
    face_colors.resize(objs[i]->get_face_count());

    for (auto &face_color : face_colors) {
      std::array<float, 3> color;
      if (value == 1) {
        color = generate_random_color();
      } else {
        color = {0.5f, 0.5f, 0.5f};
      }
      face_color = color;
    }

    objs[i]->set_face_color(face_colors);
  }

  glutPostRedisplay();
}

void RenderModeMenuCallback(int value) {
  if (value == 1) {
    current_draw_mode = DrawableObject::DrawType::POINTS;
  } else if (value == 2) {
    current_draw_mode = DrawableObject::DrawType::LINES;
  } else {
    current_draw_mode = DrawableObject::DrawType::FACES;
  }
  glutPostRedisplay();
}

void RotationModeCallback(int value) {
  auto transform = current_display_obj->get_transform();
  switch (value) {
  case 1: {
    transform->set_rotation_type(Transform::RotationType::EULER);
    break;
  }
  case 2: {
    transform->set_rotation_type(Transform::RotationType::ANYAXIS);
    break;
  }
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

// fuck the switch hell, it's so ugly, i know how to fix it, but I'm too lazy,
// some one help me to separate to self Input class.
void OnKeyBoardPress(unsigned char key, int x, int y) {
  switch (key) {
  case 'q': {
    current_control_mode = ControlMode::CAMERAMOVE;
    std::cout << "select mode, do nothing Now" << std::endl;
    break;
  }
  case 'w': {
    current_control_mode = ControlMode::TRANSLATE;
    std::cout << "Translate mode" << std::endl;
    break;
  }
  case 'e': {
    current_control_mode = ControlMode::ROTATE;
    std::cout << "Rotate mode" << std::endl;
    break;
  }
  case 'r': {
    current_control_mode = ControlMode::SCALE;
    std::cout << "Scale mode" << std::endl;
    break;
  }
  case 'f': {
    current_control_mode = ControlMode::LOOKATPOINTMOVE;
    std::cout << "Change LookAt Point" << std::endl;
    break;
  }
  }

  switch (key) {
  case 'x': {
    current_control_axis = CurrentControlAxis::X;
    break;
  }
  case 'y': {
    current_control_axis = CurrentControlAxis::Y;
    break;
  }
  case 'z': {
    current_control_axis = CurrentControlAxis::Z;
    break;
  }
  }

  if (!(key == 'a' || key == 'd' || key == ' ')) {
    return;
  }

  float value_to_add = key == 'a' ? -1.0f : 1.0f;

  if (current_display_obj == nullptr) {
    return;
  }

  auto transform = current_display_obj->get_transform();

  if (key == ' ') {
    current_display_obj->set_transform_to_target({0, 0, 0}, ortho_settings);
    transform->set_rotation_by_euler({0, 0, 0});
    transform->set_rotation_by_axis(0, {0, 0, 1});
    camera_pos = {0, 0, 10};
    camera_look_at = {0, 0, 0};
    glutPostRedisplay();
    return;
  }
  switch (current_control_mode) {
  case ControlMode::CAMERAMOVE: {
    switch (current_control_axis) {
    case CurrentControlAxis::X: {
      std::cout << "x" << std::endl;
      camera_pos[0] += value_to_add;
      break;
    }
    case CurrentControlAxis::Y: {
      camera_pos[1] += value_to_add;

      break;
    }
    case CurrentControlAxis::Z: {
      camera_pos[2] += value_to_add;
      break;
    }
    }
    break;
  }
  case ControlMode::LOOKATPOINTMOVE: {
    switch (current_control_axis) {
    case CurrentControlAxis::X: {
      std::cout << "x" << std::endl;
      camera_look_at[0] += value_to_add;
      break;
    }
    case CurrentControlAxis::Y: {
      camera_look_at[1] += value_to_add;
      break;
    }
    case CurrentControlAxis::Z: {
      camera_look_at[2] += value_to_add;
      break;
    }
    }
    break;
  }
  case ControlMode::TRANSLATE: {
    switch (current_control_axis) {
    case CurrentControlAxis::X: {
      transform->modify_location({value_to_add, 0, 0});
      break;
    }
    case CurrentControlAxis::Y: {
      transform->modify_location({0, value_to_add, 0});
      break;
    }
    case CurrentControlAxis::Z: {
      transform->modify_location({0, 0, value_to_add});
      break;
    }
    }
    break;
  }
  case ControlMode::ROTATE: {
    switch (transform->get_rotation_type()) {

    case Transform::RotationType::EULER: {
      switch (current_control_axis) {
      case CurrentControlAxis::X: {
        transform->modify_rotation({value_to_add, 0, 0});
        break;
      }
      case CurrentControlAxis::Y: {
        transform->modify_rotation({0, value_to_add, 0});
        break;
      }
      case CurrentControlAxis::Z: {
        transform->modify_rotation({0, 0, value_to_add});
        break;
      }
      }
      break;
    }
    case Transform::RotationType::ANYAXIS: {
      transform->modify_rotation_by_axis(
          value_to_add,
          {mouseViewport1WorldPos2[0] - mouseViewport1WorldPos1[0],
           mouseViewport1WorldPos2[1] - mouseViewport1WorldPos1[1],
           mouseViewport1WorldPos2[2] - mouseViewport1WorldPos1[2]});
      break;
    }
    }
    break;
  }
  case ControlMode::SCALE: {
    value_to_add *= 0.1f * current_display_obj->get_proper_scale();
    switch (current_control_axis) {
    case CurrentControlAxis::X: {
      transform->modify_scale({value_to_add, 0, 0});
      break;
    }
    case CurrentControlAxis::Y: {
      transform->modify_scale({0, value_to_add, 0});
      break;
    }
    case CurrentControlAxis::Z: {
      transform->modify_scale({0, 0, value_to_add});
      break;
    }
    }
    break;
  }
  }
  glutPostRedisplay();
}

void MousePress(int button, int state, int x, int y) {
  if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
    if (mousepoint01Status == false) {
      mouseViewport1WorldPos1[0] = (2 * ((float)x / 800) - 1) * (10);
      mouseViewport1WorldPos1[1] = (-2 * ((float)y / 800) + 1) * (10);
      mouseViewport1WorldPos1[2] = 0;
      mousepoint01Status = true;
    } else {
      mouseViewport1WorldPos2[0] = (2 * ((float)x / 800) - 1) * (10);
      mouseViewport1WorldPos2[1] = (-2 * ((float)y / 800) + 1) * (10);
      mouseViewport1WorldPos2[2] = 0;
      mousepoint01Status = false;
    }
  }
  // gluInvertMatrix();
  glutPostRedisplay();
}

void LoadObjFileAndRecreateMenu(
    std::vector<std::shared_ptr<DrawableObject>> &obj_list, int &menu_id) {

  if (menu_id != -1) {
    glutDestroyMenu(menu_id);
  }
  obj_list.clear(); // clear old date, make sure no repeated data.

  int submenu_color = glutCreateMenu(ColorModeMenuCallback);
  glutAddMenuEntry("RandomColor", 1);
  glutAddMenuEntry("SingleColor", 2);

  int submenu_renderMode = glutCreateMenu(RenderModeMenuCallback);
  glutAddMenuEntry("Points", 1);
  glutAddMenuEntry("Lines", 2);
  glutAddMenuEntry("Faces", 3);

  int submenu_rotation_mode = glutCreateMenu(RotationModeCallback);
  glutAddMenuEntry("Euler", 1);
  glutAddMenuEntry("Any Axis", 2);

  menu_id = glutCreateMenu(MenuCallback); // GLUT Not Support remove menu
                                          // item, so we need reCreate one.

  glutAddSubMenu("ColorMode", submenu_color);
  glutAddSubMenu("RenderMode", submenu_renderMode);
  glutAddSubMenu("RotationMode", submenu_rotation_mode);

  int count = 1;
  for (const auto &entry : std::filesystem::directory_iterator(RESOURCE_DIR)) {
    if (std::filesystem::is_regular_file(entry)) {
      std::cout << "file_name: " << entry.path().filename() << std::endl;
      const std::string &filePath = entry.path().string();
      obj_list.push_back(obj_loader->get_OBJ(filePath));
      std::string filename = entry.path().filename().string();
      const char *menuEntry = filename.c_str();
      glutAddMenuEntry(menuEntry, count);
      count++;
    }
  }
  glutAddMenuEntry("Reload Path", count);

  glutAttachMenu(GLUT_RIGHT_BUTTON);
}

std::array<float, 3> generate_random_color() {
  std::array<float, 3> color{};
  color[0] = (float)rand() / RAND_MAX;
  color[1] = (float)rand() / RAND_MAX;
  color[2] = (float)rand() / RAND_MAX;
  return color;
}