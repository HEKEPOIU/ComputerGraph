#include "DrawableObject.hpp"
#include <GL/gl.h>
#include <iostream>
#include <random>
#include <stdlib.h>

void DrawableObject::draw(GLenum mode, std::array<float, 3> glColor3f) {
  int loop_back = 0;
  if (mode == GL_LINES) {
    loop_back = 1;
  }

  glBegin(mode);
  for (const auto &face : _face) {
    glColor3fv(glColor3f.data());
    for (int i = 0; i < face.size() + loop_back; i++) {
      int j = i % face.size();
      glVertex3fv(&_vertex[face[j] - 1][0]);
    }
  }
  glEnd();
}