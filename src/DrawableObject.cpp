#include "DrawableObject.hpp"
#include <GL/gl.h>

void DrawableObject::draw(DrawType mode) {

  GLenum gl_mode;
  int loop_back = 0;
  if (mode == DrawType::LINES) {
    loop_back = 1;
    gl_mode = GL_LINES;
  } else if (mode == DrawType::POINTS) {
    gl_mode = GL_POINTS;
  } else {
    gl_mode = _face[0].size() == 3 ? GL_TRIANGLES : GL_POLYGON;
  }

  _transform->Apply_transform();
  glPointSize(2.0f);
  glBegin(gl_mode);
  for (int i = 0; i < _face.size(); i++) {

    glColor3f(_faceColor[i][0], _faceColor[i][1], _faceColor[i][2]);

    for (int k = 0; k < _face[i].size() + loop_back; k++) {
      int j = k % _face[i].size();
      glVertex3fv(&_vertex[_face[i][j] - 1][0]);
    }
  }
  glPointSize(1.0f);
  glEnd();

  _transform->reverse_transform();
}