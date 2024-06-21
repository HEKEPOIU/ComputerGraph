#include "DrawableObject.hpp"
#include <algorithm>
#include <array>
#include <freeglut.h>
#include <freeglut_std.h>

void DrawableObject::draw(DrawType mode) {

  GLenum gl_mode;
  int loop_back = 0;
  if (mode == DrawType::LINES) {
    loop_back = 1;
    gl_mode = GL_LINES;
  } else if (mode == DrawType::POINTS) {
    gl_mode = GL_POINTS;
  } else {
    gl_mode = _faces[0].size() == 3 ? GL_TRIANGLES : GL_POLYGON;
  }

  _transform->Apply_transform();
  glPointSize(5.0f);
  for (int i = 0; i < _faces.size(); i++) {

    glBegin(gl_mode);
    glColor3f(_faceColor[i].x, _faceColor[i].y, _faceColor[i].z);

    for (int k = 0; k < _faces[i].size() + loop_back; k++) {
      int j = k % _faces[i].size();
      glVertex3fv(&_vertex[_faces[i][j] - 1].x);
    }
    glEnd();
  }
  glPointSize(1.0f);

  _transform->reverse_transform();
}

std::array<float, 6> DrawableObject::get_bbox() {
  float min_x{_vertex[0].x};
  float max_x{_vertex[0].x};
  float min_y{_vertex[0].y};
  float max_y{_vertex[0].y};
  float min_z{_vertex[0].z};
  float max_z{_vertex[0].z};

  for (auto &vertex : _vertex) {
    if (vertex.x < min_x) {
      min_x = vertex.x;
    }
    if (vertex.x > max_x) {
      max_x = vertex.x;
    }
    if (vertex.y < min_y) {
      min_y = vertex.y;
    }
    if (vertex.y > max_y) {
      max_y = vertex.y;
    }
    if (vertex.z < min_z) {
      min_z = vertex.z;
    }
    if (vertex.z > max_z) {
      max_z = vertex.z;
    }
  }

  return {min_x, max_x, min_y, max_y, min_z, max_z};
}

void DrawableObject::set_transform_to_target(
    const Vec3 &target_position,
    const std::array<int, 6> &view_space) {
  auto bbox = get_bbox();

  std::array<int, 3> view_space_size{view_space[1] - view_space[0],
                                     view_space[3] - view_space[2],
                                     view_space[5] - view_space[4]};

  float scale_x = view_space_size[0] / (bbox[1] - bbox[0]);
  float scale_y = view_space_size[1] / (bbox[3] - bbox[2]);
  float scale_z = view_space_size[2] / (bbox[5] - bbox[4]);
  properScale = std::min(scale_x, std::min(scale_y, scale_z)) * 0.8f;
  _transform->set_scale({properScale, properScale, properScale});
  _transform->set_location(target_position);

  // because the scale Will effect the position, so we need multiply it.
  _transform->modify_location({-((bbox[1] + bbox[0]) * properScale / 2),
                               -((bbox[3] + bbox[2]) * properScale / 2),
                               -((bbox[5] + bbox[4]) * properScale / 2)});
}
