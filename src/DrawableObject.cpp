#include "DrawableObject.hpp"
#include "Texture.hpp"
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
    gl_mode = _vertex_indices[0].size() == 3 ? GL_TRIANGLES : GL_POLYGON;
  }

  if (_texture) {
    _texture->ApplyTexture();
  }
  _transform->Apply_transform();
  for (int i = 0; i < _vertex_indices.size(); i++) {
    glBegin(gl_mode);

    for (int k = 0; k < _vertex_indices[i].size() + loop_back; k++) {
      int j = k % _vertex_indices[i].size();
      glTexCoord2dv(&_tex_coords[_tex_coord_indices[i][j]].x);
      glNormal3dv(&_normals[_normal_indices[i][j]].x);
      glVertex3dv(&_vertex[_vertex_indices[i][j]].x);
    }
    glEnd();
  }
  _transform->reverse_transform();
}

std::array<double, 6> DrawableObject::get_bbox() {
  double min_x{_vertex[0].x};
  double max_x{_vertex[0].x};
  double min_y{_vertex[0].y};
  double max_y{_vertex[0].y};
  double min_z{_vertex[0].z};
  double max_z{_vertex[0].z};

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
    const Vec3 &target_position, const std::array<int, 6> &view_space) {
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

void DrawableObject::SetTexture(const std::shared_ptr<Texture> tex) {
  _texture = tex;
}