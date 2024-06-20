#include "Transform.hpp"
#include <cmath>
#include <freeglut.h>
#include <freeglut_std.h>


std::array<GLfloat, 16> Transform::get_rotate_matrix(float angle, float x,
                                                     float y, float z) {

  std::array<GLfloat, 16> rotMatrix;
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
  return rotMatrix;
}

std::array<GLfloat, 16> Transform::get_scale_matrix(float x, float y, float z) {
  std::array<GLfloat, 16> scaleMatrix = {x,    0.0f, 0.0f, 0.0f, 0.0f, y,
                                         0.0f, 0.0f, 0.0f, 0.0f, z,    0.0f,
                                         0.0f, 0.0f, 0.0f, 1.0f};
  return scaleMatrix;
}

std::array<GLfloat, 16> Transform::get_translate_matrix(GLfloat x, GLfloat y,
                                                        GLfloat z) {
  std::array<GLfloat, 16> location_matrix;
  glGetFloatv(GL_MODELVIEW_MATRIX, location_matrix.data());

  location_matrix[12] = x;
  location_matrix[13] = y;
  location_matrix[14] = z;
  return location_matrix;
}

void Transform::set_location(const std::array<float, 3> &location) {
  _location = location;
}

void Transform::set_rotation_by_euler(const std::array<float, 3> &rotation) {
  _rotation = rotation;
}

void Transform::set_rotation_by_axis(float angle,
                                     const std::array<float, 3> &axis) {
  _current_rotation_angle = angle;
  _current_rotation_axis = axis;
}

void Transform::set_scale(const std::array<float, 3> &scale) { _scale = scale; }

void Transform::Apply_transform() {

  glMultMatrixf(
      get_translate_matrix(_location[0], _location[1], _location[2]).data());
  if (_rotation_type == RotationType::EULER) {
    glMultMatrixf(get_rotate_matrix(_rotation[0], 1, 0, 0).data());
    glMultMatrixf(get_rotate_matrix(_rotation[1], 0, 1, 0).data());
    glMultMatrixf(get_rotate_matrix(_rotation[2], 0, 0, 1).data());

  } else if (_rotation_type == RotationType::ANYAXIS) {
    glMultMatrixf(
        get_rotate_matrix(_current_rotation_angle, _current_rotation_axis[0],
                          _current_rotation_axis[1], _current_rotation_axis[2])
            .data());
  }
  glMultMatrixf(get_scale_matrix(_scale[0], _scale[1], _scale[2]).data());
}

void Transform::reverse_transform() {
  glMultMatrixf(
      get_translate_matrix(-_location[0], -_location[1], -_location[2]).data());
  if (_rotation_type == RotationType::EULER) {
    glMultMatrixf(get_rotate_matrix(-_rotation[0], 1, 0, 0).data());
    glMultMatrixf(get_rotate_matrix(-_rotation[1], 0, 1, 0).data());
    glMultMatrixf(get_rotate_matrix(-_rotation[2], 0, 0, 1).data());

  } else if (_rotation_type == RotationType::ANYAXIS) {
    glMultMatrixf(
        get_rotate_matrix(-_current_rotation_angle, _current_rotation_axis[0],
                          _current_rotation_axis[1], _current_rotation_axis[2])
            .data());
  }
  glMultMatrixf(get_scale_matrix(-_scale[0], -_scale[1], -_scale[2]).data());
}