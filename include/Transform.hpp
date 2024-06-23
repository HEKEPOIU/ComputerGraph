#ifndef TRANSFORM_HPP
#define TRANSFORM_HPP

#include "CGMath.hpp"
#include <array>
#include <freeglut.h>
#include <freeglut_std.h>

#define M_PI 3.14159265358979323846

class Transform {

public:
  enum class RotationType { EULER, ANYAXIS };

  Transform() { reset_transform(); };
  void set_location(const Vec3 &location);
  void set_rotation_by_euler(const Vec3 &rotation);

  // note this function will break the ```_rotation``` variable, because I don't
  // have know how to convert rotation matrix to euler.
  void set_rotation_by_axis(float angle, const Vec3 &axis);
  void set_scale(const Vec3 &scale);
  const Vec3 &get_location() { return _location; };
  const Vec3 &get_rotation() { return _rotation; };
  const Vec3 &get_scale() { return _scale; };

  void modify_location(const Vec3 &pos) {
    Vec3 new_location = get_location();
    new_location = new_location + pos;
    set_location(new_location);
  }

  void modify_rotation_by_axis(float angle, const Vec3 &axis) {
    float new_angle = _current_rotation_angle;
    new_angle += angle;
    set_rotation_by_axis(new_angle, axis);
  }
  void modify_rotation(const Vec3 &axis) {
    Vec3 new_rotation = get_rotation();
    new_rotation = new_rotation + axis;
    set_rotation_by_euler(new_rotation);
  }

  void modify_scale(const Vec3 &scale) {
    Vec3 new_scale = get_scale();
    new_scale = new_scale + scale;
    set_scale(new_scale);
  }

  void Apply_transform();
  void reverse_transform();

  void reset_transform() {
    _location = {0, 0, 0};
    _rotation = {0, 0, 0};
    _scale = {1.0, 1.0, 1.0};
    _current_rotation_axis = {0, 0, 1};
    _current_rotation_angle = 0;
  }

  void set_rotation_type(RotationType type) { _rotation_type = type; }
  RotationType get_rotation_type() { return _rotation_type; }

private:
  RotationType _rotation_type = RotationType::EULER;

  std::array<GLfloat, 16> get_rotate_matrix(float angle, float x, float y,
                                            float z);
  std::array<GLfloat, 16> get_scale_matrix(float x, float y, float z);
  std::array<GLfloat, 16> get_translate_matrix(GLfloat x, GLfloat y, GLfloat z);

  Vec3 _location{0, 0, 0};
  Vec3 _rotation{0, 0, 0};
  Vec3 _current_rotation_axis{0, 0, 1};
  float _current_rotation_angle = 0;
  Vec3 _scale{1.0, 1.0, 1.0};
};

#endif
