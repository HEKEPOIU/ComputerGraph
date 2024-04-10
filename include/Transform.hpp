#ifndef TRANSFORM_HPP
#define TRANSFORM_HPP

#include <GL/gl.h>
#include <array>
class Transform {

public:
  enum class RotationType { EULER, ANYAXIS };

  Transform() { reset_transform(); };
  void set_location(const std::array<float, 3> &location);
  void set_rotation_by_euler(const std::array<float, 3> &rotation);

  // note this function will break the ```_rotation``` variable, because I don't
  // have know how to convert rotation matrix to euler.
  void set_rotation_by_axis(float angle, const std::array<float, 3> &axis);
  void set_scale(const std::array<float, 3> &scale);
  const std::array<float, 3> &get_location() { return _location; };
  const std::array<float, 3> &get_rotation() { return _rotation; };
  const std::array<float, 3> &get_scale() { return _scale; };

  void modify_location(const std::array<float, 3> &axis) {
    std::array<float, 3> new_location = get_location();
    new_location[0] += axis[0];
    new_location[1] += axis[1];
    new_location[2] += axis[2];
    set_location(new_location);
  }

  void modify_rotation_by_axis(float angle, const std::array<float, 3> &axis) {
    float new_angle = _current_rotation_angle;
    new_angle += angle;
    set_rotation_by_axis(new_angle, axis);
  }
  void modify_rotation(const std::array<float, 3> &axis) {
    std::array<float, 3> new_rotation = get_rotation();
    new_rotation[0] += axis[0];
    new_rotation[1] += axis[1];
    new_rotation[2] += axis[2];

    set_rotation_by_euler(new_rotation);
  }

  void modify_scale(const std::array<float, 3> &axis) {
    std::array<float, 3> new_scale = get_scale();
    new_scale[0] += axis[0];
    new_scale[1] += axis[1];
    new_scale[2] += axis[2];
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

  std::array<float, 3> _location{0, 0, 0};
  std::array<float, 3> _rotation{0, 0, 0};
  std::array<float, 3> _current_rotation_axis{0, 0, 1};
  float _current_rotation_angle = 0;
  std::array<float, 3> _scale{1.0, 1.0, 1.0};
};

#endif