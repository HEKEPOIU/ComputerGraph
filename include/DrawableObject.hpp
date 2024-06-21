#ifndef DRAWABLEOBJECT_HPP
#define DRAWABLEOBJECT_HPP

#include <array>

#include "Transform.hpp"
#include <memory>
#include <vector>

class DrawableObject {

public:
  enum class DrawType {
    POINTS,
    LINES,
    FACES,
  };

  DrawableObject(const std::vector<Vec3> &vertex,
                 const std::vector<std::vector<int>> &face)
      : _vertex(vertex), _faces(face) {
    _faceColor.resize(face.size());
    for (int i = 0; i < _faceColor.size(); i++) {
      _faceColor[i] = {0.5f, 0.5f, 0.5f};
    }
  }

  void draw(DrawType mode);
  int get_face_count() { return _faces.size(); }
  void set_face_color(const std::vector<Vec3> &faceColor) {
    _faceColor = faceColor;
  }

  float get_proper_scale() { return properScale; }

  // this function will adject the transform to the target postion (include
  // position, rotation, scale)
  void set_transform_to_target(const Vec3 &target_position,
                               const std::array<int, 6> &view_space);

  std::array<float, 6> get_bbox();

  std::shared_ptr<Transform> &get_transform() { return _transform; }

private:
  std::vector<Vec3> _vertex;

  // the face have change more than 3 vertex.
  std::vector<std::vector<int>> _faces;
  std::vector<Vec3> _faceColor;

  std::shared_ptr<Transform> _transform = std::make_shared<Transform>();
  float properScale{1.0f};

  // We can read more data form the file, but now we don't need.
};

#endif
