#ifndef DRAWABLEOBJECT_HPP
#define DRAWABLEOBJECT_HPP

#include <array>

#include "CGMath.hpp"
#include "Texture.hpp"
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

  DrawableObject(const std::vector<Vec3> &vertices,
                 const std::vector<Vec2> &tex_coords,
                 const std::vector<Vec3> &normals,
                 const std::vector<std::vector<int>> &vertex_indices,
                 const std::vector<std::vector<int>> &tex_coord_indices,
                 const std::vector<std::vector<int>> &normal_indices)
      : _vertex(vertices), _tex_coords(tex_coords), _normals(normals),
        _vertex_indices(vertex_indices), _tex_coord_indices(tex_coord_indices),
        _normal_indices(normal_indices) {}

  void draw(DrawType mode);
  int get_face_count() { return _vertex_indices.size(); }

  float get_proper_scale() { return properScale; }

  // this function will adject the transform to the target postion (include
  // position, rotation, scale)
  void set_transform_to_target(const Vec3 &target_position,
                               const std::array<int, 6> &view_space);

  std::array<double, 6> get_bbox();

  std::shared_ptr<Transform> &get_transform() { return _transform; }
  void SetTexture(const std::shared_ptr<Texture> tex);

private:
  std::vector<Vec3> _vertex;
  std::vector<Vec2> _tex_coords;
  std::vector<Vec3> _normals;
  // the face have change more than 3 vertex.
  std::vector<std::vector<int>> _vertex_indices;
  std::vector<std::vector<int>> _tex_coord_indices;
  std::vector<std::vector<int>> _normal_indices;
  std::shared_ptr<Texture> _texture;

  std::shared_ptr<Transform> _transform = std::make_shared<Transform>();
  float properScale{1.0f};

  // We can read more data form the file, but now we don't need.
};

#endif
