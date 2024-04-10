#ifndef DrawableObject_hpp
#define DrawableObject_hpp

#include <GL/gl.h>
#include <array>
#include <vector>
class DrawableObject {

public:
  DrawableObject(const std::vector<std::array<float, 3>> &vertex,
                 const std::vector<std::vector<int>> &face)
      : _vertex(vertex), _face(face) {}

  void draw(GLenum mode, std::array<float, 3> glColor3f);

private:
  std::vector<std::array<float, 3>> _vertex;

  // the face have change more than 3 vertex.
  std::vector<std::vector<int>> _face;

  // We can read more data form the file, but now we don't need.
};

#endif