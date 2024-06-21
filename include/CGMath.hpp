#ifndef CHMATH_HPP
#define CHMATH_HPP

struct Vec3 {
  float x;
  float y;
  float z;

  Vec3(float x = 0, float y = 0, float z = 0) : x(x), y(y), z(z) {}

  Vec3 operator+(const Vec3 &a) const { return Vec3(a.x + x); }
};

#endif
