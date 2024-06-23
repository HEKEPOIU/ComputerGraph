#ifndef CHMATH_HPP
#define CHMATH_HPP

struct Vec3 {
  double x;
  double y;
  double z;

  Vec3(double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z) {}

  Vec3 operator+(const Vec3 &a) const {
    return Vec3(a.x + x, a.y + y, a.z + z);
  }

  Vec3 operator-(const Vec3 &a) const {
    return Vec3(a.x - x, a.y - y, a.z - z);
  }
};

struct Vec2 {
  double x;
  double y;

  Vec2(double x = 0, double y = 0) : x(x), y(y) {}

  Vec2 operator+(const Vec3 &a) const { return Vec2(a.x + x, a.y + y); }

  Vec2 operator-(const Vec2 &a) const { return Vec2(a.x - x, a.y - y); }
};

#endif
