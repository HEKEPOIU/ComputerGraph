#include "freeglut_std.h"
#include <GL/gl.h>
#include <X11/Xlib.h>
#include <algorithm>
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include <cstring>
#include <iostream>
#include <math.h>
#define PI 3.14159265359f

typedef float M3DVector2f[2];  // 3D points = 3D Vectors, but we need a
typedef double M3DVector2d[2]; // 2D representations sometimes... (x,y) order

typedef float M3DVector3f[3];  // Vector of three floats (x, y, z)
typedef double M3DVector3d[3]; // Vector of three doubles (x, y, z)

typedef float M3DVector4f[4]; // Lesser used... Do we really need these?
typedef double
    M3DVector4d[4]; // Yes, occasionaly we do need a trailing w component

// 3x3 matrix - column major. X vector is 0, 1, 2, etc.
//		0	3	6
//		1	4	7
//		2	5	8
typedef float
    M3DMatrix33f[9]; // A 3 x 3 matrix, column major (floats) - OpenGL Style
typedef double
    M3DMatrix33d[9]; // A 3 x 3 matrix, column major (doubles) - OpenGL Style

// 4x4 matrix - column major. X vector is 0, 1, 2, etc.
//	0	4	8	12
//	1	5	9	13
//	2	6	10	14
//	3	7	11	15
typedef float
    M3DMatrix44f[16]; // A 4 X 4 matrix, column major (floats) - OpenGL style
typedef double
    M3DMatrix44d[16]; // A 4 x 4 matrix, column major (doubles) - OpenGL style

// Transformation matrix to project shadow
inline void m3dScaleVector3(M3DVector3f v, const float scale) {
  v[0] *= scale;
  v[1] *= scale;
  v[2] *= scale;
}
inline void m3dCrossProduct3(M3DVector3f result, const M3DVector3f u,
                             const M3DVector3f v) {
  result[0] = u[1] * v[2] - v[1] * u[2];
  result[1] = -u[0] * v[2] + v[0] * u[2];
  result[2] = u[0] * v[1] - v[0] * u[1];
}
inline float m3dGetVectorLengthSquared3(const M3DVector3f u) {
  return (u[0] * u[0]) + (u[1] * u[1]) + (u[2] * u[2]);
}
inline float m3dGetVectorLength3(const M3DVector3f u) {
  return sqrtf(m3dGetVectorLengthSquared3(u));
}
inline void m3dNormalizeVector3(M3DVector3f u) {
  m3dScaleVector3(u, 1.0f / m3dGetVectorLength3(u));
}

void m3dGetPlaneEquation(M3DVector4f planeEq, const M3DVector3f p1,
                         const M3DVector3f p2, const M3DVector3f p3) {
  // Get two vectors... do the cross product
  M3DVector3f v1, v2;

  // V1 = p3 - p1
  v1[0] = p3[0] - p1[0];
  v1[1] = p3[1] - p1[1];
  v1[2] = p3[2] - p1[2];

  // V2 = P2 - p1
  v2[0] = p2[0] - p1[0];
  v2[1] = p2[1] - p1[1];
  v2[2] = p2[2] - p1[2];

  // Unit normal to plane - Not sure which is the best way here
  m3dCrossProduct3(planeEq, v1, v2);
  m3dNormalizeVector3(planeEq);

  // Back substitute to get D
  planeEq[3] = -(planeEq[0] * p3[0] + planeEq[1] * p3[1] + planeEq[2] * p3[2]);
}

void m3dMakePlanarShadowMatrix(M3DMatrix44f proj, const M3DVector4f planeEq,
                               const M3DVector3f vLightPos) {
  // These just make the code below easier to read. They will be
  // removed by the optimizer.
  float a = planeEq[0];
  float b = planeEq[1];
  float c = planeEq[2];
  float d = planeEq[3];

  float dx = -vLightPos[0];
  float dy = -vLightPos[1];
  float dz = -vLightPos[2];

  // Now build the projection matrix
  proj[0] = b * dy + c * dz;
  proj[1] = -a * dy;
  proj[2] = -a * dz;
  proj[3] = 0.0;

  proj[4] = -b * dx;
  proj[5] = a * dx + c * dz;
  proj[6] = -b * dz;
  proj[7] = 0.0;

  proj[8] = -c * dx;
  proj[9] = -c * dy;
  proj[10] = a * dx + b * dy;
  proj[11] = 0.0;

  proj[12] = -d * dx;
  proj[13] = -d * dy;
  proj[14] = -d * dz;
  proj[15] = a * dx + b * dy + c * dz;
  // Shadow matrix ready
}
void m3dFindNormal(M3DVector3f result, const M3DVector3f point1,
                   const M3DVector3f point2, const M3DVector3f point3) {
  M3DVector3f v1, v2; // Temporary vectors

  // Calculate two vectors from the three points. Assumes counter clockwise
  // winding!
  v1[0] = point1[0] - point2[0];
  v1[1] = point1[1] - point2[1];
  v1[2] = point1[2] - point2[2];

  v2[0] = point2[0] - point3[0];
  v2[1] = point2[1] - point3[1];
  v2[2] = point2[2] - point3[2];

  // Take the cross product of the two vectors to get
  // the normal vector.
  m3dCrossProduct3(result, v1, v2);
}
inline void m3dSetMatrixColumn44(M3DMatrix44f dst, const M3DVector4f src,
                                 const int column) {
  memcpy(dst + (4 * column), src, sizeof(float) * 4);
}

inline void m3dSetMatrixColumn44(M3DMatrix44d dst, const M3DVector4d src,
                                 const int column) {
  memcpy(dst + (4 * column), src, sizeof(double) * 4);
}
inline void m3dCopyVector3(M3DVector3f dst, const M3DVector3f src) {
  memcpy(dst, src, sizeof(M3DVector3f));
}
inline void m3dCopyVector3(M3DVector3d dst, const M3DVector3d src) {
  memcpy(dst, src, sizeof(M3DVector3d));
}
#define A(row, col) a[(col << 2) + row]
#define B(row, col) b[(col << 2) + row]
#define P(row, col) product[(col << 2) + row]
void m3dMatrixMultiply44(M3DMatrix44f product, const M3DMatrix44f a,
                         const M3DMatrix44f b) {
  for (int i = 0; i < 4; i++) {
    float ai0 = A(i, 0), ai1 = A(i, 1), ai2 = A(i, 2), ai3 = A(i, 3);
    P(i, 0) = ai0 * B(0, 0) + ai1 * B(1, 0) + ai2 * B(2, 0) + ai3 * B(3, 0);
    P(i, 1) = ai0 * B(0, 1) + ai1 * B(1, 1) + ai2 * B(2, 1) + ai3 * B(3, 1);
    P(i, 2) = ai0 * B(0, 2) + ai1 * B(1, 2) + ai2 * B(2, 2) + ai3 * B(3, 2);
    P(i, 3) = ai0 * B(0, 3) + ai1 * B(1, 3) + ai2 * B(2, 3) + ai3 * B(3, 3);
  }
}

// Ditto above, but for doubles
void m3dMatrixMultiply44(M3DMatrix44d product, const M3DMatrix44d a,
                         const M3DMatrix44d b) {
  for (int i = 0; i < 4; i++) {
    double ai0 = A(i, 0), ai1 = A(i, 1), ai2 = A(i, 2), ai3 = A(i, 3);
    P(i, 0) = ai0 * B(0, 0) + ai1 * B(1, 0) + ai2 * B(2, 0) + ai3 * B(3, 0);
    P(i, 1) = ai0 * B(0, 1) + ai1 * B(1, 1) + ai2 * B(2, 1) + ai3 * B(3, 1);
    P(i, 2) = ai0 * B(0, 2) + ai1 * B(1, 2) + ai2 * B(2, 2) + ai3 * B(3, 2);
    P(i, 3) = ai0 * B(0, 3) + ai1 * B(1, 3) + ai2 * B(2, 3) + ai3 * B(3, 3);
  }
}
void m3dLoadIdentity44(M3DMatrix44f m) {
  // Don't be fooled, this is still column major
  static M3DMatrix44f identity = {1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f,
                                  0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f,
                                  0.0f, 0.0f, 0.0f, 1.0f};

  memcpy(m, identity, sizeof(M3DMatrix44f));
}
inline void m3dCrossProduct(M3DVector3f result, const M3DVector3f u,
                            const M3DVector3f v) {
  result[0] = u[1] * v[2] - v[1] * u[2];
  result[1] = -u[0] * v[2] + v[0] * u[2];
  result[2] = u[0] * v[1] - v[0] * u[1];
}

inline void m3dCrossProduct(M3DVector3d result, const M3DVector3d u,
                            const M3DVector3d v) {
  result[0] = u[1] * v[2] - v[1] * u[2];
  result[1] = -u[0] * v[2] + v[0] * u[2];
  result[2] = u[0] * v[1] - v[0] * u[1];
}

// 4x4 double
void m3dLoadIdentity44(M3DMatrix44d m) {
  static M3DMatrix44d identity = {1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                                  0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0};

  memcpy(m, identity, sizeof(M3DMatrix44d));
}
void m3dLoadIdentity33(M3DMatrix33f m) {
  // Don't be fooled, this is still column major
  static M3DMatrix33f identity = {1.0f, 0.0f, 0.0f, 0.0f, 1.0f,
                                  0.0f, 0.0f, 0.0f, 1.0f};

  memcpy(m, identity, sizeof(M3DMatrix33f));
}

// 3x3 double
void m3dLoadIdentity33(M3DMatrix33d m) {
  // Don't be fooled, this is still column major
  static M3DMatrix33d identity = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

  memcpy(m, identity, sizeof(M3DMatrix33d));
}
#define M33(row, col) m[col * 3 + row]
void m3dRotationMatrix33(M3DMatrix33f m, float angle, float x, float y,
                         float z) {

  float mag, s, c;
  float xx, yy, zz, xy, yz, zx, xs, ys, zs, one_c;

  s = float(sin(angle));
  c = float(cos(angle));

  mag = float(sqrt(x * x + y * y + z * z));

  // Identity matrix
  if (mag == 0.0f) {
    m3dLoadIdentity33(m);
    return;
  }

  // Rotation matrix is normalized
  x /= mag;
  y /= mag;
  z /= mag;

  xx = x * x;
  yy = y * y;
  zz = z * z;
  xy = x * y;
  yz = y * z;
  zx = z * x;
  xs = x * s;
  ys = y * s;
  zs = z * s;
  one_c = 1.0f - c;

  M33(0, 0) = (one_c * xx) + c;
  M33(0, 1) = (one_c * xy) - zs;
  M33(0, 2) = (one_c * zx) + ys;

  M33(1, 0) = (one_c * xy) + zs;
  M33(1, 1) = (one_c * yy) + c;
  M33(1, 2) = (one_c * yz) - xs;

  M33(2, 0) = (one_c * zx) - ys;
  M33(2, 1) = (one_c * yz) + xs;
  M33(2, 2) = (one_c * zz) + c;
}

void m3dRotationMatrix44(M3DMatrix44f m, float angle, float x, float y,
                         float z) {
  float mag, s, c;
  float xx, yy, zz, xy, yz, zx, xs, ys, zs, one_c;

  s = float(sin(angle));
  c = float(cos(angle));

  mag = float(sqrt(x * x + y * y + z * z));

  // Identity matrix
  if (mag == 0.0f) {
    m3dLoadIdentity44(m);
    return;
  }

  // Rotation matrix is normalized
  x /= mag;
  y /= mag;
  z /= mag;

#define M(row, col) m[col * 4 + row]

  xx = x * x;
  yy = y * y;
  zz = z * z;
  xy = x * y;
  yz = y * z;
  zx = z * x;
  xs = x * s;
  ys = y * s;
  zs = z * s;
  one_c = 1.0f - c;

  M(0, 0) = (one_c * xx) + c;
  M(0, 1) = (one_c * xy) - zs;
  M(0, 2) = (one_c * zx) + ys;
  M(0, 3) = 0.0f;

  M(1, 0) = (one_c * xy) + zs;
  M(1, 1) = (one_c * yy) + c;
  M(1, 2) = (one_c * yz) - xs;
  M(1, 3) = 0.0f;

  M(2, 0) = (one_c * zx) - ys;
  M(2, 1) = (one_c * yz) + xs;
  M(2, 2) = (one_c * zz) + c;
  M(2, 3) = 0.0f;

  M(3, 0) = 0.0f;
  M(3, 1) = 0.0f;
  M(3, 2) = 0.0f;
  M(3, 3) = 1.0f;

#undef M
}
inline void m3dLoadVector3(M3DVector3f v, const float x, const float y,
                           const float z) {
  v[0] = x;
  v[1] = y;
  v[2] = z;
}
inline void m3dLoadVector3(M3DVector3d v, const double x, const double y,
                           const double z) {
  v[0] = x;
  v[1] = y;
  v[2] = z;
}
// Create a Translation matrix. Only 4x4 matrices have translation components
inline void m3dTranslationMatrix44(M3DMatrix44f m, float x, float y, float z) {
  m3dLoadIdentity44(m);
  m[12] = x;
  m[13] = y;
  m[14] = z;
}
__inline void m3dRotateVector(M3DVector3f vOut, const M3DVector3f p,
                              const M3DMatrix33f m) {
  vOut[0] = m[0] * p[0] + m[3] * p[1] + m[6] * p[2];
  vOut[1] = m[1] * p[0] + m[4] * p[1] + m[7] * p[2];
  vOut[2] = m[2] * p[0] + m[5] * p[1] + m[8] * p[2];
}

// Ditto above, but for doubles
__inline void m3dRotateVector(M3DVector3d vOut, const M3DVector3d p,
                              const M3DMatrix33d m) {
  vOut[0] = m[0] * p[0] + m[3] * p[1] + m[6] * p[2];
  vOut[1] = m[1] * p[0] + m[4] * p[1] + m[7] * p[2];
  vOut[2] = m[2] * p[0] + m[5] * p[1] + m[8] * p[2];
}

inline void m3dTranslationMatrix44(M3DMatrix44d m, double x, double y,
                                   double z) {
  m3dLoadIdentity44(m);
  m[12] = x;
  m[13] = y;
  m[14] = z;
}
class GLFrame {
protected:
  M3DVector3f vOrigin;  // Where am I?
  M3DVector3f vForward; // Where am I going?
  M3DVector3f vUp;      // Which way is up?

public:
  // Default position and orientation. At the origin, looking
  // down the positive Z axis (right handed coordinate system).
  GLFrame(void) {
    // At origin
    vOrigin[0] = 0.0f;
    vOrigin[1] = 0.0f;
    vOrigin[2] = 0.0f;

    // Up is up (+Y)
    vUp[0] = 0.0f;
    vUp[1] = 1.0f;
    vUp[2] = 0.0f;

    // Forward is -Z (default OpenGL)
    vForward[0] = 0.0f;
    vForward[1] = 0.0f;
    vForward[2] = -1.0f;
  }

  /////////////////////////////////////////////////////////////
  // Set Location
  inline void SetOrigin(const M3DVector3f vPoint) {
    m3dCopyVector3(vOrigin, vPoint);
  }

  inline void SetOrigin(float x, float y, float z) {
    vOrigin[0] = x;
    vOrigin[1] = y;
    vOrigin[2] = z;
  }
  void ApplyActorTransform(bool bRotationOnly = false) {
    M3DMatrix44f rotMat;

    GetMatrix(rotMat, bRotationOnly);

    // Apply rotation to the current matrix
    glMultMatrixf(rotMat);
  }
  inline void GetCameraOrientation(M3DMatrix44f m) {
    M3DVector3f x, z;

    // Make rotation matrix
    // Z vector is reversed
    z[0] = -vForward[0];
    z[1] = -vForward[1];
    z[2] = -vForward[2];

    // X vector = Y cross Z
    m3dCrossProduct(x, vUp, z);

// Matrix has no translation information and is
// transposed.... (rows instead of columns)
#define M(row, col) m[col * 4 + row]
    M(0, 0) = x[0];
    M(0, 1) = x[1];
    M(0, 2) = x[2];
    M(0, 3) = 0.0;
    M(1, 0) = vUp[0];
    M(1, 1) = vUp[1];
    M(1, 2) = vUp[2];
    M(1, 3) = 0.0;
    M(2, 0) = z[0];
    M(2, 1) = z[1];
    M(2, 2) = z[2];
    M(2, 3) = 0.0;
    M(3, 0) = 0.0;
    M(3, 1) = 0.0;
    M(3, 2) = 0.0;
    M(3, 3) = 1.0;
#undef M
  }
  inline void ApplyCameraTransform(bool bRotOnly = false) {
    M3DMatrix44f m;

    GetCameraOrientation(m);

    // Camera Transform
    glMultMatrixf(m);

    // If Rotation only, then do not do the translation
    if (!bRotOnly)
      glTranslatef(-vOrigin[0], -vOrigin[1], -vOrigin[2]);

    /*gluLookAt(vOrigin[0], vOrigin[1], vOrigin[2],
                            vOrigin[0] + vForward[0],
                            vOrigin[1] + vForward[1],
                            vOrigin[2] + vForward[2],
                            vUp[0], vUp[1], vUp[2]);
    */
  }
  inline void GetOrigin(M3DVector3f vPoint) { m3dCopyVector3(vPoint, vOrigin); }

  inline float GetOriginX(void) { return vOrigin[0]; }
  inline float GetOriginY(void) { return vOrigin[1]; }
  inline float GetOriginZ(void) { return vOrigin[2]; }

  /////////////////////////////////////////////////////////////
  // Set Forward Direction
  inline void SetForwardVector(const M3DVector3f vDirection) {
    m3dCopyVector3(vForward, vDirection);
  }

  inline void SetForwardVector(float x, float y, float z) {
    vForward[0] = x;
    vForward[1] = y;
    vForward[2] = z;
  }

  inline void GetForwardVector(M3DVector3f vVector) {
    m3dCopyVector3(vVector, vForward);
  }

  /////////////////////////////////////////////////////////////
  // Set Up Direction
  inline void SetUpVector(const M3DVector3f vDirection) {
    m3dCopyVector3(vUp, vDirection);
  }

  inline void SetUpVector(float x, float y, float z) {
    vUp[0] = x;
    vUp[1] = y;
    vUp[2] = z;
  }

  inline void GetUpVector(M3DVector3f vVector) { m3dCopyVector3(vVector, vUp); }

  /////////////////////////////////////////////////////////////
  // Get Axes
  inline void GetZAxis(M3DVector3f vVector) { GetForwardVector(vVector); }
  inline void GetYAxis(M3DVector3f vVector) { GetUpVector(vVector); }
  inline void GetXAxis(M3DVector3f vVector) {
    m3dCrossProduct3(vVector, vUp, vForward);
  }

  /////////////////////////////////////////////////////////////
  // Translate along orthonormal axis... world or local
  inline void TranslateWorld(float x, float y, float z) {
    vOrigin[0] += x;
    vOrigin[1] += y;
    vOrigin[2] += z;
  }

  inline void TranslateLocal(float x, float y, float z) {
    MoveForward(z);
    MoveUp(y);
    MoveRight(x);
  }

  /////////////////////////////////////////////////////////////
  // Move Forward (along Z axis)
  inline void MoveForward(float fDelta) {
    // Move along direction of front direction
    vOrigin[0] += vForward[0] * fDelta;
    vOrigin[1] += vForward[1] * fDelta;
    vOrigin[2] += vForward[2] * fDelta;
  }

  // Move along Y axis
  inline void MoveUp(float fDelta) {
    // Move along direction of up direction
    vOrigin[0] += vUp[0] * fDelta;
    vOrigin[1] += vUp[1] * fDelta;
    vOrigin[2] += vUp[2] * fDelta;
  }

  // Move along X axis
  inline void MoveRight(float fDelta) {
    // Move along direction of right vector
    M3DVector3f vCross;
    m3dCrossProduct3(vCross, vUp, vForward);

    vOrigin[0] += vCross[0] * fDelta;
    vOrigin[1] += vCross[1] * fDelta;
    vOrigin[2] += vCross[2] * fDelta;
  }

  ///////////////////////////////////////////////////////////////////////
  // Just assemble the matrix
  void GetMatrix(M3DMatrix44f matrix, bool bRotationOnly = false) {
    // Calculate the right side (x) vector, drop it right into the matrix
    M3DVector3f vXAxis;
    m3dCrossProduct3(vXAxis, vUp, vForward);

    // Set matrix column does not fill in the fourth value...
    m3dSetMatrixColumn44(matrix, vXAxis, 0);
    matrix[3] = 0.0f;

    // Y Column
    m3dSetMatrixColumn44(matrix, vUp, 1);
    matrix[7] = 0.0f;

    // Z Column
    m3dSetMatrixColumn44(matrix, vForward, 2);
    matrix[11] = 0.0f;

    // Translation (already done)
    if (bRotationOnly == true) {
      matrix[12] = 0.0f;
      matrix[13] = 0.0f;
      matrix[14] = 0.0f;
    } else
      m3dSetMatrixColumn44(matrix, vOrigin, 3);

    matrix[15] = 1.0f;
  }

  ////////////////////////////////////////////////////////////////////////
  // Assemble the camera matrix
  void GetCameraMatrix(M3DMatrix44f m, bool bRotationOnly = false) {
    M3DVector3f x, z;

    // Make rotation matrix
    // Z vector is reversed
    z[0] = -vForward[0];
    z[1] = -vForward[1];
    z[2] = -vForward[2];

    // X vector = Y cross Z
    m3dCrossProduct3(x, vUp, z);

// Matrix has no translation information and is
// transposed.... (rows instead of columns)
#define M(row, col) m[col * 4 + row]
    M(0, 0) = x[0];
    M(0, 1) = x[1];
    M(0, 2) = x[2];
    M(0, 3) = 0.0;
    M(1, 0) = vUp[0];
    M(1, 1) = vUp[1];
    M(1, 2) = vUp[2];
    M(1, 3) = 0.0;
    M(2, 0) = z[0];
    M(2, 1) = z[1];
    M(2, 2) = z[2];
    M(2, 3) = 0.0;
    M(3, 0) = 0.0;
    M(3, 1) = 0.0;
    M(3, 2) = 0.0;
    M(3, 3) = 1.0;
#undef M

    if (bRotationOnly)
      return;

    // Apply translation too
    M3DMatrix44f trans, M;
    m3dTranslationMatrix44(trans, -vOrigin[0], -vOrigin[1], -vOrigin[2]);

    m3dMatrixMultiply44(M, m, trans);

    // Copy result back into m
    memcpy(m, M, sizeof(float) * 16);
  }
  static float DetIJ(const M3DMatrix44f m, const int i, const int j) {
    int x, y, ii, jj;
    float ret, mat[3][3];

    x = 0;
    for (ii = 0; ii < 4; ii++) {
      if (ii == i)
        continue;
      y = 0;
      for (jj = 0; jj < 4; jj++) {
        if (jj == j)
          continue;
        mat[x][y] = m[(ii * 4) + jj];
        y++;
      }
      x++;
    }

    ret = mat[0][0] * (mat[1][1] * mat[2][2] - mat[2][1] * mat[1][2]);
    ret -= mat[0][1] * (mat[1][0] * mat[2][2] - mat[2][0] * mat[1][2]);
    ret += mat[0][2] * (mat[1][0] * mat[2][1] - mat[2][0] * mat[1][1]);

    return ret;
  }
  static double DetIJ(const M3DMatrix44d m, const int i, const int j) {
    int x, y, ii, jj;
    double ret, mat[3][3];

    x = 0;
    for (ii = 0; ii < 4; ii++) {
      if (ii == i)
        continue;
      y = 0;
      for (jj = 0; jj < 4; jj++) {
        if (jj == j)
          continue;
        mat[x][y] = m[(ii * 4) + jj];
        y++;
      }
      x++;
    }

    ret = mat[0][0] * (mat[1][1] * mat[2][2] - mat[2][1] * mat[1][2]);
    ret -= mat[0][1] * (mat[1][0] * mat[2][2] - mat[2][0] * mat[1][2]);
    ret += mat[0][2] * (mat[1][0] * mat[2][1] - mat[2][0] * mat[1][1]);

    return ret;
  }
  void m3dInvertMatrix44(M3DMatrix44f mInverse, const M3DMatrix44f m) {
    int i, j;
    float det, detij;

    // calculate 4x4 determinant
    det = 0.0f;
    for (i = 0; i < 4; i++) {
      det += (i & 0x1) ? (-m[i] * DetIJ(m, 0, i)) : (m[i] * DetIJ(m, 0, i));
    }
    det = 1.0f / det;

    // calculate inverse
    for (i = 0; i < 4; i++) {
      for (j = 0; j < 4; j++) {
        detij = DetIJ(m, j, i);
        mInverse[(i * 4) + j] =
            ((i + j) & 0x1) ? (-detij * det) : (detij * det);
      }
    }
  }

  ////////////////////////////////////////////////////////////////////////////
  ///
  // Invert matrix
  void m3dInvertMatrix44(M3DMatrix44d mInverse, const M3DMatrix44d m) {
    int i, j;
    double det, detij;

    // calculate 4x4 determinant
    det = 0.0;
    for (i = 0; i < 4; i++) {
      det += (i & 0x1) ? (-m[i] * DetIJ(m, 0, i)) : (m[i] * DetIJ(m, 0, i));
    }
    det = 1.0 / det;

    // calculate inverse
    for (i = 0; i < 4; i++) {
      for (j = 0; j < 4; j++) {
        detij = DetIJ(m, j, i);
        mInverse[(i * 4) + j] =
            ((i + j) & 0x1) ? (-detij * det) : (detij * det);
      }
    }
  }
  // Rotate around local Y
  void RotateLocalY(float fAngle) {
    M3DMatrix44f rotMat;

    // Just Rotate around the up vector
    // Create a rotation matrix around my Up (Y) vector
    m3dRotationMatrix44(rotMat, fAngle, vUp[0], vUp[1], vUp[2]);

    M3DVector3f newVect;

    // Rotate forward pointing vector (inlined 3x3 transform)
    newVect[0] = rotMat[0] * vForward[0] + rotMat[4] * vForward[1] +
                 rotMat[8] * vForward[2];
    newVect[1] = rotMat[1] * vForward[0] + rotMat[5] * vForward[1] +
                 rotMat[9] * vForward[2];
    newVect[2] = rotMat[2] * vForward[0] + rotMat[6] * vForward[1] +
                 rotMat[10] * vForward[2];
    m3dCopyVector3(vForward, newVect);
  }

  // Rotate around local Z
  void RotateLocalZ(float fAngle) {
    M3DMatrix44f rotMat;

    // Only the up vector needs to be rotated
    m3dRotationMatrix44(rotMat, fAngle, vForward[0], vForward[1], vForward[2]);

    M3DVector3f newVect;
    newVect[0] = rotMat[0] * vUp[0] + rotMat[4] * vUp[1] + rotMat[8] * vUp[2];
    newVect[1] = rotMat[1] * vUp[0] + rotMat[5] * vUp[1] + rotMat[9] * vUp[2];
    newVect[2] = rotMat[2] * vUp[0] + rotMat[6] * vUp[1] + rotMat[10] * vUp[2];
    m3dCopyVector3(vUp, newVect);
  }

  void RotateLocalX(float fAngle) {
    M3DMatrix33f rotMat;
    M3DVector3f localX;
    M3DVector3f rotVec;

    // Get the local X axis
    m3dCrossProduct3(localX, vUp, vForward);

    // Make a Rotation Matrix
    m3dRotationMatrix33(rotMat, fAngle, localX[0], localX[1], localX[2]);

    // Rotate Y, and Z
    m3dRotateVector(rotVec, vUp, rotMat);
    m3dCopyVector3(vUp, rotVec);

    m3dRotateVector(rotVec, vForward, rotMat);
    m3dCopyVector3(vForward, rotVec);
  }

  // Reset axes to make sure they are orthonormal. This should be called on
  // occasion if the matrix is long-lived and frequently transformed.
  void Normalize(void) {
    M3DVector3f vCross;

    // Calculate cross product of up and forward vectors
    m3dCrossProduct3(vCross, vUp, vForward);

    // Use result to recalculate forward vector
    m3dCrossProduct3(vForward, vCross, vUp);

    // Also check for unit length...
    m3dNormalizeVector3(vUp);
    m3dNormalizeVector3(vForward);
  }

  // Rotate in world coordinates...
  void RotateWorld(float fAngle, float x, float y, float z) {
    M3DMatrix44f rotMat;

    // Create the Rotation matrix
    m3dRotationMatrix44(rotMat, fAngle, x, y, z);

    M3DVector3f newVect;

    // Transform the up axis (inlined 3x3 rotation)
    newVect[0] = rotMat[0] * vUp[0] + rotMat[4] * vUp[1] + rotMat[8] * vUp[2];
    newVect[1] = rotMat[1] * vUp[0] + rotMat[5] * vUp[1] + rotMat[9] * vUp[2];
    newVect[2] = rotMat[2] * vUp[0] + rotMat[6] * vUp[1] + rotMat[10] * vUp[2];
    m3dCopyVector3(vUp, newVect);

    // Transform the forward axis
    newVect[0] = rotMat[0] * vForward[0] + rotMat[4] * vForward[1] +
                 rotMat[8] * vForward[2];
    newVect[1] = rotMat[1] * vForward[0] + rotMat[5] * vForward[1] +
                 rotMat[9] * vForward[2];
    newVect[2] = rotMat[2] * vForward[0] + rotMat[6] * vForward[1] +
                 rotMat[10] * vForward[2];
    m3dCopyVector3(vForward, newVect);
  }

  // Rotate around a local axis
  void RotateLocal(float fAngle, float x, float y, float z) {
    M3DVector3f vWorldVect;
    M3DVector3f vLocalVect;
    m3dLoadVector3(vLocalVect, x, y, z);

    LocalToWorld(vLocalVect, vWorldVect, true);
    RotateWorld(fAngle, vWorldVect[0], vWorldVect[1], vWorldVect[2]);
  }

  // Convert Coordinate Systems
  // This is pretty much, do the transformation represented by the rotation
  // and position on the point
  // Is it better to stick to the convention that the destination always comes
  // first, or use the conventions that "sounds" like the function...
  void LocalToWorld(const M3DVector3f vLocal, M3DVector3f vWorld,
                    bool bRotOnly = false) {
    // Create the rotation matrix based on the vectors
    M3DMatrix44f rotMat;

    GetMatrix(rotMat, true);

    // Do the rotation (inline it, and remove 4th column...)
    vWorld[0] =
        rotMat[0] * vLocal[0] + rotMat[4] * vLocal[1] + rotMat[8] * vLocal[2];
    vWorld[1] =
        rotMat[1] * vLocal[0] + rotMat[5] * vLocal[1] + rotMat[9] * vLocal[2];
    vWorld[2] =
        rotMat[2] * vLocal[0] + rotMat[6] * vLocal[1] + rotMat[10] * vLocal[2];

    // Translate the point
    if (!bRotOnly) {
      vWorld[0] += vOrigin[0];
      vWorld[1] += vOrigin[1];
      vWorld[2] += vOrigin[2];
    }
  }

  // Change world coordinates into "local" coordinates
  void WorldToLocal(const M3DVector3f vWorld, M3DVector3f vLocal) {
    ////////////////////////////////////////////////
    // Translate the origin
    M3DVector3f vNewWorld;
    vNewWorld[0] = vWorld[0] - vOrigin[0];
    vNewWorld[1] = vWorld[1] - vOrigin[1];
    vNewWorld[2] = vWorld[2] - vOrigin[2];

    // Create the rotation matrix based on the vectors
    M3DMatrix44f rotMat;
    M3DMatrix44f invMat;
    GetMatrix(rotMat, true);

    // Do the rotation based on inverted matrix
    m3dInvertMatrix44(invMat, rotMat);

    vLocal[0] = invMat[0] * vNewWorld[0] + invMat[4] * vNewWorld[1] +
                invMat[8] * vNewWorld[2];
    vLocal[1] = invMat[1] * vNewWorld[0] + invMat[5] * vNewWorld[1] +
                invMat[9] * vNewWorld[2];
    vLocal[2] = invMat[2] * vNewWorld[0] + invMat[6] * vNewWorld[1] +
                invMat[10] * vNewWorld[2];
  }

  /////////////////////////////////////////////////////////////////////////////
  // Transform a point by frame matrix
  void TransformPoint(M3DVector3f vPointSrc, M3DVector3f vPointDst) {
    M3DMatrix44f m;
    GetMatrix(m, false); // Rotate and translate
    vPointDst[0] = m[0] * vPointSrc[0] + m[4] * vPointSrc[1] +
                   m[8] * vPointSrc[2] + m[12]; // * v[3];
    vPointDst[1] = m[1] * vPointSrc[0] + m[5] * vPointSrc[1] +
                   m[9] * vPointSrc[2] + m[13]; // * v[3];
    vPointDst[2] = m[2] * vPointSrc[0] + m[6] * vPointSrc[1] +
                   m[10] * vPointSrc[2] + m[14]; // * v[3];
  }

  ////////////////////////////////////////////////////////////////////////////
  // Rotate a vector by frame matrix
  void RotateVector(M3DVector3f vVectorSrc, M3DVector3f vVectorDst) {
    M3DMatrix44f m;
    GetMatrix(m, true); // Rotate only

    vVectorDst[0] =
        m[0] * vVectorSrc[0] + m[4] * vVectorSrc[1] + m[8] * vVectorSrc[2];
    vVectorDst[1] =
        m[1] * vVectorSrc[0] + m[5] * vVectorSrc[1] + m[9] * vVectorSrc[2];
    vVectorDst[2] =
        m[2] * vVectorSrc[0] + m[6] * vVectorSrc[1] + m[10] * vVectorSrc[2];
  }
};

#define NUM_SPHERES 30
GLFrame spheres[NUM_SPHERES];
GLFrame frameCamera;

// Light and material Data
GLfloat fLightPos[4] = {-100.0f, 100.0f, 50.0f, 1.0f}; // Point source
GLfloat fNoLight[] = {0.0f, 0.0f, 0.0f, 0.0f};
GLfloat fLowLight[] = {0.25f, 0.25f, 0.25f, 1.0f};
GLfloat fBrightLight[] = {1.0f, 1.0f, 1.0f, 1.0f};

M3DMatrix44f mShadowMatrix;

#define GROUND_TEXTURE 0
#define TORUS_TEXTURE 1
#define SPHERE_TEXTURE 2
#define NUM_TEXTURES 3
GLuint textureObjects[NUM_TEXTURES];

const char *szTextureFiles[] = {RESOURCE_DIR "/grass.tga",
                                RESOURCE_DIR "/wood.tga",
                                RESOURCE_DIR "/orb.tga"};
bool timerActive = false;
bool isSpaceDown = false;
float headAngle = 0.0f;
float LHandAngle = 0.0f;
float RHandAngle = 0.0f;
float LHandTAngle = 0.0f;
float RHandTAngle = 0.0f;

float LLegAngle = 0.0f;
float RLegAngle = 0.0f;
float LLegTAngle = 0.0f;
float RLegTAngle = 0.0f;
//////////////////////////////////////////////////////////////////
// This function does any needed initialization on the rendering
// context.

#pragma pack(push, 1) // 開啟緊湊字節對齊
typedef struct {
  GLbyte idLength;
  GLbyte colorMapType;
  GLbyte imageType;
  GLshort colorMapOrigin;
  GLshort colorMapLength;
  GLbyte colorMapDepth;
  GLshort xOrigin;
  GLshort yOrigin;
  GLshort width;
  GLshort height;
  GLbyte bitsPerPixel;
  GLbyte imageDescriptor;
} TGAHEADER;
#pragma pack(pop) // 關閉緊湊字節對齊

GLbyte *gltLoadTGA(const char *filename, GLint *width, GLint *height,
                   GLint *components, GLenum *format) {
  FILE *file;
  TGAHEADER tgaHeader;
  unsigned long imageSize;
  short depth;
  GLbyte *pixels = NULL;

  // 打開文件
  file = fopen(filename, "rb");
  if (file == NULL)
    return NULL;

  // 讀取 TGA 文件頭
  fread(&tgaHeader, sizeof(TGAHEADER), 1, file);

  // 取得寬度、高度和色深
  *width = tgaHeader.width;
  *height = tgaHeader.height;
  depth = tgaHeader.bitsPerPixel / 8; // Color depth in bytes per pixel

  // 確定組件數和格式
  switch (depth) {
  case 3: // 24-bit
    *components = GL_RGB;
    *format = GL_BGR;
    break;
  case 4: // 32-bit
    *components = GL_RGBA;
    *format = GL_BGRA;
    break;
  default:
    fclose(file);
    return NULL;
  }

  // 計算圖像大小（字節數）
  imageSize = *width * *height * depth;

  // 分配內存以存儲像素數據
  pixels = (GLbyte *)malloc(imageSize);

  // 讀取像素數據
  if (fread(pixels, imageSize, 1, file) != 1) {
    free(pixels);
    fclose(file);
    return NULL;
  }

  // 關閉文件
  fclose(file);

  return pixels;
}
void SetupRC() {
  M3DVector3f vPoints[3] = {
      {0.0f, -0.4f, 0.0f}, {10.0f, -0.4f, 0.0f}, {5.0f, -0.4f, -5.0f}};
  int iSphere;
  int i;

  // Grayish background
  glClearColor(fLowLight[0], fLowLight[1], fLowLight[2], fLowLight[3]);

  // Clear stencil buffer with zero, increment by one whenever anybody
  // draws into it. When stencil function is enabled, only write where
  // stencil value is zero. This prevents the transparent shadow from drawing
  // over itself
  glStencilOp(GL_INCR, GL_INCR, GL_INCR);
  glClearStencil(0);
  glStencilFunc(GL_EQUAL, 0x0, 0x01);

  // Cull backs of polygons
  glCullFace(GL_BACK);
  glFrontFace(GL_CCW);
  glEnable(GL_CULL_FACE);
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_MULTISAMPLE_ARB);

  // Setup light parameters
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, fNoLight);
  glLightfv(GL_LIGHT0, GL_AMBIENT, fLowLight);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, fBrightLight);
  glLightfv(GL_LIGHT0, GL_SPECULAR, fBrightLight);
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);

  // Calculate shadow matrix
  M3DVector4f pPlane;
  m3dGetPlaneEquation(pPlane, vPoints[0], vPoints[1], vPoints[2]);
  m3dMakePlanarShadowMatrix(mShadowMatrix, pPlane, fLightPos);

  // Mostly use material tracking
  glEnable(GL_COLOR_MATERIAL);
  glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
  glMaterialfv(GL_FRONT, GL_SPECULAR, fBrightLight);
  glMateriali(GL_FRONT, GL_SHININESS, 128);

  // Randomly place the sphere inhabitants
  for (iSphere = 0; iSphere < NUM_SPHERES; iSphere++) {
    // Pick a random location between -20 and 20 at .1 increments
    spheres[iSphere].SetOrigin(((float)((rand() % 400) - 200) * 0.1f), 0.0,
                               (float)((rand() % 400) - 200) * 0.1f);
  }

  // Set up texture maps
  glEnable(GL_TEXTURE_2D);
  glGenTextures(NUM_TEXTURES, textureObjects);
  glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

  for (i = 0; i < NUM_TEXTURES; i++) {
    GLbyte *pBytes;
    GLint iWidth, iHeight, iComponents;
    GLenum eFormat;

    glBindTexture(GL_TEXTURE_2D, textureObjects[i]);

    // Load this texture map
    pBytes = gltLoadTGA(szTextureFiles[i], &iWidth, &iHeight, &iComponents,
                        &eFormat);
    gluBuild2DMipmaps(GL_TEXTURE_2D, iComponents, iWidth, iHeight, eFormat,
                      GL_UNSIGNED_BYTE, pBytes);
    free(pBytes);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
                    GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  }
}

////////////////////////////////////////////////////////////////////////
// Do shutdown for the rendering context
void ShutdownRC(void) {
  // Delete the textures
  glDeleteTextures(NUM_TEXTURES, textureObjects);
}

///////////////////////////////////////////////////////////
// Draw the ground as a series of triangle strips
void DrawGround(void) {
  GLfloat fExtent = 20.0f;
  GLfloat fStep = 1.0f;
  GLfloat y = -0.4f;
  GLfloat iStrip, iRun;
  GLfloat s = 0.0f;
  GLfloat t = 0.0f;
  GLfloat texStep = 1.0f / (fExtent * .075f);

  glBindTexture(GL_TEXTURE_2D, textureObjects[GROUND_TEXTURE]);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

  for (iStrip = -fExtent; iStrip <= fExtent; iStrip += fStep) {
    t = 0.0f;
    glBegin(GL_TRIANGLE_STRIP);

    for (iRun = fExtent; iRun >= -fExtent; iRun -= fStep) {
      glTexCoord2f(s, t);
      glNormal3f(0.0f, 1.0f, 0.0f); // All Point up
      glVertex3f(iStrip, y, iRun);

      glTexCoord2f(s + texStep, t);
      glNormal3f(0.0f, 1.0f, 0.0f); // All Point up
      glVertex3f(iStrip + fStep, y, iRun);

      t += texStep;
    }
    glEnd();
    s += texStep;
  }
}
void DrawBox(float x, float y, float z, float w, float h, float d) {
  glPushMatrix();
  float ScaleDelta = 0.05f;
  glScalef(w * 0.5, h * 0.5, d * 0.5);
  float TranslateDleta = 0.2f;
  glTranslatef(x, y, z);
  glBegin(GL_QUADS);
  // Front Face
  glTexCoord2f(0.0f, 0.0f);
  glVertex3f(-1.0f, -1.0f, 1.0f);
  glTexCoord2f(1.0f, 0.0f);
  glVertex3f(1.0f, -1.0f, 1.0f);
  glTexCoord2f(1.0f, 1.0f);
  glVertex3f(1.0f, 1.0f, 1.0f);
  glTexCoord2f(0.0f, 1.0f);
  glVertex3f(-1.0f, 1.0f, 1.0f);

  // Back Face
  glTexCoord2f(1.0f, 0.0f);
  glVertex3f(-1.0f, -1.0f, -1.0f);
  glTexCoord2f(1.0f, 1.0f);
  glVertex3f(-1.0f, 1.0f, -1.0f);
  glTexCoord2f(0.0f, 1.0f);
  glVertex3f(1.0f, 1.0f, -1.0f);
  glTexCoord2f(0.0f, 0.0f);
  glVertex3f(1.0f, -1.0f, -1.0f);

  // Top Face
  glTexCoord2f(0.0f, 1.0f);
  glVertex3f(-1.0f, 1.0f, -1.0f);
  glTexCoord2f(0.0f, 0.0f);
  glVertex3f(-1.0f, 1.0f, 1.0f);
  glTexCoord2f(1.0f, 0.0f);
  glVertex3f(1.0f, 1.0f, 1.0f);
  glTexCoord2f(1.0f, 1.0f);
  glVertex3f(1.0f, 1.0f, -1.0f);

  // Bottom Face
  glTexCoord2f(1.0f, 1.0f);
  glVertex3f(-1.0f, -1.0f, -1.0f);
  glTexCoord2f(0.0f, 1.0f);
  glVertex3f(1.0f, -1.0f, -1.0f);
  glTexCoord2f(0.0f, 0.0f);
  glVertex3f(1.0f, -1.0f, 1.0f);
  glTexCoord2f(1.0f, 0.0f);
  glVertex3f(-1.0f, -1.0f, 1.0f);

  // Right face
  glTexCoord2f(1.0f, 0.0f);
  glVertex3f(1.0f, -1.0f, -1.0f);
  glTexCoord2f(1.0f, 1.0f);
  glVertex3f(1.0f, 1.0f, -1.0f);
  glTexCoord2f(0.0f, 1.0f);
  glVertex3f(1.0f, 1.0f, 1.0f);
  glTexCoord2f(0.0f, 0.0f);
  glVertex3f(1.0f, -1.0f, 1.0f);

  // Left Face
  glTexCoord2f(0.0f, 0.0f);
  glVertex3f(-1.0f, -1.0f, -1.0f);
  glTexCoord2f(1.0f, 0.0f);
  glVertex3f(-1.0f, -1.0f, 1.0f);
  glTexCoord2f(1.0f, 1.0f);
  glVertex3f(-1.0f, 1.0f, 1.0f);
  glTexCoord2f(0.0f, 1.0f);
  glVertex3f(-1.0f, 1.0f, -1.0f);

  glEnd();
  glPopMatrix();
}
void gltDrawSphere(GLfloat fRadius, GLint iSlices, GLint iStacks) {
  GLfloat drho = (GLfloat)(3.141592653589) / (GLfloat)iStacks;
  GLfloat dtheta = 2.0f * (GLfloat)(3.141592653589) / (GLfloat)iSlices;
  GLfloat ds = 1.0f / (GLfloat)iSlices;
  GLfloat dt = 1.0f / (GLfloat)iStacks;
  GLfloat t = 1.0f;
  GLfloat s = 0.0f;
  GLint i, j; // Looping variables

  for (i = 0; i < iStacks; i++) {
    GLfloat rho = (GLfloat)i * drho;
    GLfloat srho = (GLfloat)(sin(rho));
    GLfloat crho = (GLfloat)(cos(rho));
    GLfloat srhodrho = (GLfloat)(sin(rho + drho));
    GLfloat crhodrho = (GLfloat)(cos(rho + drho));

    // Many sources of OpenGL sphere drawing code uses a triangle fan
    // for the caps of the sphere. This however introduces texturing
    // artifacts at the poles on some OpenGL implementations
    glBegin(GL_TRIANGLE_STRIP);
    s = 0.0f;
    for (j = 0; j <= iSlices; j++) {
      GLfloat theta = (j == iSlices) ? 0.0f : j * dtheta;
      GLfloat stheta = (GLfloat)(-sin(theta));
      GLfloat ctheta = (GLfloat)(cos(theta));

      GLfloat x = stheta * srho;
      GLfloat y = ctheta * srho;
      GLfloat z = crho;

      glTexCoord2f(s, t);
      glNormal3f(x, y, z);
      glVertex3f(x * fRadius, y * fRadius, z * fRadius);

      x = stheta * srhodrho;
      y = ctheta * srhodrho;
      z = crhodrho;
      glTexCoord2f(s, t - dt);
      s += ds;
      glNormal3f(x, y, z);
      glVertex3f(x * fRadius, y * fRadius, z * fRadius);
    }
    glEnd();

    t -= dt;
  }
}

void setupTextureMatrix(float scaleX, float scaleY) {
  glMatrixMode(GL_TEXTURE);
  glLoadIdentity();
  glScalef(scaleX, scaleY, 1.0f); // 設置紋理的縮放因子
  glMatrixMode(GL_MODELVIEW);     // 恢復到模型視圖矩陣
}

///////////////////////////////////////////////////////////////////////
// Draw random inhabitants and the rotating torus/sphere duo
void DrawInhabitants(GLint nShadow) {
  static GLfloat yRot = 0.0f; // Rotation angle for animation
  GLint i;

  if (nShadow == 0) {
    yRot += 0.5f;
    glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
  } else
    glColor4f(0.00f, 0.00f, 0.00f, .6f); // Shadow color

  // Draw the randomly located spheres
  for (i = 0; i < NUM_SPHERES; i++) {
    glPushMatrix();
    spheres[i].ApplyActorTransform();
    gltDrawSphere(0.3, 21, 11);
    glPopMatrix();
  }

  glPushMatrix();
  glTranslatef(0.0f, 0.1f, -2.5f);
  glPushMatrix();
  glRotatef(-yRot * 2.0f, 0.0f, 1.0f, 0.0f);
  glTranslatef(1.0f, 0.0f, 0.0f);
  gltDrawSphere(0.3, 21, 11);

  glPopMatrix();

  if (nShadow == 0) {
    // Torus alone will be specular
    glMaterialfv(GL_FRONT, GL_SPECULAR, fBrightLight);
  }

  glRotatef(yRot, 0.0f, 1.0f, 0.0f);
  glBindTexture(GL_TEXTURE_2D, textureObjects[TORUS_TEXTURE]);
  //   --------------------------

  glPushMatrix();
  glScalef(0.005, 0.005, 0.005);
  glPushMatrix();
  glTranslatef(0, 30, 0);
  // glColor3ub(255, 255, 0);

  DrawBox(0, 0, 0, 35, 50, 10); // body

  glPushMatrix();
  glRotatef(headAngle, 0, 0, 1); // Head rotation
  glTranslatef(0, 40, 0);
  DrawBox(0, 0, 0, 20, 20, 10); // Head
  glPopMatrix();

  glPushMatrix();
  glRotatef(LHandAngle, 0, 0, 1); // LHand rotation
  glTranslatef(-40, 10, 0);
  DrawBox(0, 0, 0, 40, 10, 10); // LHand

  glPushMatrix();
  glRotatef(LHandTAngle, 0, 0, 1); // LHand_T rotation
  glTranslatef(-50, 0, 0);
  DrawBox(0, 0, 0, 40, 10, 10); // LHand_T
  glPopMatrix();

  glPopMatrix();

  glPushMatrix();
  glRotatef(RHandAngle, 0, 0, 1); // RHand rotation
  glTranslatef(40, 10, 0);
  DrawBox(0, 0, 0, 40, 10, 10); // RHand

  glPushMatrix();
  glRotatef(RHandTAngle, 0, 0, 1); // RHand_T rotation
  glTranslatef(50, 0, 0);
  DrawBox(0, 0, 0, 40, 10, 10); // RHand_T
  glPopMatrix();

  glPopMatrix();

  glPushMatrix();
  glRotatef(LLegAngle, 0, 0, 1); // LLeg rotation
  glTranslatef(-10, -60, 0);
  DrawBox(0, 0, 0, 10, 40, 10); // LLeg

  glPushMatrix();
  glRotatef(LLegTAngle, 0, 0, 1); // LLeg_T rotation
  glTranslatef(0, -50, 0);
  DrawBox(0, 0, 0, 10, 40, 10); // LLeg_T
  glPopMatrix();

  glPopMatrix();

  glPushMatrix();
  glRotatef(RLegAngle, 0, 0, 1); // RLeg_T rotation
  glTranslatef(10, -60, 0);
  DrawBox(0, 0, 0, 10, 40, 10); // RLeg

  glPushMatrix();
  glRotatef(RLegTAngle, 0, 0, 1); // RLeg_T rotation
  glTranslatef(0, -50, 0);
  DrawBox(0, 0, 0, 10, 40, 10); // RLeg_T
  glPopMatrix();

  glPopMatrix();
  glPopMatrix();
  glPopMatrix();

  //------------------------------------------
  //   glMaterialfv(GL_FRONT, GL_SPECULAR, fNoLight);
  glPopMatrix();
}

// Called to draw scene
void RenderScene(void) {
  // Clear the window with current clearing color
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

  glPushMatrix();
  frameCamera.ApplyCameraTransform();

  // Position light before any other transformations
  glLightfv(GL_LIGHT0, GL_POSITION, fLightPos);

  // Draw the ground
  glColor3f(1.0f, 1.0f, 1.0f);
  DrawGround();

  // Draw shadows first
  glDisable(GL_DEPTH_TEST);
  glDisable(GL_LIGHTING);
  glDisable(GL_TEXTURE_2D);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glEnable(GL_STENCIL_TEST);
  glPushMatrix();
  glMultMatrixf(mShadowMatrix);
  DrawInhabitants(1);
  glPopMatrix();
  glDisable(GL_STENCIL_TEST);
  glDisable(GL_BLEND);
  glEnable(GL_LIGHTING);
  glEnable(GL_TEXTURE_2D);
  glEnable(GL_DEPTH_TEST);

  // Draw inhabitants normally
  DrawInhabitants(0);

  glPopMatrix();

  // Do the buffer Swap
  glutSwapBuffers();
}

// Respond to arrow keys by moving the camera frame of reference
void SpecialKeys(int key, int x, int y) {
  if (key == GLUT_KEY_UP)
    frameCamera.MoveForward(0.1f);

  if (key == GLUT_KEY_DOWN)
    frameCamera.MoveForward(-0.1f);

  if (key == GLUT_KEY_LEFT)
    frameCamera.RotateLocalY(0.1f);

  if (key == GLUT_KEY_RIGHT)
    frameCamera.RotateLocalY(-0.1f);

  // Refresh the Window
  glutPostRedisplay();
}

///////////////////////////////////////////////////////////
// Called by GLUT library when idle (window not being
// resized or moved)
void TimerFunction(int value) {
  // Redraw the scene with new coordinates
  if (timerActive) {

    glutPostRedisplay();
  }
  glutTimerFunc(3, TimerFunction, 1);
}

void ChangeSize(int w, int h) {
  GLfloat fAspect;

  // Prevent a divide by zero, when window is too short
  // (you cant make a window of zero width).
  if (h == 0)
    h = 1;

  glViewport(0, 0, w, h);

  fAspect = (GLfloat)w / (GLfloat)h;

  // Reset the coordinate system before modifying
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  // Set the clipping volume
  gluPerspective(35.0f, fAspect, 1.0f, 50.0f);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
}
void KeyboardFunc(unsigned char key, int x, int y) {
  switch (key) {
  case 'z':
    LLegAngle += 2;
    break;
  case 'x':
    LLegTAngle += 2;
    break;
  case 'c':
    RLegAngle += 2;
    break;
  case 'v':
    RLegTAngle += 2;
    break;
  case 'a':
    LHandAngle += 2;
    break;
  case 's':
    LHandTAngle += 2;
    break;
  case 'd':
    RHandAngle += 2;
    break;
  case 'f':
    RHandTAngle += 2;
    break;
  case 'w':
    headAngle += 2;
    break;

  case 'b':
    LLegAngle -= 2;
    break;
  case 'n':
    LLegTAngle -= 2;
    break;
  case 'm':
    RLegAngle -= 2;
    break;
  case ',':
    RLegTAngle -= 2;
    break;
  case 'h':
    LHandAngle -= 2;
    break;
  case 'j':
    LHandTAngle -= 2;
    break;
  case 'k':
    RHandAngle -= 2;
    break;
  case 'l':
    RHandTAngle -= 2;
    break;
  case 'i':
    headAngle -= 2;
    break;
  }
  switch (key) {

  case ' ':
    bool prev = isSpaceDown;
    isSpaceDown = !isSpaceDown;
    if ((isSpaceDown == true && prev == false) ||
        (isSpaceDown == false && prev == true)) {
      std::cout << "change" << std::endl;
      timerActive = !timerActive;
    }
    break;
  }
}
int main(int argc, char *argv[]) {
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | GLUT_STENCIL);
  glutInitWindowSize(800, 600);
  glutCreateWindow("OpenGL SphereWorld Demo + Texture Maps");
  glutReshapeFunc(ChangeSize);
  glutDisplayFunc(RenderScene);
  glutSpecialFunc(SpecialKeys);
  glutKeyboardFunc(KeyboardFunc);
  SetupRC();
  glutTimerFunc(33, TimerFunction, 1);

  glutMainLoop();

  ShutdownRC();

  return 0;
}
