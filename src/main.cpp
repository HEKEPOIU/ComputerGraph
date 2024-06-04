#include "freeglut_std.h"
#include <array>
#include <math.h>

// Rotation amounts
static GLfloat xRot = 0.0f;
static GLfloat yRot = 0.0f;

// These values need to be available globally
// Light values and coordinates
GLfloat ambientLight[] = {0.3f, 0.3f, 0.3f, 1.0f};
GLfloat diffuseLight[] = {0.7f, 0.7f, 0.7f, 1.0f};
GLfloat specular[] = {1.0f, 1.0f, 1.0f, 1.0f};
std::array<std::array<float, 4>, 4> lightPos = {
    {{{-75.0f, 150.0f, -50.0f, 0.0f}},
     {{-250.0f, 150.0f, -50.0f, 0.0f}},
     {{-150.0f, 150.0f, -50.0f, 0.0f}},
     {{305.0f, 250.0f, -50.0f, 0.0f}}}};
GLfloat specref[] = {1.0f, 1.0f, 1.0f, 1.0f};
std::array<bool, 4> bLight = {true, false, false, false};

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
M3DMatrix44f shadowMat;
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

void SetupRC();

void KeyBoradFunc(unsigned char key, int x, int y) {
  switch (key) {
  case ' ':
    xRot = 0.0f;
    yRot = 0.0f;
    break;
  case '1':
    for (int i = 0; i < bLight.size(); i++) {
      if (i == 0) {
        bLight[0] = true;
      } else {
        bLight[i] = false;
      }
    }
    SetupRC();
    break;
  case '2':

    for (int i = 0; i < bLight.size(); i++) {
      if (i == 1) {
        bLight[i] = true;
      } else {
        bLight[i] = false;
      }
    }
    SetupRC();
    break;
  case '3':

    for (int i = 0; i < bLight.size(); i++) {
      if (i == 2) {
        bLight[i] = true;
      } else {
        bLight[i] = false;
      }
    }
    SetupRC();
    break;
  case '4':

    for (int i = 0; i < bLight.size(); i++) {
      if (i == 3) {
        bLight[i] = true;
      } else {
        bLight[i] = false;
      }
    }
    SetupRC();
    break;
  default:
    break;
  }
  glutPostRedisplay();
}

////////////////////////////////////////////////
// This function just specifically draws the jet
void DrawJet(int nShadow) {
  M3DVector3f vNormal; // Storeage for calculated surface normal

  // Nose Cone /////////////////////////////
  // Set material color, note we only have to set to black
  // for the shadow once
  if (nShadow == 0)
    glColor3ub(255, 0, 0);
  else
    glColor3ub(50, 0, 0);

  // Nose Cone - Points straight down
  // Set material color
  glBegin(GL_TRIANGLES);
  glNormal3f(0.0f, -1.0f, 0.0f);
  glNormal3f(0.0f, -1.0f, 0.0f);
  glVertex3f(0.0f, 0.0f, 60.0f);
  glVertex3f(-15.0f, 0.0f, 30.0f);
  glVertex3f(15.0f, 0.0f, 30.0f);

  // Verticies for this panel
  {
    M3DVector3f vPoints[3] = {
        {15.0f, 0.0f, 30.0f}, {0.0f, 15.0f, 30.0f}, {0.0f, 0.0f, 60.0f}};

    // Calculate the normal for the plane
    m3dFindNormal(vNormal, vPoints[0], vPoints[1], vPoints[2]);
    glNormal3fv(vNormal);
    glVertex3fv(vPoints[0]);
    glVertex3fv(vPoints[1]);
    glVertex3fv(vPoints[2]);
  }

  {
    M3DVector3f vPoints[3] = {
        {0.0f, 0.0f, 60.0f}, {0.0f, 15.0f, 30.0f}, {-15.0f, 0.0f, 30.0f}};

    m3dFindNormal(vNormal, vPoints[0], vPoints[1], vPoints[2]);
    glNormal3fv(vNormal);
    glVertex3fv(vPoints[0]);
    glVertex3fv(vPoints[1]);
    glVertex3fv(vPoints[2]);
  }

  // Body of the Plane ////////////////////////
  {
    M3DVector3f vPoints[3] = {
        {-15.0f, 0.0f, 30.0f}, {0.0f, 15.0f, 30.0f}, {0.0f, 0.0f, -56.0f}};

    m3dFindNormal(vNormal, vPoints[0], vPoints[1], vPoints[2]);
    glNormal3fv(vNormal);
    glVertex3fv(vPoints[0]);
    glVertex3fv(vPoints[1]);
    glVertex3fv(vPoints[2]);
  }

  {
    M3DVector3f vPoints[3] = {
        {0.0f, 0.0f, -56.0f}, {0.0f, 15.0f, 30.0f}, {15.0f, 0.0f, 30.0f}};

    m3dFindNormal(vNormal, vPoints[0], vPoints[1], vPoints[2]);
    glNormal3fv(vNormal);
    glVertex3fv(vPoints[0]);
    glVertex3fv(vPoints[1]);
    glVertex3fv(vPoints[2]);
  }

  glNormal3f(0.0f, -1.0f, 0.0f);
  glVertex3f(15.0f, 0.0f, 30.0f);
  glVertex3f(-15.0f, 0.0f, 30.0f);
  glVertex3f(0.0f, 0.0f, -56.0f);

  ///////////////////////////////////////////////
  // Left wing
  // Large triangle for bottom of wing
  {
    M3DVector3f vPoints[3] = {
        {0.0f, 2.0f, 27.0f}, {-60.0f, 2.0f, -8.0f}, {60.0f, 2.0f, -8.0f}};

    m3dFindNormal(vNormal, vPoints[0], vPoints[1], vPoints[2]);
    glNormal3fv(vNormal);
    glVertex3fv(vPoints[0]);
    glVertex3fv(vPoints[1]);
    glVertex3fv(vPoints[2]);
  }

  {
    M3DVector3f vPoints[3] = {
        {60.0f, 2.0f, -8.0f}, {0.0f, 7.0f, -8.0f}, {0.0f, 2.0f, 27.0f}};

    m3dFindNormal(vNormal, vPoints[0], vPoints[1], vPoints[2]);
    glNormal3fv(vNormal);
    glVertex3fv(vPoints[0]);
    glVertex3fv(vPoints[1]);
    glVertex3fv(vPoints[2]);
  }

  {
    M3DVector3f vPoints[3] = {
        {60.0f, 2.0f, -8.0f}, {-60.0f, 2.0f, -8.0f}, {0.0f, 7.0f, -8.0f}};

    m3dFindNormal(vNormal, vPoints[0], vPoints[1], vPoints[2]);
    glNormal3fv(vNormal);
    glVertex3fv(vPoints[0]);
    glVertex3fv(vPoints[1]);
    glVertex3fv(vPoints[2]);
  }

  {
    M3DVector3f vPoints[3] = {
        {0.0f, 2.0f, 27.0f}, {0.0f, 7.0f, -8.0f}, {-60.0f, 2.0f, -8.0f}};

    m3dFindNormal(vNormal, vPoints[0], vPoints[1], vPoints[2]);
    glNormal3fv(vNormal);
    glVertex3fv(vPoints[0]);
    glVertex3fv(vPoints[1]);
    glVertex3fv(vPoints[2]);
  }

  // Tail section///////////////////////////////
  // Bottom of back fin
  glNormal3f(0.0f, -1.0f, 0.0f);
  glVertex3f(-30.0f, -0.50f, -57.0f);
  glVertex3f(30.0f, -0.50f, -57.0f);
  glVertex3f(0.0f, -0.50f, -40.0f);

  {
    M3DVector3f vPoints[3] = {
        {0.0f, -0.5f, -40.0f}, {30.0f, -0.5f, -57.0f}, {0.0f, 4.0f, -57.0f}};

    m3dFindNormal(vNormal, vPoints[0], vPoints[1], vPoints[2]);
    glNormal3fv(vNormal);
    glVertex3fv(vPoints[0]);
    glVertex3fv(vPoints[1]);
    glVertex3fv(vPoints[2]);
  }

  {
    M3DVector3f vPoints[3] = {
        {0.0f, 4.0f, -57.0f}, {-30.0f, -0.5f, -57.0f}, {0.0f, -0.5f, -40.0f}};

    m3dFindNormal(vNormal, vPoints[0], vPoints[1], vPoints[2]);
    glNormal3fv(vNormal);
    glVertex3fv(vPoints[0]);
    glVertex3fv(vPoints[1]);
    glVertex3fv(vPoints[2]);
  }

  {
    M3DVector3f vPoints[3] = {
        {30.0f, -0.5f, -57.0f}, {-30.0f, -0.5f, -57.0f}, {0.0f, 4.0f, -57.0f}};

    m3dFindNormal(vNormal, vPoints[0], vPoints[1], vPoints[2]);
    glNormal3fv(vNormal);
    glVertex3fv(vPoints[0]);
    glVertex3fv(vPoints[1]);
    glVertex3fv(vPoints[2]);
  }

  {
    M3DVector3f vPoints[3] = {
        {0.0f, 0.5f, -40.0f}, {3.0f, 0.5f, -57.0f}, {0.0f, 25.0f, -65.0f}};

    m3dFindNormal(vNormal, vPoints[0], vPoints[1], vPoints[2]);
    glNormal3fv(vNormal);
    glVertex3fv(vPoints[0]);
    glVertex3fv(vPoints[1]);
    glVertex3fv(vPoints[2]);
  }

  {
    M3DVector3f vPoints[3] = {
        {0.0f, 25.0f, -65.0f}, {-3.0f, 0.5f, -57.0f}, {0.0f, 0.5f, -40.0f}};

    m3dFindNormal(vNormal, vPoints[0], vPoints[1], vPoints[2]);
    glNormal3fv(vNormal);
    glVertex3fv(vPoints[0]);
    glVertex3fv(vPoints[1]);
    glVertex3fv(vPoints[2]);
  }

  {
    M3DVector3f vPoints[3] = {
        {3.0f, 0.5f, -57.0f}, {-3.0f, 0.5f, -57.0f}, {0.0f, 25.0f, -65.0f}};

    m3dFindNormal(vNormal, vPoints[0], vPoints[1], vPoints[2]);
    glNormal3fv(vNormal);
    glVertex3fv(vPoints[0]);
    glVertex3fv(vPoints[1]);
    glVertex3fv(vPoints[2]);
  }

  glEnd();
}

// Called to draw scene
void RenderScene(void) {
  // Clear the window with current clearing color
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // Draw the ground, we do manual shading to a darker green
  // in the background to give the illusion of depth
  glBegin(GL_QUADS);
  glColor3ub(0, 32, 0); // light green ground
  glVertex3f(400.0f, -150.0f, -200.0f);
  glVertex3f(-400.0f, -150.0f, -200.0f);
  glColor3ub(0, 255, 0); // make it in green gradient
  glVertex3f(-400.0f, -150.0f, 200.0f);
  glVertex3f(400.0f, -150.0f, 200.0f);
  glEnd();

  // Save the matrix state and do the rotations
  glPushMatrix();

  // Draw jet at new orientation, put light in correct position
  // before rotating the jet
  glEnable(GL_LIGHTING);
  for (int i = 0; i < bLight.size(); i++) {
    if (bLight[i]) {
      glEnable(GL_LIGHT0 + i);
      glLightfv(GL_LIGHT0 + i, GL_POSITION, lightPos[i].data());
      glRotatef(xRot, 1.0f, 0.0f, 0.0f);
      glRotatef(yRot, 0.0f, 1.0f, 0.0f);
    } else {
      glDisable(GL_LIGHT0 + i);
    }
  }

  DrawJet(0);

  // Restore original matrix state
  glPopMatrix();

  // Get ready to draw the shadow and the ground
  // First disable lighting and save the projection state
  glDisable(GL_DEPTH_TEST);
  glDisable(GL_LIGHTING);
  glPushMatrix();

  // Multiply by shadow projection matrix
  glMultMatrixf((GLfloat *)shadowMat);

  // Now rotate the jet around in the new flattend space
  glRotatef(xRot, 1.0f, 0.0f, 0.0f);
  glRotatef(yRot, 0.0f, 1.0f, 0.0f);

  // Pass true to indicate drawing shadow
  DrawJet(1);

  // Restore the projection to normal
  glPopMatrix();

  // Draw the light source
  for (int i = 0; i < bLight.size(); i++) {
    if (bLight[i]) {
      glPushMatrix();
      glTranslatef(lightPos[i][0], lightPos[i][1], lightPos[i][2]);
      glColor3ub(255, 255, 0);
      glutSolidSphere(5.0f, 10, 10);
      glPopMatrix();
    }
  }

  // Restore lighting state variables
  glEnable(GL_DEPTH_TEST);

  // Display the results
  glutSwapBuffers();
}

// This function does any needed initialization on the rendering
// context.
void SetupRC() {
  // Any three points on the ground (counter clockwise order)
  M3DVector3f points[3] = {{-30.0f, -149.0f, -20.0f},
                           {-30.0f, -149.0f, 20.0f},
                           {40.0f, -149.0f, 20.0f}};

  glEnable(GL_DEPTH_TEST); // Hidden surface removal
  glFrontFace(GL_CCW);     // Counter clock-wise polygons face out
  glEnable(GL_CULL_FACE);  // Do not calculate inside of jet

  // Setup and enable light 0
  glLightfv(GL_LIGHT0, GL_AMBIENT, ambientLight);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseLight);
  glLightfv(GL_LIGHT0, GL_SPECULAR, specular);
  glLightfv(GL_LIGHT0, GL_POSITION, lightPos[0].data());
  glEnable(GL_LIGHT0);

  glLightfv(GL_LIGHT1, GL_AMBIENT, ambientLight);
  glLightfv(GL_LIGHT1, GL_DIFFUSE, diffuseLight);
  glLightfv(GL_LIGHT1, GL_SPECULAR, specular);
  glLightfv(GL_LIGHT1, GL_POSITION, lightPos[1].data());
  // glEnable(GL_LIGHT1);

  glLightfv(GL_LIGHT2, GL_AMBIENT, ambientLight);
  glLightfv(GL_LIGHT2, GL_DIFFUSE, diffuseLight);
  glLightfv(GL_LIGHT2, GL_SPECULAR, specular);
  glLightfv(GL_LIGHT2, GL_POSITION, lightPos[2].data());
  // glEnable(GL_LIGHT2);

  glLightfv(GL_LIGHT3, GL_AMBIENT, ambientLight);
  glLightfv(GL_LIGHT3, GL_DIFFUSE, diffuseLight);
  glLightfv(GL_LIGHT3, GL_SPECULAR, specular);
  glLightfv(GL_LIGHT3, GL_POSITION, lightPos[3].data());
  // glEnable(GL_LIGHT3);

  // Enable color tracking
  glEnable(GL_COLOR_MATERIAL);

  // Set Material properties to follow glColor values
  glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);

  // All materials hereafter have full specular reflectivity
  // with a high shine
  glMaterialfv(GL_FRONT, GL_SPECULAR, specref);
  glMateriali(GL_FRONT, GL_SHININESS, 128);

  // Light blue background
  glClearColor(0.0f, 0.0f, 1.0f, 1.0f);

  // Get the plane equation from three points on the ground
  M3DVector4f vPlaneEquation;
  m3dGetPlaneEquation(vPlaneEquation, points[0], points[1], points[2]);

  // Calculate projection matrix to draw shadow on the ground

  for (int i = 0; i < bLight.size(); i++) {
    if (bLight[i]) {
      m3dMakePlanarShadowMatrix(shadowMat, vPlaneEquation, lightPos[i].data());

    } else {
    }
  }

  glEnable(GL_NORMALIZE);
}

void SpecialKeys(int key, int x, int y) {
  if (key == GLUT_KEY_UP)
    xRot -= 5.0f;

  if (key == GLUT_KEY_DOWN)
    xRot += 5.0f;

  if (key == GLUT_KEY_LEFT)
    yRot -= 5.0f;

  if (key == GLUT_KEY_RIGHT)
    yRot += 5.0f;

  if (key > 356.0f)
    xRot = 0.0f;

  if (key < -1.0f)
    xRot = 355.0f;

  if (key > 356.0f)
    yRot = 0.0f;

  if (key < -1.0f)
    yRot = 355.0f;

  // Refresh the Window
  glutPostRedisplay();
}

void ChangeSize(int w, int h) {
  GLfloat fAspect;

  // Prevent a divide by zero
  if (h == 0)
    h = 1;

  // Set Viewport to window dimensions
  glViewport(0, 0, w, h);

  // Reset coordinate system
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  fAspect = (GLfloat)w / (GLfloat)h;
  gluPerspective(60.0f, fAspect, 200.0, 500.0);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  // Move out Z axis so we can see everything
  glTranslatef(0.0f, 0.0f, -400.0f);
  glLightfv(GL_LIGHT0, GL_POSITION, lightPos[0].data());
  glLightfv(GL_LIGHT1, GL_POSITION, lightPos[1].data());
  glLightfv(GL_LIGHT2, GL_POSITION, lightPos[2].data());
  glLightfv(GL_LIGHT3, GL_POSITION, lightPos[3].data());
}

int main(int argc, char *argv[]) {
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  glutInitWindowSize(800, 600);
  glutCreateWindow("Shadow");
  glutReshapeFunc(ChangeSize);
  glutSpecialFunc(SpecialKeys);
  glutKeyboardFunc(KeyBoradFunc);
  glutDisplayFunc(RenderScene);
  SetupRC();
  glutMainLoop();

  return 0;
}
