#include "CGMath.hpp"
#include "Transform.hpp"
#include <GL/gl.h>
#include <iostream>
#include <ostream>
#define STB_IMAGE_IMPLEMENTATION
#include "Dp.hpp"
#include "DrawableObject.hpp"
#include "OBJLoader.hpp"
#include "Texture.hpp"
/*** freeglut***/
#include <freeglut.h>
#include <freeglut_std.h>
#include <memory>
#include <string>
#include <vector>
GLFrame frameCamera;

// Light and material Data
GLfloat fLightPos[4] = {-100.0f, 100.0f, 50.0f, 1.0f}; // Point source
GLfloat fNoLight[] = {0.f, 0.f, 0.f, 0.f};
GLfloat fLowLight[] = {0.25f, 0.25f, 0.25f, 1.0f};
GLfloat fBrightLight[] = {1.0f, 1.0f, 1.0f, 1.0f};

// This Plane point.
M3DVector3f vPoints[3] = {
    {0.0f, -0.4f, 0.0f}, {10.0f, -0.4f, 0.0f}, {5.0f, -0.4f, -5.0f}};

M3DMatrix44f mShadowMatrix;

std::shared_ptr<Texture> groundTexture;
#define TRUNKWOOD 0
#define CITY 1
#define UFO 2
#define HORSE 3
#define PENGUIN 4
std::unique_ptr<OBJLoader> _objLoader = std::make_unique<OBJLoader>();
std::vector<std::string> _modelPaths = {"trunk_wood", "city", "ufo", "horse",
                                        "penguin"};
Vec3 preMousePos = Vec3(0, 0);
std::shared_ptr<DrawableObject> _skybox;
std::vector<std::shared_ptr<DrawableObject>> _models;
float horse_angle = 0;
float horse_float_time = 0;
float lastX = 3;
float lastZ = 0.0f;
float lastDirectionAngle = 0.0f;
void SetupRC() {

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

  // Mostly use material tracking
  glEnable(GL_COLOR_MATERIAL);
  glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
  glMaterialfv(GL_FRONT, GL_SPECULAR, fBrightLight);
  glMateriali(GL_FRONT, GL_SHININESS, 128);

  // Set up texture maps, only need load once,
  glEnable(GL_TEXTURE_2D);
  glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

  groundTexture = std::make_shared<Texture>(RESOURCE_DIR "/texture/rock.tga");
  for (auto modelpath : _modelPaths) {
    auto model = _objLoader->get_OBJ(RESOURCE_DIR "/" + modelpath + ".obj");
    _models.push_back(model);
    auto texture = std::make_shared<Texture>(RESOURCE_DIR "/texture/" +
                                             modelpath + "_diff.tga");
    model->SetTexture(texture);
  }
  _models[TRUNKWOOD]->set_transform_to_target({0, 0, 0},
                                              {-1, 1, -1, 1, -100, 100});
  _models[TRUNKWOOD]->get_transform()->set_scale({3, 3, 3});
  _models[TRUNKWOOD]->get_transform()->set_location({-.5, -0.5, -5});
  _models[CITY]->set_transform_to_target({0, 0, 0},
                                         {-15, 15, -15, 15, -500, 500});
  _models[CITY]->get_transform()->set_location({0, -0.5, -20});
  _models[UFO]->set_transform_to_target({0, 0, 0}, {-1, 1, -1, 1, -100, 100});
  _models[UFO]->get_transform()->set_location({0, 5, -15});
  _models[HORSE]->set_transform_to_target({0, 0, 0}, {-1, 1, -1, 1, -100, 100});
  _skybox = _objLoader->get_OBJ(RESOURCE_DIR "/"
                                             "skybox"
                                             ".obj");
  auto texture = std::make_shared<Texture>(RESOURCE_DIR "/texture/"
                                                        "skybox"
                                                        "_diff.tga");
  _skybox->SetTexture(texture);
  _skybox->get_transform()->set_scale({3, 3, 3});
}

////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////
// Draw the ground as a series of triangle strips
void DrawGround(void) {
  GLfloat fExtent = 100.0f;
  GLfloat fStep = 1.0f;
  GLfloat y = -0.4f;
  GLfloat iStrip, iRun;
  GLfloat s = 0.0f;
  GLfloat t = 0.0f;
  GLfloat texStep = 1.0f / (fExtent * .075f);

  groundTexture->ApplyTexture();
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

///////////////////////////////////////////////////////////////////////
// Draw random inhabitants and the rotating torus/sphere duo
void DrawInhabitants(GLint nShadow) {
  if (nShadow == 0) {
    glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
  } else
    glColor4f(0.00f, 0.00f, 0.00f, .6f); // Shadow color

  if (nShadow == 0) {
    // Torus alone will be specular
    glMaterialfv(GL_FRONT, GL_SPECULAR, fBrightLight);
  }

  //   --------------------------
  //   Draw the Object below.
  glPushMatrix();
  { _models[TRUNKWOOD]->draw(DrawableObject::DrawType::FACES); }
  glPopMatrix();

  glPushMatrix();
  float horseX = cos(horse_angle);
  glScalef(1, horseX + 1, 1);
  { _models[PENGUIN]->draw(DrawableObject::DrawType::FACES); }
  glPopMatrix();

  glPushMatrix();
  { _models[CITY]->draw(DrawableObject::DrawType::FACES); }
  glPopMatrix();
  glPushMatrix();
  {
    float horseX = 3 * cos(horse_angle);
    float horseZ = 3 * sin(horse_angle);
    float horseY = 3 * cos(horse_float_time);
    // 計算馬的朝向角度
    float deltaX = horseX - lastX;
    float deltaZ = horseZ - lastZ;
    float directionAngle = atan2(deltaZ, deltaX) * 180 / M_PI; // 弧度轉度數
    float forwardX = cos(horse_angle);
    float forwardZ = sin(horse_angle);
    float facingAngle = atan2(forwardZ, forwardX) * 180 / M_PI;

    float smoothingFactor = 0.1f; // 平滑因子
    if (fabs(directionAngle - lastDirectionAngle) < 180) {
      directionAngle = lastDirectionAngle * (1.0f - smoothingFactor) +
                       directionAngle * smoothingFactor;
    }

    lastX = horseX;
    lastZ = horseZ;
    lastDirectionAngle = directionAngle;
    glTranslatef(0, horseY, 0);
    _models[UFO]->draw(DrawableObject::DrawType::FACES);
    glPushMatrix();
    {
      glTranslatef(horseX, horseY, horseZ);
      glRotatef(facingAngle, 0, 1, 0);
      _models[HORSE]->draw(DrawableObject::DrawType::FACES);
    }
    glPopMatrix();
  }
  glPopMatrix();
}

// Called to draw scene
void RenderScene(void) {
  // Clear the window with current clearing color
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

  glPushMatrix();
  frameCamera.ApplyCameraTransform();

  _skybox->draw(DrawableObject::DrawType::FACES);
  // Position light before any other transformations
  glLightfv(GL_LIGHT0, GL_POSITION, fLightPos);

  // Draw the ground
  glColor3f(1.0f, 1.0f, 1.0f);
  DrawGround();

  // Calculate shadow matrix
  M3DVector4f pPlane;
  m3dGetPlaneEquation(pPlane, vPoints[0], vPoints[1], vPoints[2]);
  m3dMakePlanarShadowMatrix(mShadowMatrix, pPlane, fLightPos);

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
void SpecialKeys(unsigned char key, int x, int y) {
  if (key == 'w')
    frameCamera.MoveForward(0.1f);

  if (key == 's')
    frameCamera.MoveForward(-0.1f);

  // Refresh the Window
}

///////////////////////////////////////////////////////////
// Called by GLUT library when idle (window not being
// resized or moved)
void TimerFunction(int value) {
  _models[TRUNKWOOD]->get_transform()->modify_rotation({0, 1.f, 0});
  horse_angle += .025f;
  horse_float_time += 0.05f;
  glutTimerFunc(3, TimerFunction, 1);
  glutPostRedisplay();
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
  gluPerspective(35.0f, fAspect, 1.0f, 1000.0f);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
}
void MouseFunc(int x, int y) {
  Vec3 nowMousePos = {float(y), float(x), 0};
  Vec3 delta = nowMousePos - preMousePos;
  frameCamera.RotateLocal(0.01f, -delta.x, delta.y, delta.z);
  preMousePos = {float(y), float(x), 0};
}
int main(int argc, char *argv[]) {
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | GLUT_STENCIL);
  glutInitWindowSize(800, 600);
  glutCreateWindow("ComputerGraph Final Project");
  glutReshapeFunc(ChangeSize);
  glutDisplayFunc(RenderScene);
  glutKeyboardFunc(SpecialKeys);
  glutPassiveMotionFunc(MouseFunc);
  SetupRC();
  glutTimerFunc(33, TimerFunction, 1);

  frameCamera.SetOrigin(0, 2, 0);
  glutMainLoop();

  return 0;
}
