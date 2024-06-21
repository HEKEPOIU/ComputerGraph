#include "Dp.hpp"
#include "DrawableObject.hpp"
#include "ImageLoader.hpp"
/*** freeglut***/
#include <algorithm>
#include <freeglut.h>
#include <freeglut_std.h>
#include <iostream>
#include <memory>
#include <ostream>

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

const char *szTextureFiles[] = {RESOURCE_DIR "/texture/grass.tga",
                                RESOURCE_DIR "/texture/wood.tga",
                                RESOURCE_DIR "/texture/orb.tga"};
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

//
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

    glBindTexture(GL_TEXTURE_2D, textureObjects[i]);
    GLenum eFormat;

    auto imageData = std::make_unique<ImageLoader>(szTextureFiles[i]);
    if (imageData->channels == 1) {
      eFormat = GL_RED;
    } else if (imageData->channels == 3) {
      eFormat = GL_RGB;
    } else if (imageData->channels == 4) {
      eFormat = GL_RGBA;
    } else {
      std::cerr << "Unsupported image format!" << std::endl;
      return;
    }

    gluBuild2DMipmaps(GL_TEXTURE_2D, imageData->channels, imageData->width,
                      imageData->height, eFormat, GL_UNSIGNED_BYTE,
                      imageData->data.get());

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
  glBindTexture(GL_TEXTURE_2D, textureObjects[SPHERE_TEXTURE]);

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
