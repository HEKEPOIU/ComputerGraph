#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string>
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include <freeglut.h>

// Lighting data
GLfloat lightAmbient[] = {0.2f, 0.2f, 0.2f, 1.0f};
GLfloat lightDiffuse[] = {0.7f, 0.7f, 0.7f, 1.0f};
GLfloat lightSpecular[] = {0.9f, 0.9f, 0.9f};
GLfloat materialColor[] = {0.8f, 0.0f, 0.0f};
GLfloat vLightPos[] = {-80.0f, 120.0f, 100.0f, 0.0f};
GLfloat ground[3][3] = {
    {0.0f, -25.0f, 0.0f}, {10.0f, -25.0f, 0.0f}, {10.0f, -25.0f, -10.0f}};

GLuint textures[4];

int nStep = 0;

void MyKeyboard(unsigned char key, int x, int y) {
  switch (key) {
  case 'r':
    if (nStep < 4)
      nStep++;
    else
      nStep = 0;
    break;
  default:
    break;
  }

  glutPostRedisplay();
}

// Called to draw scene
void RenderScene(void) {
  // Clear the window with current clearing color
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glShadeModel(GL_SMOOTH);
  glEnable(GL_NORMALIZE);

  glPushMatrix();

  // Draw plane that the cube rests on
  glDisable(GL_LIGHTING);
  glColor3ub(255, 255, 255);
  if (nStep >= 1) {
    glEnable(GL_TEXTURE_2D); // 啟動openGL的2D材質填充模式

    // 將textures[0]中所儲存的材質貼在四邊形上
    glBindTexture(GL_TEXTURE_2D, textures[0]);
    glBegin(GL_QUADS);
    glTexCoord2f(0.0f, 0.0f);
    glVertex3f(-100.0f, -25.3f, -100.0f);
    glTexCoord2f(0.0f, 1.0f);
    glVertex3f(-100.0f, -25.3f, 100.0f);
    glTexCoord2f(1.0f, 1.0f);
    glVertex3f(100.0f, -25.3f, 100.0f);
    glTexCoord2f(1.0f, 0.0f);
    glVertex3f(100.0f, -25.3f, -100.0f);
    glEnd();

    glDisable(GL_TEXTURE_2D);
  } else {
    // glColor3f(0.0f, 0.0f, 0.90f); // Blue
    glBegin(GL_QUADS);
    glVertex3f(-100.0f, -25.3f, -100.0f);
    glVertex3f(-100.0f, -25.3f, 100.0f);
    glVertex3f(100.0f, -25.3f, 100.0f);
    glVertex3f(100.0f, -25.3f, -100.0f);
    glEnd();
  }

  // Set drawing color to Red
  glColor3f(1.0f, 0.0f, 0.0f);

  // Move the cube slightly forward and to the left
  glTranslatef(-10.0f, 0.0f, 10.0f);

  glColor3ub(255, 255, 255);

  if (nStep >= 2) {
    glEnable(GL_TEXTURE_2D);
  }
  // Front Face (before rotation)
  glBindTexture(GL_TEXTURE_2D, textures[1]);
  glBegin(GL_QUADS);
  glTexCoord2f(1.0f, 1.0f);
  glVertex3f(25.0f, 25.0f, 25.0f);
  glTexCoord2f(1.0f, 0.0f);
  glVertex3f(25.0f, -25.0f, 25.0f);
  glTexCoord2f(0.0f, 0.0f);
  glVertex3f(-25.0f, -25.0f, 25.0f);
  glTexCoord2f(0.0f, 1.0f);
  glVertex3f(-25.0f, 25.0f, 25.0f);
  glEnd();

  glDisable(GL_TEXTURE_2D);

  if (nStep >= 3) {
    glEnable(GL_TEXTURE_2D);
  }
  // Top of cube
  glBindTexture(GL_TEXTURE_2D, textures[2]);
  glBegin(GL_QUADS);
  // Front Face
  glTexCoord2f(0.0f, 0.0f);
  glVertex3f(25.0f, 25.0f, 25.0f);
  glTexCoord2f(1.0f, 0.0f);
  glVertex3f(25.0f, 25.0f, -25.0f);
  glTexCoord2f(1.0f, 1.0f);
  glVertex3f(-25.0f, 25.0f, -25.0f);
  glTexCoord2f(0.0f, 1.0f);
  glVertex3f(-25.0f, 25.0f, 25.0f);
  glEnd();

  glDisable(GL_TEXTURE_2D);

  if (nStep >= 4) {
    glEnable(GL_TEXTURE_2D);
  }
  // Last two segments for effect
  glBindTexture(GL_TEXTURE_2D, textures[3]);
  glBegin(GL_QUADS);
  glTexCoord2f(1.0f, 1.0f);
  glVertex3f(25.0f, 25.0f, -25.0f);
  glTexCoord2f(1.0f, 0.0f);
  glVertex3f(25.0f, -25.0f, -25.0f);
  glTexCoord2f(0.0f, 0.0f);
  glVertex3f(25.0f, -25.0f, 25.0f);
  glTexCoord2f(0.0f, 1.0f);
  glVertex3f(25.0f, 25.0f, 25.0f);
  glEnd();

  glDisable(GL_TEXTURE_2D);

  glTranslatef(-10.0f, 0.0f, 10.0f);

  glPopMatrix();

  // Flush drawing commands
  glutSwapBuffers();
}

void LoadTexture(const std::string &path, int index) {
  int width4, height4, channels4;
  unsigned char *image4 =
      stbi_load(path.c_str(), &width4, &height4, &channels4, 0);
  if (!image4) {
    std::cout << "image4 empty\n";
  } else {
    // 將讀取進來的圖片檔案當作材質存進textures中
    glGenTextures(1, &textures[index]);
    glBindTexture(GL_TEXTURE_2D, textures[index]);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);

    GLenum format;
    if (channels4 == 1)
      format = GL_RED;
    else if (channels4 == 3)
      format = GL_RGB;
    else if (channels4 == 4)
      format = GL_RGBA;
    else
      format = 0; // 处理不支持的通道数

    if (format) {
      glTexImage2D(GL_TEXTURE_2D, 0, format, width4, height4, 0, format,
                   GL_UNSIGNED_BYTE, image4);
    } else {
      std::cout << "Unsupported image format\n";
    }

    stbi_image_free(image4);
  }
}

// This function does any needed initialization on the rendering
// context.
void SetupRC() {
  GLbyte *pBytes;
  GLint nWidth, nHeight, nComponents;
  GLenum format;

  // Black background
  glClearColor(0.0f, 0.0f, 0.0f, 1.0f);

  glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE,
            GL_MODULATE); // 設定openGL材質紋理的參數和材質的組合模式
  glGenTextures(
      4, textures); // 註冊一個大小為4的陣列讓openGL儲存材質，名稱為textures

  LoadTexture(RESOURCE_DIR "//floor.jpg", 0);
  LoadTexture(RESOURCE_DIR "/Block4.jpg", 1);
  LoadTexture(RESOURCE_DIR "/Block5.jpg", 2);
  LoadTexture(RESOURCE_DIR "/Block6.jpg", 3);

  // 下面要設定剩下的三種材質
  // ...
}

void ChangeSize(int w, int h) {
  // Calculate new clipping volume
  GLfloat windowWidth = 100.f;
  GLfloat windowHeight = 100.f;

  // Set the viewport to be the entire window
  glViewport(0, 0, w, h);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  // Set the clipping volume
  glOrtho(-100.0f, windowWidth, -100.0f, windowHeight, -200.0f, 200.0f);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glLightfv(GL_LIGHT0, GL_POSITION, vLightPos);

  glRotatef(30.0f, 1.0f, 0.0f, 0.0f);
  glRotatef(330.0f, 0.0f, 1.0f, 0.0f);
}

int main(int argc, char *argv[]) {
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
  glutInitWindowSize(800, 600);
  glutCreateWindow("Textures");
  SetupRC();
  glutReshapeFunc(ChangeSize);
  glutDisplayFunc(RenderScene);
  glutKeyboardFunc(MyKeyboard);

  glutMainLoop();
  glDeleteTextures(4, textures);
  return 0;
}
