
#include <X11/Xlib.h>
#include <array>
#include <cmath>
#include <iostream>
#include <math.h>
#include <ostream>
#include <queue>
#include <random>
#include <stdbool.h>
#include <stdio.h>
/*** freeglut***/
#include "freeglut_std.h"
#include <freeglut.h>
#include <utility>
#include <vector>

#define PI 3.141592653589793

void ChangeSize(int, int);
void RenderScene(void);
void Timer(int);
void MenuCallback(int);
void OnKeyBoardPress(unsigned char, int, int);
void TranslateMatrix(GLfloat, GLfloat, GLfloat);
void RotateMatrix(float angle, float x, float y, float z);
void ScaleMatrix(float x, float y, float z);
void MousePress(int button, int state, int x, int y);
void ChangeGridPointSize(std::vector<std::vector<int>> &grid_point, int x,
                         int y);
void DrawOutLine(int x, int y);
void DrawGrid(std::vector<std::vector<int>> &grid_point);
void DrawGridLine(int x0, int y0, int x1, int y1);
void GridSpaceToGridVectorSpace(int gridX, int gridY, int &x, int &y);
void ScreenPosToGridPoint(float x, float y, int &gridX, int &gridY);
void HalfSpaceFillCell(const std::vector<std::array<float, 2>> &mousePoints,
                       const std::vector<std::array<float, 2>> &,
                       const std::array<float, 3> &color);
void CrowFillCell(const std::vector<std::array<float, 2>> &gridPoint);

std::array<int, 2> windowSize{800, 800};
std::vector<std::vector<int>> _gridPoint{};
std::vector<std::pair<bool, std::array<float, 3>>> _gridDraw{};
std::queue<std::pair<int, std::array<float, 3>>> _waitToDraw{};
std::vector<std::array<float, 2>> mousePoints;
int size_x = 21; // -10 ~ 10
int size_y = 21;

int main(int argc, char **argv) {
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
  glutInitWindowSize(800, 800);
  glutInitWindowPosition(600, 80);
  glutCreateWindow("Simple Triangle");
  glutCreateMenu(MenuCallback);
  glutAddMenuEntry("10", 1);
  glutAddMenuEntry("15", 2);
  glutAddMenuEntry("20", 3);
  ChangeGridPointSize(_gridPoint, 21, 21);
  glutAttachMenu(GLUT_RIGHT_BUTTON);
  glutReshapeFunc(ChangeSize);
  glutKeyboardFunc(OnKeyBoardPress);
  glutMouseFunc(MousePress);
  glutDisplayFunc(RenderScene);
  glutTimerFunc(0, Timer, 0);
  glutMainLoop(); // http://www.programmer-club.com.tw/ShowSameTitleN/opengl/2288.html
  return 0;
}
void ChangeSize(int w, int h) {
  printf("Window Size= %d X %d\n", w, h);
  glViewport(0, 0, w, h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(-400, 400, -400, 400, -50, 50);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
}
void RenderScene(void) {
  glClearColor(0, 0, 0, 1.0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(0, 0, 10.0f, 0, 0, 0, 0, 1, 0);
  glEnable(GL_DEPTH_TEST);
  glViewport(0, 0, windowSize[0], windowSize[1]);

  DrawGrid(_gridPoint);

  glutSwapBuffers();
}

void MenuCallback(int value) {
  // TODO: change grid size.
  switch (value) {
  case 1:
    ChangeGridPointSize(_gridPoint, 21, 21);
    mousePoints.clear();
    break;
  case 2:
    ChangeGridPointSize(_gridPoint, 31, 31);
    mousePoints.clear();
    break;
  case 3:
    ChangeGridPointSize(_gridPoint, 41, 41);
    mousePoints.clear();
    break;
  }
  glutPostRedisplay();
}

void print4mat(GLfloat array[16]) {
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      printf("%-3f ", array[4 * j + i]);
    }
    printf("\n");
  }
}

void TranslateMatrix(GLfloat x, GLfloat y, GLfloat z) {
  GLfloat rotMatrix[16];
  glGetFloatv(GL_MODELVIEW_MATRIX, rotMatrix);

  rotMatrix[12] = x;
  rotMatrix[13] = y;
  rotMatrix[14] = z;
  glMultMatrixf(rotMatrix);
}

void RotateMatrix(float angle, float x, float y, float z) {
  GLfloat rotMatrix[16];
  float angleDeg2Rad = angle * (M_PI / 180.0f);

  float cosAngle = cos(angleDeg2Rad);
  float sinAngle = sin(angleDeg2Rad);

  float len = sqrt(x * x + y * y + z * z);

  if (len != 0) {
    x /= len;
    y /= len;
    z /= len;
  }

  float omc = 1.0f - cosAngle;

  rotMatrix[0] = x * x * omc + cosAngle;
  rotMatrix[1] = y * x * omc + z * sinAngle;
  rotMatrix[2] = x * z * omc - y * sinAngle;
  rotMatrix[3] = 0.0f;

  rotMatrix[4] = x * y * omc - z * sinAngle;
  rotMatrix[5] = y * y * omc + cosAngle;
  rotMatrix[6] = y * z * omc + x * sinAngle;
  rotMatrix[7] = 0.0f;

  rotMatrix[8] = x * z * omc + y * sinAngle;
  rotMatrix[9] = y * z * omc - x * sinAngle;
  rotMatrix[10] = z * z * omc + cosAngle;
  rotMatrix[11] = 0.0f;

  rotMatrix[12] = 0.0f;
  rotMatrix[13] = 0.0f;
  rotMatrix[14] = 0.0f;
  rotMatrix[15] = 1.0f;
  glMultMatrixf(rotMatrix);
}

void ScaleMatrix(float x, float y, float z) {
  GLfloat scaleMatrix[16]{};
  scaleMatrix[0] = x;
  scaleMatrix[5] = y;
  scaleMatrix[10] = z;
  scaleMatrix[15] = 1.0f;
  glMultMatrixf(scaleMatrix);
}

void ChangeGridPointSize(std::vector<std::vector<int>> &grid_point, int totalX,
                         int totalY) {
  size_x = totalX;
  size_y = totalY;
  grid_point.clear();
  _gridDraw.clear();
  _waitToDraw = {};

  int halfX = totalX / 2;
  int halfY = totalY / 2;
  for (int x = -halfX; x <= halfX; x++) {
    for (int y = -halfY; y <= halfY; y++) {
      grid_point.push_back({x, y});
      _gridDraw.push_back({false, {0.0f, 0.0f, 0.0f}});
    }
  }
}

void DrawOutLine(int x, int y) {
  int w = windowSize[0];
  int h = windowSize[1];
  int sizeX = w / size_x;
  int sizeY = h / size_y;
  glColor3f(1.0f, 1.0f, 1.0f);

  glBegin(GL_LINES);

  glVertex3f(x - sizeX / 2, y + sizeY / 2, 0);
  glVertex3f(x - sizeX / 2, y - sizeY / 2, 0);

  glVertex3f(x - sizeX / 2, y - sizeY / 2, 0);
  glVertex3f(x + sizeX / 2, y - sizeY / 2, 0);

  glVertex3f(x + sizeX / 2, y - sizeY / 2, 0);
  glVertex3f(x + sizeX / 2, y + sizeY / 2, 0);

  glVertex3f(x + sizeX / 2, y + sizeY / 2, 0);
  glVertex3f(x - sizeX / 2, y + sizeY / 2, 0);

  glEnd();
}

void DrawFillCell(int x, int y, std::array<float, 3> color) {
  int w = windowSize[0];
  int h = windowSize[1];
  int sizeX = w / size_x;
  int sizeY = h / size_y;

  glColor3fv(color.data());
  glBegin(GL_QUADS);

  glVertex3f(x - sizeX / 2, y + sizeY / 2, 0);
  glVertex3f(x - sizeX / 2, y - sizeY / 2, 0);

  glVertex3f(x + sizeX / 2, y - sizeY / 2, 0);

  glVertex3f(x + sizeX / 2, y + sizeY / 2, 0);

  glEnd();
}

void DrawGrid(std::vector<std::vector<int>> &grid_point) {
  int spacingX = windowSize[0] / (size_x);
  int spacingY = windowSize[1] / (size_y);
  for (int i = 0; i < grid_point.size(); i++) {
    if (_gridDraw[i].first == true) {
      DrawFillCell(grid_point[i][0] * spacingX, grid_point[i][1] * spacingY,
                   _gridDraw[i].second);
    } else {
      DrawOutLine(grid_point[i][0] * spacingX, grid_point[i][1] * spacingY);
    }
  }
}

void Timer(int value) {
  if (!_waitToDraw.empty()) {
    auto value = _waitToDraw.front();
    _gridDraw[value.first].first = true;
    _gridDraw[value.first].second = value.second;
    _waitToDraw.pop();
  }
  glutPostRedisplay();
  glutTimerFunc(50, Timer, 0);
}

void OnKeyBoardPress(unsigned char key, int x, int y) {
  switch (key) {
  case ' ':
    mousePoints.clear();
    break;
  default:
    break;
  }
  glutPostRedisplay();
}
void ScreenPosToGridPoint(float x, float y, int &gridX, int &gridY) {
  int w = windowSize[0];
  int h = windowSize[1];
  int sizeX = w / (size_x);
  int sizeY = h / (size_y);
  gridX = std::round(x / sizeX);
  gridY = std::round(y / sizeY);
}

void MousePress(int button, int state, int x, int y) {
  if (state == GLUT_DOWN) {
    float mousePosX = (2 * ((float)x / 800) - 1) * (400);
    float mousePosY = (-2 * ((float)y / 800) + 1) * (400);

    if (mousePoints.size() >= 3) {
      // mousePoints.clear();
      ChangeGridPointSize(_gridPoint, size_x, size_y);
      // glutPostRedisplay();
      // return;
    }

    mousePoints.push_back({mousePosX, mousePosY});

    if (mousePoints.size() >= 2) {
      std::vector<std::array<float, 2>> mpoint;
      // glColor3f(1.0f, 0.f, 0.f);
      for (int i = 0; i < mousePoints.size(); i++) {
        // Debug Use
        // glVertex3f(mousePoints[i][0], mousePoints[i][1], 0);
        // glVertex3f(mousePoints[(i + 1) % mousePoints.size()][0],
        //            mousePoints[(i + 1) % mousePoints.size()][1], 0);

        int gridX0, gridY0, gridX1, gridY1;
        ScreenPosToGridPoint(mousePoints[i][0], mousePoints[i][1], gridX0,
                             gridY0);
        ScreenPosToGridPoint(mousePoints[(i + 1) % mousePoints.size()][0],
                             mousePoints[(i + 1) % mousePoints.size()][1],
                             gridX1, gridY1);
        GridSpaceToGridVectorSpace(gridX0, gridY0, gridX0, gridY0);
        GridSpaceToGridVectorSpace(gridX1, gridY1, gridX1, gridY1);

        DrawGridLine(gridX0, gridY0, gridX1, gridY1);
        mpoint.push_back({(float)gridX0, (float)gridY0});
      }

      if (mousePoints.size() >= 3) {
        CrowFillCell(mpoint);
      }
    }

    glutPostRedisplay();
  }
  // gluInvertMatrix();
}

void GridSpaceToGridVectorSpace(int gridX, int gridY, int &x, int &y) {
  x = gridX + size_x / 2;
  y = gridY + size_y / 2;
}

float RandomF(float minY, float maxY) {
  std::random_device rd;  // 使用随机设备生成种子
  std::mt19937 gen(rd()); // 使用梅森旋转算法的随机数生成器

  // 定义0到1之间的浮点数分布
  std::uniform_real_distribution<> dis(minY, maxY);

  // 生成并输出一个随机浮点数
  return dis(gen);
}

void DrawGridLine(int x0, int y0, int x1, int y1) {
  int dx = abs(x1 - x0);
  int dy = abs(y1 - y0);
  int sx = x0 < x1 ? 1 : -1;
  int sy = y0 < y1 ? 1 : -1;
  int err = dx - dy;

  float len = sqrtf(powf(abs(x1 - x0), 2) + powf(abs(y1 - y0), 2));

  std::array<std::array<float, 3>, 2> pColor;
  pColor[0] = {RandomF(0, 1), RandomF(0, 1), RandomF(0, 1)};
  pColor[1] = {RandomF(0, 1), RandomF(0, 1), RandomF(0, 1)};
  if (_gridDraw[y0 + x0 * size_y].first) {
    pColor[0] = _gridDraw[y0 + x0 * size_y].second;
  }
  if (_gridDraw[y1 + x1 * size_y].first) {
    pColor[1] = _gridDraw[y1 + x1 * size_y].second;
  }

  _gridDraw[y0 + x0 * size_y].first = true;
  _gridDraw[y0 + x0 * size_y].second = pColor[0];

  // std::cout << "Start Points: " << x0 - size_x / 2 << " " << y0 - size_x / 2
  //           << std::endl;
  while (x0 != x1 || y0 != y1) {
    _gridDraw[y0 + x0 * size_y].first = true;
    // std::cout << x0 - size_x / 2 << " " << y0 - size_x / 2 << std::endl;
    int e2 = 2 * err;
    if (e2 > -dy) {
      err -= dy;
      x0 += sx;
    }
    if (e2 < dx) {
      err += dx;
      y0 += sy;
    }
    float lenToSecond = sqrtf(powf(abs(x1 - x0), 2) + powf(abs(y1 - y0), 2));
    float t = (len - lenToSecond) / len;
    _gridDraw[y0 + x0 * size_y].second = {
        pColor[0][0] * (1 - t) + pColor[1][0] * (t),
        pColor[0][1] * (1 - t) + pColor[1][1] * (t),
        pColor[0][2] * (1 - t) + pColor[1][2] * (t)};
  }
  _gridDraw[y1 + x1 * size_y].first = true;
  _gridDraw[y1 + x1 * size_y].second = pColor[1];
}

std::array<float, 4>
GetBbox(const std::vector<std::array<float, 2>> &FillPoints) {

  float xMin = FillPoints[0][0];
  float xMax = FillPoints[0][0];
  float yMin = FillPoints[0][1];
  float yMax = FillPoints[0][1];

  for (auto point : FillPoints) {
    if (point[0] < xMin) {
      xMin = point[0];
    }
    if (point[0] > xMax) {
      xMax = point[0];
    }
    if (point[1] < yMin) {
      yMin = point[1];
    }
    if (point[1] > yMax) {
      yMax = point[1];
    }
  }
  return {xMin, xMax, yMin, yMax};
}
std::array<int, 4> BboxToCellBbox(const std::array<float, 4> &bbox) {
  std::array<int, 4> cellBbox;
  ScreenPosToGridPoint(bbox[0], bbox[2], cellBbox[0], cellBbox[2]);
  ScreenPosToGridPoint(bbox[1], bbox[3], cellBbox[1], cellBbox[3]);
  GridSpaceToGridVectorSpace(cellBbox[0], cellBbox[2], cellBbox[0],
                             cellBbox[2]);
  GridSpaceToGridVectorSpace(cellBbox[1], cellBbox[3], cellBbox[1],
                             cellBbox[3]);
  return cellBbox;
}
std::vector<std::array<int, 2>> GridPointToGridVectorPoint(
    const std::vector<std::array<float, 2>> &ScreenPoints) {

  std::vector<std::array<int, 2>> grid_vector_point;
  for (auto point : ScreenPoints) {
    std::array<int, 2> girdPoint;
    ScreenPosToGridPoint(point[0], point[1], girdPoint[0], girdPoint[1]);
    GridSpaceToGridVectorSpace(girdPoint[0], girdPoint[1], girdPoint[0],
                               girdPoint[1]);
    grid_vector_point.push_back(girdPoint);
    std::cout << girdPoint[0] << " " << girdPoint[1] << std::endl;
  }
  return grid_vector_point;
}

void findEdge(const std::vector<std::array<float, 2>> &vList, int &i, int n,
              int &y, std::array<float, 2> &d, std::array<float, 2> &e,
              int &rem, int lr) {

  i = (i + lr) % n;
  rem--;

  e = vList[i];
  int prev = ((i - lr) + n) % n;
  d[0] = (vList[prev][0] - e[0]) / (vList[prev][1] - e[1]);
  d[1] = vList[prev][1] - e[1];
}

void CrowFillCell(const std::vector<std::array<float, 2>> &gridPoint) {

  int li, ri;                  // left & right upper endpoint indices
  float ly, ry;                // left & right upper endpoint y values
  std::array<float, 2> dl, dr; // current left edge and delta
  std::array<float, 2> l, r;   // current left and right edge
  int rem;                     // number of remaining vertices
  int y;                       // current scanline
  float lx, rx;
  int n = gridPoint.size();
  std::vector<std::pair<int, std::array<float, 3>>> originPointColor;

  int minYindex = 0;
  int minY = gridPoint[0][1];
  int minX = gridPoint[0][0];
  for (int i = 0; i < gridPoint.size(); i++) {
    originPointColor.push_back(
        {i, _gridDraw[int(round(gridPoint[i][1])) +
                      int(round(gridPoint[i][0])) * size_y]
                .second});
    if (minY > gridPoint[i][1]) {
      minX = gridPoint[i][0];
      minY = gridPoint[i][1];
      minYindex = i;
    }
  }

  li = ri = minYindex;
  ly = ry = y = round(gridPoint[minYindex][1]);
  lx = rx = round(gridPoint[minYindex][0]);
  rem = n;

  auto drawSpan = [](int y, int xl, int xr) {
    int len = xr - xl;
    std::array<std::array<float, 3>, 2> pColor;

    if (_gridDraw[y + xl * size_y].first == true) {
      pColor[0] = _gridDraw[y + xl * size_y].second;
    } else if (_gridDraw[y + (xl + 1) * size_y].first == true) {
      pColor[0] = _gridDraw[y + (xl + 1) * size_y].second;
    } else if (_gridDraw[y + (xl - 1) * size_y].first == true) {
      pColor[0] = _gridDraw[y + (xl - 1) * size_y].second;
    }

    if (_gridDraw[y + xr * size_y].first == true) {
      pColor[1] = _gridDraw[y + xr * size_y].second;
    } else if (_gridDraw[y + (xr + 1) * size_y].first == true) {
      pColor[1] = _gridDraw[y + (xr + 1) * size_y].second;
    } else if (_gridDraw[y + (xr - 1) * size_y].first == true) {
      pColor[1] = _gridDraw[y + (xr - 1) * size_y].second;
    }

    // std::cout << "L_Point: " << pColor[0][0] << "," << pColor[0][1] << ","
    //           << pColor[0][2] << std::endl;
    // std::cout << "R_Point: " << pColor[1][0] << "," << pColor[1][1] << ","
    //           << pColor[1][2] << std::endl;
    std::cout << "len" << len << std::endl;
    for (int x = xl; x <= xr; ++x) {
      if (_gridDraw[y + x * size_y].first) {
        continue;
      }
      float lenToL = x - xl;
      float t = lenToL / len;
      _waitToDraw.push({y + x * size_y,
                        {pColor[0][0] * (1 - t) + pColor[1][0] * (t),

                         pColor[0][1] * (1 - t) + pColor[1][1] * (t),
                         pColor[0][2] * (1 - t) + pColor[1][2] * (t)}});
    }
  };

  findEdge(gridPoint, ri, n, y, dr, r, rem, -1);
  findEdge(gridPoint, li, n, y, dl, l, rem, 1);

  _waitToDraw.push({minY + minX * size_y, {1, 0, 0}});
  _waitToDraw.push(
      {int(round(gridPoint[ri][1])) + int(round(gridPoint[ri][0])) * size_y,
       {1, 0, 0}});
  _waitToDraw.push(
      {int(round(gridPoint[li][1])) + int(round(gridPoint[li][0])) * size_y,
       {1, 0, 0}});

  while (rem > 0) {

    bool isUpdateEdge = false;
    std::vector<int> currentRoundUpdateEdgePoint;
    // Find the appropriate left edge
    if (ceil(l[1]) <= y) {
      findEdge(gridPoint, li, n, y, dl, l, rem, 1);
      isUpdateEdge = true;
    }

    // Find the appropriate right edge
    if (ceil(r[1]) <= y) {
      findEdge(gridPoint, ri, n, y, dr, r, rem, -1);
      isUpdateEdge = true;
    }

    if (isUpdateEdge) {
      int rprev = ((ri + 1) + n) % n;

      currentRoundUpdateEdgePoint.push_back(rprev);
      currentRoundUpdateEdgePoint.push_back(ri);
      _waitToDraw.push({int(round(gridPoint[rprev][1])) +
                            int(round(gridPoint[rprev][0])) * size_y,
                        {1, 0, 0}});
      _waitToDraw.push(
          {int(round(gridPoint[ri][1])) + int(round(gridPoint[ri][0])) * size_y,
           {1, 0, 0}});

      int lprev = ((li - 1) + n) % n;

      currentRoundUpdateEdgePoint.push_back(lprev);
      currentRoundUpdateEdgePoint.push_back(li);

      _waitToDraw.push({int(round(gridPoint[lprev][1])) +
                            int(round(gridPoint[lprev][0])) * size_y,
                        {1, 0, 0}});
      _waitToDraw.push(
          {int(round(gridPoint[li][1])) + int(round(gridPoint[li][0])) * size_y,
           {1, 0, 0}});

      for (int i = 0; i < gridPoint.size(); i++) {
        bool isChange = false;
        for (auto changePoint : currentRoundUpdateEdgePoint) {
          if (changePoint == i) {
            isChange = true;
          }
        }
        if (!isChange) {
          std::cout << i << std::endl;
          _waitToDraw.push({int(round(gridPoint[i][1])) +
                                int(round(gridPoint[i][0])) * size_y,
                            originPointColor[i].second});
        }
      }
    }

    // While l & r span y (the current scanline)
    // std::cout << minYindex << std::endl;

    while (round(l[1]) > y && round(r[1]) > y) {
      // std::cout << r[0] << std::endl;
      // std::cout << "ly " << l[1] << " ry " << r[1] << std::endl;

      drawSpan(y, round(lx), round(rx));
      lx += dl[0];
      rx += dr[0];
      y++;
    }

    for (int i = 0; i < gridPoint.size(); i++) {
      _waitToDraw.push(
          {int(round(gridPoint[i][1])) + int(round(gridPoint[i][0])) * size_y,
           originPointColor[i].second});
    }
  }
}

void HalfSpaceFillCell(const std::vector<std::array<float, 2>> &FillPoints,
                       const std::vector<std::array<float, 2>> &gridPoint,
                       const std::array<float, 3> &color) {
  std::array<float, 4> BBox = GetBbox(FillPoints);
  std::array<int, 4> cellBBox = BboxToCellBbox(BBox);

  for (int x = cellBBox[0]; x < cellBBox[1]; x++) {
    for (int y = cellBBox[2]; y < cellBBox[3]; y++) {
      bool inside = true;

      for (int i = 0; i < gridPoint.size(); i++) {
        double determinant =
            gridPoint[i][0] * (gridPoint[(i + 1) % gridPoint.size()][1] - y) +
            gridPoint[(i + 1) % gridPoint.size()][0] * (y - gridPoint[i][1]) +
            x * (gridPoint[i][1] - gridPoint[(i + 1) % gridPoint.size()][1]);

        if (determinant < 0.0) {
          inside = false;
        }
      }

      if (inside) {
        if (_gridDraw[y + x * size_y].first) {
          continue;
        }
        _waitToDraw.push({y + x * size_y, {color[0], color[1], color[2]}});
      }
    }
  }
}
