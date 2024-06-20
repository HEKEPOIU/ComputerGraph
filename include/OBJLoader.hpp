#ifndef OBJLOADER_HPP
#define OBJLOADER_HPP

#include "DrawableObject.hpp"
#include <array>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdio.h>
#include <string>
#include <unordered_map>
#include <vector>

class OBJLoader {
public:
  OBJLoader() : _cache({}){};

  std::shared_ptr<DrawableObject> &get_OBJ(const std::string &path) {
    if (_cache.find(path) != _cache.end()) {
      return _cache[path];
    }
    std::ifstream obj_file(path);
    if (!obj_file) {
      throw std::runtime_error("Can't open file: " + path);
    }
    std::vector<std::array<float, 3>> vertex{};
    std::vector<std::vector<int>> faces{};
    std::string line{};

    while (getline(obj_file, line)) {
      if (line.empty()) {
        continue;
      }
      switch (line.front()) {
      case 'v': {
        float x, y, z;
        sscanf(line.c_str(), "v %f %f %f", &x, &y, &z);
        vertex.push_back({x, y, z});
        break;
      }
      case 'f': {
        std::istringstream iss(line);
        char f;
        iss >> f; // skip 'f'
        int index;
        std::vector<int> face;
        while (iss >> index) {
          face.push_back(index);
        }
        faces.push_back(face);
        break;
      }
      default:
        break;
      }
    }

    _cache.insert({ path,std::make_shared<DrawableObject>(vertex, faces) });
    return _cache[path];
  }

private:
  std::unordered_map<std::string, std::shared_ptr<DrawableObject>> _cache;
};

#endif