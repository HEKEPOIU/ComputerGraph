#ifndef OBJLOADER_HPP
#define OBJLOADER_HPP

#include "CGMath.hpp"
#include "DrawableObject.hpp"
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
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
    std::vector<Vec3> vertex{};
    std::vector<Vec2> tex_coords;
    std::vector<Vec3> normals;
    std::vector<std::vector<int>> faces{};
    std::vector<std::vector<int>> tex_coord_indices;
    std::vector<std::vector<int>> normal_indices;
    std::string line{};

    while (getline(obj_file, line)) {
      if (line.empty() || line.front() == '#') {
        continue; // 忽略空行和注释行
      }
      std::istringstream iss(line);
      std::string prefix;
      iss >> prefix;

      if (prefix == "v") {
        double x, y, z;
        iss >> x >> y >> z;
        vertex.push_back({x, y, z});
      } else if (prefix == "vt") {
        double u, v;
        iss >> u >> v;
        tex_coords.push_back({u, v});
      } else if (prefix == "vn") {
        double nx, ny, nz;
        iss >> nx >> ny >> nz;
        normals.push_back({nx, ny, nz});
      } else if (prefix == "f") {
        std::vector<int> vertex_index;
        std::vector<int> tex_coord_index;
        std::vector<int> normal_index;

        std::string vertex_data;
        while (iss >> vertex_data) {
          std::istringstream vertex_iss(vertex_data);
          std::string vertex_part;
          int idx[3] = {0, 0, 0}; // 初始化为0，表示无索引

          for (int i = 0; std::getline(vertex_iss, vertex_part, '/'); ++i) {
            if (!vertex_part.empty()) {
              idx[i] = std::stoi(vertex_part);
            }
          }
          vertex_index.push_back(idx[0] - 1); // 索引从0开始
          if (idx[1] != 0) {
            tex_coord_index.push_back(idx[1] - 1);
          }
          if (idx[2] != 0) {
            normal_index.push_back(idx[2] - 1);
          }
        }
        faces.push_back(vertex_index);
        if (!tex_coord_index.empty()) {
          tex_coord_indices.push_back(tex_coord_index);
        }
        if (!normal_index.empty()) {
          normal_indices.push_back(normal_index);
        }
      }
    }
    _cache.insert({path, std::make_shared<DrawableObject>(
                             vertex, tex_coords, normals, faces,
                             tex_coord_indices, normal_indices)});
    return _cache[path];
  }

private:
  std::unordered_map<std::string, std::shared_ptr<DrawableObject>> _cache;
};

#endif
