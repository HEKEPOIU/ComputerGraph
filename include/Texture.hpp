#ifndef TEXTURE_HPP
#define TEXTURE_HPP

#include "ImageLoader.hpp"
#include <memory>
#include <string>
#define GL_CLAMP_TO_EDGE 0x812F
class Texture {
public:
  Texture(const std::string &filePath = "");
  void SetUpTexture(const std::string &filePath);
  void ApplyTexture();
  ~Texture();

private:
  unsigned int textureId;
  std::unique_ptr<ImageLoader> imageData;
};

#endif
