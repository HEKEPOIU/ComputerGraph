#include "Texture.hpp"
#include "ImageLoader.hpp"
#include "freeglut.h"
#include "freeglut_std.h"
#include <memory>
Texture::Texture(const std::string &filePath) {
  if (filePath != "") {
    SetUpTexture(filePath);
  }
}
Texture::~Texture() { glDeleteTextures(1, &textureId); }

void Texture::SetUpTexture(const std::string &filePath) {

  glGenTextures(1, &textureId);
  glBindTexture(GL_TEXTURE_2D, textureId);
  GLenum eFormat;

  imageData = std::make_unique<ImageLoader>(filePath.c_str());
  if (imageData->data == nullptr)
    return;
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
                    imageData->data);

  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
                  GL_LINEAR_MIPMAP_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
}

void Texture::ApplyTexture() { glBindTexture(GL_TEXTURE_2D, textureId); }
