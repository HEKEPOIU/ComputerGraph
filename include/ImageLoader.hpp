#ifndef IMAGELOADER_HPP
#define IMAGELOADER_HPP

#include "stb_image.h"
#include <iostream>
#include <ostream>
#include <string>

// Simple Wrapper for stb_image, If Failed to load Image, data will be null.
class ImageLoader {
public:
  unsigned char *data;
  int width;
  int height;
  int channels;
  ImageLoader(const std::string &filename) {
    data = stbi_load(filename.c_str(), &width, &height, &channels, 0);

    if (!data) {
      std::cerr << "Failed to load image: " << filename << std::endl;
    }
  }
  ~ImageLoader() { stbi_image_free(data); }
};

#endif
