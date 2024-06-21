#ifndef IMAGELOADER_HPP
#define IMAGELOADER_HPP

#include <memory>
#include <stdexcept>
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

class ImageLoader {
private:
  struct StbiImageDeleter {
    void operator()(unsigned char *ptr) const { stbi_image_free(ptr); }
  };

public:
  std::unique_ptr<unsigned char, StbiImageDeleter> data;
  int width;
  int height;
  int channels;
  ImageLoader(const char *filename) {
    data = std::unique_ptr<unsigned char, StbiImageDeleter>(
        stbi_load(filename, &width, &height, &channels, 0), StbiImageDeleter());

    if (!data) {
      throw std::runtime_error(std::string("Failed to load image: ") + filename);
    }
  }
};

#endif
