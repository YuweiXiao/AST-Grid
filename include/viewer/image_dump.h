#pragma once
#include <string>
#include "general.h"

#if defined(_WIN64)
	#define STBI_MSC_SECURE_CRT
#endif
#include <stb/stb_image_write.h>
// #include "omni_grid3.h"
// #include "omni_grid2.h"

namespace Omni {
namespace Viewer {

// convertGray - convert T into Real [0,1]
template<class T>
void grid2PNGImage(std::string filename, T* data, Size2 size, std::function<Real(T)> convertGray) {
	char *cdata = new char[size.x() * size.y()];
	for (int j = size.y() - 1; j >= 0; j--) {
		for (int i = 0; i < size.x(); i++) {
			Real value = convertGray(data[j*size.x() + i]);
			cdata[(size.y()-1-j)*size.x() + i] = (int)(255 * value);
		}
	}
	stbi_write_png(filename.c_str(), size.x(), size.y(), 1, cdata, size.x());
	return;
}

void grid2PNGImage(std::string filename, Size2 size, std::function<Real(int xIdx, int yIdx)> func, bool parallel=false);

void grid2PNGColorImage(std::string filename, Size2 size, std::function<Vector3f(int xIdx, int yIdx)> func);

} // end of namespace Viewer
} // end of namespace Omni