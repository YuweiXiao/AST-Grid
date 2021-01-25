#include <fstream>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "viewer/image_dump.h"
#include "thread_pool.h"
#include "util.h"

namespace Omni {
namespace Viewer {

void grid2PNGImage(std::string filename, Size2 size, std::function<Real(int xIdx, int yIdx)> func, bool parallel) {
	char *cdata = new char[size.x() * size.y()];
	if(parallel) {
		ThreadPool::parallelForTF(0, size.y(), [&](int j){
			for (int i = 0; i < size.x(); i++) {
				REAL value = func(i, j);
				cdata[(size.y() - 1 - j)*size.x() + i] = (int)(255 * value);
			}
		});
	} else {
		for (int j = size.y() - 1; j >= 0; j--) {
			for (int i = 0; i < size.x(); i++) {
				REAL value = func(i, j);
				cdata[(size.y() - 1 - j)*size.x() + i] = (int)(255 * value);
			}
		}
	}
	if (!isDirExist(filename))
		createDir(filename);
	stbi_write_png(filename.c_str(), size.x(), size.y(), 1, cdata, size.x());
	delete[] cdata;
	return;
}

void grid2PNGColorImage(std::string filename, Size2 size, std::function<Vector3f(int xIdx, int yIdx)> func)
{
	unsigned char *cdata = new unsigned char[size.x() * size.y() * 3];
	ThreadPool::parallelForTF(0, int(size.y()), [&](int j){
	// for (int j = size.y() - 1; j >= 0; j--) {
		for (int i = 0; i < size.x(); i++) {
			Vector3f colorVec = func(i, j);
			cdata[3 *(size.y() - 1 - j) * size.x() + 3 * i] = (unsigned char)(255 * colorVec[0]);
			cdata[3 * (size.y() - 1 - j)* size.x() + 3 * i + 1] = (unsigned char)(255 * colorVec[1]);
			cdata[3 * (size.y() - 1 - j)* size.x() + 3 * i + 2] = (unsigned char)(255 * colorVec[2]);
		}
	});
	if (!isDirExist(filename))
		createDir(filename);
	stbi_write_png(filename.c_str(), size.x(), size.y(), 3, cdata, size.x() * 3);
	delete[] cdata;
	return;
}

} // end of namespace Viewer
} // end of namespace Omni