#include "viewer/grid_dump.h"
#include <fstream>
#include "util.h"
#define STB_IMAGE_IMPLEMENTATION
#include <stb/stb_image.h>
#include "jet/fmm_levelset_solver2.h"

namespace Omni {

void dumpSdf3D(const std::string &filename, CellCenteredScalarGrid3Ptr grid) {
	if (!isDirExist(filename))
		createDir(filename);
	std::ofstream out(filename, std::ios::binary);
	if (!out) {
		throw std::runtime_error("[dumpOctreeOmniGrid3D] can not open file");
	}
	Vector3f domain(grid->resolution().x() * grid->gridSpacing(),
		grid->resolution().y() * grid->gridSpacing(),
		grid->resolution().z() * grid->gridSpacing());
	Size3 res = grid->resolution();
	// REAL xstep = domain.x() / res.x();
	// REAL ystep = domain.y() / res.y();
	// REAL zstep = domain.z() / res.z();
	out.write((char *)&res.x(), sizeof(int));
	out.write((char *)&res.y(), sizeof(int));
	out.write((char *)&res.z(), sizeof(int));
	out.write((char *)&domain.x(), sizeof(REAL));
	for (int z = 0; z < res.z(); z++) {
		for (int y = 0; y < res.y(); y++) {
			for (int x = 0; x < res.x(); x++) {
				REAL d = grid->get(x, y, z);
				out.write((char *)&d, sizeof(REAL));
			}
		}
	}
}


bool readSdf3D(const std::string & filename, CellCenteredScalarGrid3Ptr grid)
{
	std::ifstream in(filename, std::ios::binary);
	if (!in) {
		return false;
	}
	int xres, yres, zres;
	REAL delta_h, v;
	in.read((char *)&xres, sizeof(int));
	in.read((char *)&yres, sizeof(int));
	in.read((char *)&zres, sizeof(int));
	in.read((char *)&delta_h, sizeof(REAL));
	
	if (xres != grid->resolution().x() || yres != grid->resolution().y() || zres != grid->resolution().z()) {
		return false;
	}
	for (int z = 0; z < zres; z++) {
		for (int y = 0; y < yres; y++) {
			for (int x = 0; x < xres; x++) {
				in.read((char *)&v, sizeof(REAL));
				grid->set(x, y, z, v);
			}
		}
	}
	return true;
}

void dumpOctreeOmniGrid3DBinary(const std::string & filename, Size3 res, OctreeOmniScalarGrid3 &grid)
{
	if (!isDirExist(filename))
		createDir(filename);
	std::ofstream out(filename, std::ios::binary);
	if (!out) {
		throw std::runtime_error("[dumpOctreeOmniGrid3D] can not open file");
	}
	Vector3f domain(grid.resolution().x() * grid.gridSpacing(),
		grid.resolution().y() * grid.gridSpacing(),
		grid.resolution().z() * grid.gridSpacing());
	REAL xstep = domain.x() / res.x();
	REAL ystep = domain.y() / res.y();
	REAL zstep = domain.z() / res.z();
	out.write((char *)&res.x(), sizeof(int));
	out.write((char *)&res.y(), sizeof(int));
	out.write((char *)&res.z(), sizeof(int));
	out << "\n";
	out.write((char *)&domain.x(), sizeof(REAL));
	out.write((char *)&domain.y(), sizeof(REAL));
	out.write((char *)&domain.z(), sizeof(REAL));
	out << "\n";
	for (int z = 0; z < res.z(); z++) {
		for (int y = 0; y < res.y(); y++) {
			for (int x = 0; x < res.x(); x++) {
				REAL d = grid.sample((x + 0.5) * xstep,	(y + 0.5) * ystep, (z + 0.5) * zstep);
				out.write((char *)&d, sizeof(REAL));
			}
		}
	}
}

void dumpOctreeGrid3DBinary(const std::string &filename, const Size3& res, OctreeCellCenteredScalarGrid3 &grid) {
	if (!isDirExist(filename))
		createDir(filename);
	std::ofstream out(filename, std::ios::binary);
	if (!out) {
		throw std::runtime_error("[dumpOctreeOmniGrid3D] can not open file");
	}
	Vector3f domain = grid.resolution().template cast<Real>() * grid.gridSpacing();
	Real xstep = domain.x() / res.x();
	Real ystep = domain.y() / res.y();
	Real zstep = domain.z() / res.z();
	out.write((char *)&res.x(), sizeof(int));
	out.write((char *)&res.y(), sizeof(int));
	out.write((char *)&res.z(), sizeof(int));
	out << "\n";
	out.write((char *)&domain.x(), sizeof(Real));
	out.write((char *)&domain.y(), sizeof(Real));
	out.write((char *)&domain.z(), sizeof(Real));
	out << "\n";

	std::vector<Real> data(res.x() * res.y() * res.z());
	int resXY = res.x() * res.y();
	ThreadPool::parallelForTF(0, res.z(), [&](int z) {
		for (int y = 0; y < res.y(); y++) {
			for (int x = 0; x < res.x(); x++) {
				data[z * resXY + y * res.x() + x] = grid.sample((x + 0.5) * xstep,	(y + 0.5) * ystep, (z + 0.5) * zstep);
			}
		}
	});
	out.write((char*)&data[0], res.x() * res.y() * res.z() * sizeof(Real));
	out.close();
}

void dumpLevelSetGrid3DBinary(const std::string &filename, CellCenteredScalarGrid3 &grid) {
	if (!isDirExist(filename))
		createDir(filename);
	std::ofstream out(filename, std::ios::binary);
	if (!out) {
		throw std::runtime_error("[dumpOctreeOmniGrid3D] can not open file");
	}
	Vector3f domain(grid.resolution().x() * grid.gridSpacing(),
		grid.resolution().y() * grid.gridSpacing(),
		grid.resolution().z() * grid.gridSpacing());
	Size3 res = grid.resolution();
	// REAL xstep = domain.x() / res.x();
	// REAL ystep = domain.y() / res.y();
	// REAL zstep = domain.z() / res.z();
	out.write((char *)&res.x(), sizeof(int));
	out.write((char *)&res.y(), sizeof(int));
	out.write((char *)&res.z(), sizeof(int));
	out.write((char *)&domain.x(), sizeof(REAL));
	out.write((char *)&domain.y(), sizeof(REAL));
	out.write((char *)&domain.z(), sizeof(REAL));
	for (int z = 0; z < res.z(); z++) {
		for (int y = 0; y < res.y(); y++) {
			for (int x = 0; x < res.x(); x++) {
				float d = static_cast<float>(grid.get(x, y, z));
				out.write((char *)&d, sizeof(float));
			}
		}
	}
}

bool readLevelSetGrid3DBinary(const std::string &filename, CellCenteredScalarGrid3 &grid) {
	std::ifstream in(filename, std::ios::binary);
	if (!in) {
		spdlog::error("readLevelSetGrid3DBinary::can not open file: {}", filename);
		return false;
	}
	int resx, resy, resz;
	REAL delta_h;
	float v;

	in.read((char *)&resx, sizeof(int));
	in.read((char *)&resy, sizeof(int));
	in.read((char *)&resz, sizeof(int));
	if(resx != grid.resolution().x() || resy != grid.resolution().y() || resz != grid.resolution().z()) {
		spdlog::error("readLevelSetGrid3DBinary::data format differ: {}", filename);
		return false;
	}
	in.read((char *)&delta_h, sizeof(REAL));
	in.read((char *)&delta_h, sizeof(REAL));
	in.read((char *)&delta_h, sizeof(REAL));
	delta_h = delta_h / resx;
	for (int z = 0; z < resz; z++) {
		for (int y = 0; y < resy; y++) {
			for (int x = 0; x < resx; x++) {
				in.read((char *)&v, sizeof(float));
				grid.set(x, y, z, static_cast<REAL>(v));
			}
		}
	}
	return true;
}

void dumpLevelSetGrid3DBinary(const std::string &filename, NarrowbandLevelSet &grid) {
	if (!isDirExist(filename))
		createDir(filename);
	std::ofstream out(filename, std::ios::binary);
	if (!out) {
		throw std::runtime_error("[dumpOctreeOmniGrid3D] can not open file");
	}
	Size3 res = grid.resolution();
	out.write((char *)&res.x(), sizeof(int));
	out.write((char *)&res.y(), sizeof(int));
	out.write((char *)&res.z(), sizeof(int));
	REAL delta_h = grid.gridSpacing();
	out.write((char *)&delta_h, sizeof(REAL));
	auto inOutMarker = grid.getInOutMarker();
	auto narrowbandMarker = grid.getNarrowbandMarker();
	int dSize = inOutMarker->getData().size();
	out.write((char *)&dSize, sizeof(int));
	for(int i:inOutMarker->getData()) {
		out.write((char *)&i, sizeof(int));
	}
	dSize = narrowbandMarker->getData().size();
	out.write((char *)&dSize, sizeof(int));
	for(int i:narrowbandMarker->getData()) {
		out.write((char *)&i, sizeof(int));
	}

	for (int z = 0; z < res.z(); z++) {
		for (int y = 0; y < res.y(); y++) {
			for (int x = 0; x < res.x(); x++) {
				if(grid.isExist(x, y, z)) {
					float d = static_cast<float>(grid.get(x, y, z));
					out.write((char *)&d, sizeof(float));
				}
			}
		}
	}
}

bool readLevelSetGrid3DBinary(const std::string &filename, NarrowbandLevelSet &grid) {
	std::ifstream in(filename, std::ios::binary);
	if (!in) {
		spdlog::error("readLevelSetGrid3DBinary::can not open file: {}", filename);
		return false;
	}
	int resx, resy, resz;
	REAL delta_h;
	float v;

	in.read((char *)&resx, sizeof(int));
	in.read((char *)&resy, sizeof(int));
	in.read((char *)&resz, sizeof(int));
	if(resx != grid.resolution().x() || resy != grid.resolution().y() || resz != grid.resolution().z()) {
		spdlog::error("readLevelSetGrid3DBinary::data format differ: {}", filename);
		return false;
	}
	in.read((char *)&delta_h, sizeof(REAL));

	auto inOutMarker = make_shared<BitMap>(Size3(resx, resy, resz), delta_h);
	auto narrowbandMarker = make_shared<BitMap>(Size3(resx, resy, resz), delta_h);
	int size; //, tItem;
	in.read((char *)&size, sizeof(int));
	for(int i=0; i < size; ++i) {
		in.read((char *)&inOutMarker->getData()[i], sizeof(int));
	}
	in.read((char *)&size, sizeof(int));
	for(int i=0; i < size; ++i) {
		in.read((char *)&narrowbandMarker->getData()[i], sizeof(int));
	}

	grid.setInOutMarker(inOutMarker);
	grid.setNarrowbandMarker(narrowbandMarker);
	auto levelSet = make_shared<BitMapCellCenteredScalarGrid3>(Size3(resx, resy, resz), delta_h, narrowbandMarker);
	grid.setLevelSet(levelSet);
	for (int z = 0; z < resz; z++) {
		for (int y = 0; y < resy; y++) {
			for (int x = 0; x < resx; x++) {
				if(grid.isExist(x, y, z)) {
					in.read((char *)&v, sizeof(float));
					grid.set(x, y, z, v);
				}
			}
		}
	}
	return true;
}


void dumpOctreeLayout3dBinary(const std::string & filename, OctreeGridLayout3 & grid)
{
	if (!isDirExist(filename))
		createDir(filename);
	std::ofstream out(filename, std::ios::binary);
	if (!out) {
		throw std::runtime_error("[dumpOctreeLayoutGrid] can not open file");
	}

	//dump layer
	int level = grid.level();
	int size = 0;
	out.write((char *)&level, sizeof(int));
	for (int i = 0; i < level; i++) {
		auto layer = grid.getLayerPtr(i);
		auto& map = layer->getMapTable();
		size = map.size();
		out.write((char *)&size, sizeof(int));
		for (auto iter = map.begin(); iter != map.end(); ++iter) {
			LCOOR_T coor = iter->first;
			LayerNode3 node = iter->second;
			out.write((char *)&coor, sizeof(LCOOR_T));
			out.write((char *)&node, sizeof(LayerNode3));
		}
		auto& ghost = layer->getGhostNodeIdx();
		size = ghost.size();
		out.write((char *)&size, sizeof(int));
		for (int j = 0; j < size; j++) {
			LCOOR_T coor = ghost[j];
			out.write((char *)&coor, sizeof(LCOOR_T));
		}
		auto& node = layer->getNodeIdx();
		size = node.size();
		out.write((char *)&size, sizeof(int));
		for (int j = 0; j < size; j++) {
			LCOOR_T coor = node[j];
			out.write((char *)&coor, sizeof(LCOOR_T));
		}
	}
	//dump halfGhostIdx
	auto& halfGhost = grid.getHalfGhostIdx();
	size = halfGhost.size(); 
	out.write((char *)&size, sizeof(int));
	for (int j = 0; j < size; j++) {
		LCOOR_T coor = halfGhost[j];
		out.write((char *)&coor, sizeof(LCOOR_T));
	}
}

bool readOctreeLayout3dBinary(const std::string & filename, OctreeGridLayout3 & grid)
{
	std::ifstream in(filename, std::ios::binary);
	if (!in) {
		return false;
	}

	int level = 0;
	int size = 0;
	in.read((char *)&level, sizeof(int));
	for (int i = 0; i < level; i++) {
		auto layer = grid.getLayerPtr(i);
		auto& map = layer->getMapTable();
		in.read((char *)&size, sizeof(int));
		for (int j = 0; j < size; j++) {
			LCOOR_T coor;
			LayerNode3 node;
			in.read((char *)&coor, sizeof(LCOOR_T));
			in.read((char *)&node, sizeof(LayerNode3));
			map.insert(std::pair<LCOOR_T, LayerNode3>(coor, node));
		}
		auto& ghost = layer->getGhostNodeIdx();
		in.read((char *)&size, sizeof(int));
		ghost.resize(size);
		for (int j = 0; j < size; j++) {
			LCOOR_T coor;
			in.read((char *)&coor, sizeof(LCOOR_T));
			ghost[j] = coor;
		}
		auto& node = layer->getNodeIdx();
		in.read((char *)&size, sizeof(int));
		node.resize(size);
		for (int j = 0; j < size; j++) {
			LCOOR_T coor;
			in.read((char *)&coor, sizeof(LCOOR_T));
			node[j] = coor;
		}
	}
	auto& halfGhost = grid.getHalfGhostIdx();
	in.read((char *)&size, sizeof(int));
	halfGhost.resize(size);
	for (int j = 0; j < size; j++) {
		LCOOR_T coor;
		in.read((char *)&coor, sizeof(LCOOR_T));
		halfGhost[j] = coor;
	}

	return true;
}

void dumpOctreeOmniScalar3dBinary(const std::string & filename, OctreeOmniScalarGrid3 & grid)
{
	if (!isDirExist(filename))
		createDir(filename);
	std::ofstream out(filename, std::ios::binary);
	if (!out) {
		throw std::runtime_error("[dumpOctreeLayoutGrid3D] can not open file");
	}
	int level;
	auto& octagonData = grid.getOctagonGrid()->getData();
	level = octagonData.size();
	out.write((char *)&level, sizeof(int));
	for (int i = 0; i < level; i++) {
		int size = octagonData[i].size();
		out.write((char *)&size, sizeof(int));
		for (int j = 0; j < size; j++) {
			auto value = octagonData[i][j];
			out.write((char *)&value, sizeof(REAL));
		}
	}
	auto& tiltData = grid.getTiltGrid()->getData();
	level = tiltData.size();
	out.write((char *)&level, sizeof(int));
	for (int i = 0; i < level; i++) {
		int size = tiltData[i].size();
		out.write((char *)&size, sizeof(int));
		for (int j = 0; j < size; j++) {
			auto value = tiltData[i][j];
			out.write((char *)&value, sizeof(REAL));
		}
	}
}

bool readOctreeOmniScalar3dBinary(const std::string & filename, OctreeOmniScalarGrid3 & grid)
{
	std::ifstream in(filename, std::ios::binary);
	if (!in) {
		return false;
	}
	int level = 0;
	auto& octagonData = grid.getOctagonGrid()->getData();
	in.read((char *)&level, sizeof(int));
	octagonData.resize(level);
	for (int i = 0; i < level; i++) {
		int size = 0;
		in.read((char *)&size, sizeof(int));
		octagonData[i].resize(size);
		for (int j = 0; j < size; j++) {
			REAL value = 0;
			in.read((char *)&value, sizeof(REAL));
			octagonData[i][j] = value;
		}
	}
	auto& tiltData = grid.getTiltGrid()->getData();
	in.read((char *)&level, sizeof(int));
	tiltData.resize(level);
	for (int i = 0; i < level; i++) {
		int size = 0;
		in.read((char *)&size, sizeof(int));
		tiltData[i].resize(size);
		for (int j = 0; j < size; j++) {
			REAL value = 0;
			in.read((char *)&value, sizeof(REAL));
			tiltData[i][j] = value;
		}
	}
	return true;
}

void dumpOctreeOmniTiltE3dBinary(const std::string & filename, OctreeTiltENodeGrid3 & grid)
{
	if (!isDirExist(filename))
		createDir(filename);
	std::ofstream out(filename, std::ios::binary);
	if (!out) {
		throw std::runtime_error("[dumpOctreeTiltEGrid3D] can not open file");
	}
	auto& data = grid.getData();
	int level = data.size();
	out.write((char *)&level, sizeof(int));
	for (int i = 0; i < level; i++) {
		int size = data[i].size();
		out.write((char *)&size, sizeof(int));
		for (int j = 0; j < size; j++) {
			auto tilt = data[i][j];
			out.write((char *)&tilt, sizeof(struct TiltENode));
		}
	}
}

bool readOctreeOmniTiltE3dBinary(const std::string & filename, OctreeTiltENodeGrid3 & grid)
{
	std::ifstream in(filename, std::ios::binary);
	if (!in) {
		return false;
	}
	auto& data = grid.getData();
	int level = 0;
	in.read((char *)&level, sizeof(int));
	data.resize(level);
	for (int i = 0; i < level; i++) {
		int size = 0;
		in.read((char *)&size, sizeof(int));
		data[i].resize(size);
		for (int j = 0; j < size; j++) {
			TiltENode tilt;
			in.read((char *)&tilt, sizeof(struct TiltENode));
			data[i][j] = tilt;
		}
	}
	return true;
}

void dumpOctreeOmniFaceCentered3dBinary(const std::string & filename, OctreeOmniFaceCenteredGrid3 & grid)
{
	if (!isDirExist(filename))
		createDir(filename);
	std::ofstream out(filename, std::ios::binary);
	if (!out) {
		throw std::runtime_error("[dumpOctreeTiltEGrid3D] can not open file");
	}
	int level = 0;
	auto& ux = grid.getMACGrid()->getUx();
	level = ux.size();
	out.write((char *)&level, sizeof(int));
	for (int i = 0; i < level; i++) {
		int size = ux[i].size();
		out.write((char *)&size, sizeof(int));
		for (int j = 0; j < size; j++) {
			auto value = ux[i][j];
			out.write((char *)&value, sizeof(REAL));
		}
	}
	auto& uy = grid.getMACGrid()->getUy();
	level = uy.size();
	out.write((char *)&level, sizeof(int));
	for (int i = 0; i < level; i++) {
		int size = uy[i].size();
		out.write((char *)&size, sizeof(int));
		for (int j = 0; j < size; j++) {
			auto value = uy[i][j];
			out.write((char *)&value, sizeof(REAL));
		}
	}

	auto& uz = grid.getMACGrid()->getUz();
	level = uz.size();
	out.write((char *)&level, sizeof(int));
	for (int i = 0; i < level; i++) {
		int size = uz[i].size();
		out.write((char *)&size, sizeof(int));
		for (int j = 0; j < size; j++) {
			auto value = uz[i][j];
			out.write((char *)&value, sizeof(REAL));
		}
	}

	auto& tilt = grid.getTiltGrid()->getGridPtr()->getData();
	level = tilt.size();
	out.write((char *)&level, sizeof(int));
	for (int i = 0; i < level; i++) {
		int size = tilt[i].size();
		out.write((char *)&size, sizeof(int));
		for (int j = 0; j < size; j++) {
			Vector8f value = tilt[i][j];
			for (int k = 0; k < 8; k++) {
				out.write((char *)&(value[k]), sizeof(REAL));
			}
		}
	}
}

bool readOctreeOmniFaceCentered3dBinary(const std::string & filename, OctreeOmniFaceCenteredGrid3 & grid)
{
	std::ifstream in(filename, std::ios::binary);
	if (!in) {
		return false;
	}
	int level = 0;
	auto& ux = grid.getMACGrid()->getUx();
	in.read((char *)&level, sizeof(int));
	ux.resize(level);
	for (int i = 0; i < level; i++) {
		int size = 0;
		in.read((char *)&size, sizeof(int));
		ux[i].resize(size);
		for (int j = 0; j < size; j++) {
			REAL value = 0;
			in.read((char *)&value, sizeof(REAL));
			ux[i][j] = value;
		}
	}
	auto& uy = grid.getMACGrid()->getUy();
	in.read((char *)&level, sizeof(int));
	uy.resize(level);
	for (int i = 0; i < level; i++) {
		int size = 0;
		in.read((char *)&size, sizeof(int));
		uy[i].resize(size);
		for (int j = 0; j < size; j++) {
			REAL value = 0;
			in.read((char *)&value, sizeof(REAL));
			uy[i][j] = value;
		}
	}
	auto& uz = grid.getMACGrid()->getUz();
	in.read((char *)&level, sizeof(int));
	uz.resize(level);
	for (int i = 0; i < level; i++) {
		int size = 0;
		in.read((char *)&size, sizeof(int));
		uz[i].resize(size);
		for (int j = 0; j < size; j++) {
			REAL value = 0;
			in.read((char *)&value, sizeof(REAL));
			uz[i][j] = value;
		}
	}

	auto& tilt = grid.getTiltGrid()->getGridPtr()->getData();
	in.read((char *)&level, sizeof(int));
	tilt.resize(level);
	for (int i = 0; i < level; i++) {
		int size = 0;
		in.read((char *)&size, sizeof(int));
		tilt[i].resize(size);
		for (int j = 0; j < size; j++) {
			Vector8f value = Vector8f::Zero();
			for (int k = 0; k < 8; k++) {
				in.read((char *)&(value[k]), sizeof(REAL));
			}
			tilt[i][j] = value;
		}
	}
	return true;
}

}