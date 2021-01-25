#pragma once
#include "general.h"
#include "solver.h"
#include "octree_grid_layout3.h"
#include "octree_omni_grid3.h"
#include "octree_omni_face_centered_grid3.h"
#include "octree_omni_tilt_grid3.h"
#include "pcgsolver/pcg_solver.h"
#include "octree_omni_boundary_solver3.h"
#include "omni_boundary_solver3.h"
#include "narrowband_levelset.h"
#include "atomic_counter.h"
#include "octree_omni_grid_layout.h"
#ifdef CUDA_SOLVER
#include "cuda_solver.cuh"
#endif

#define ENABLE_SURFACE_TENSION 
#define ENABLE_GHOST_FLUID

namespace Omni {

class OctreeOmniGridSolver3 : public Solver {
public:
    OctreeOmniGridSolver3(const Size3& res, REAL spac, const std::string& name="OctreeOmniGridSolver3");
    void initialize() override;
    int numberOfSubTimeSteps(REAL dt) const override;
    virtual int cfl(REAL dt) const;
    void setMaxCfl(REAL cfl) { maxCfl = cfl; }
    void setRho(REAL r) { rho = r; }
    virtual NarrowbandLevelSet* getLevelSet() {return nullptr;}
    AtomicCounter& getACounter() {return aCounter;}

    // scene setup function
    std::function<void(OctreeOmniGridSolver3* solver)> sceneSetup;
    std::function<void(OctreeOmniGridSolver3* solver)> splitMapSetup;
    std::function<void(OctreeOmniGridSolver3* solver)> forceModifyCB;
    std::function<void(OctreeOmniGridSolver3* solver)> constrainVelCB;
    OctreeOmniScalarGrid3Ptr& getDensity() { return density; }
    OctreeOmniScalarGrid3Ptr& getHeat() { return heat; }
    OctreeOmniFaceCenteredGrid3Ptr& getVelocity() {return velocity;}
    // OctreeOmniFaceCenteredGrid3Ptr& getForceField() {return force;}
    OctreeGridLayout3Ptr& getLayout() {return layout;}
    OctreeTiltENodeGrid3Ptr& getTiltEGrid() {return tiltE;}
    CellCenteredFlagGrid3Ptr& getSplitMap() {return splitMap;}

    bool isDumpDensity = false;
	bool isRestart = false;
    bool saveStatus = false;

	bool restart(int frame, const string& layoutfile, const string& boundaryfile, const string& tiltfile,
		const string& densityfile, const string& heatfile, const string& velocityfile);

    bool enableGravity = true;
    int numOctreeLayer = 4;
    bool enableTiltAdaptivity = true;

protected:
    void onUpdate(REAL timeStep) override;
    void saveResult(int index) override;
    void initData() override;
    void initObstacle() override;
    virtual void initCuttingCellWeights();
    virtual void advect(REAL dt);
    virtual void addForce(REAL dt);
    virtual int updateIndex();
    virtual int updateIndexParallel();
    void adjustE(REAL dThreshold, REAL dt, REAL vThreshold = EPSILON, REAL maxThreshold = 0.8);
    virtual void updateGhostCells() {}

    void extrapolateVelocity(OctreeFaceCenteredGrid3Ptr& velocity);

    void project(REAL dt);
#ifdef CUDA_SOLVER
    void updateVelocityByPressure(const std::vector<REAL>& p, REAL dt);
#else
    void updateVelocityByPressure(const VectorXd& p, REAL dt);
#endif
	void wipeBoundary();
    
    // parameter
    REAL rho = 2;
    REAL beta = 3;
    REAL maxCfl = 3.0;
    // REAL gamma = 0.073;
    REAL gamma = 0.2;
    Size3 resolution;
    REAL delta_h;

    // simulation data
    OctreeTiltENodeGrid3Ptr tiltE;
    OctreeGridLayout3Ptr layout;
    CellCenteredFlagGrid3Ptr splitMap;
    OctreeOmniFaceCenteredGrid3Ptr velocity;
    OctreeOmniFaceCenteredGrid3Ptr weightsCC;
    OctreeOmniScalarGrid3Ptr heat;
    OctreeOmniScalarGrid3Ptr density;
    OctreeOmniScalarGrid3Ptr pressure;
    OctreeOmniFlagGrid3Ptr pIndex;
    OctreeOmniGridLayoutPtr<3> tilt_layout;
#ifdef CUDA_SOLVER
    std::vector<REAL> b, vRes;
    FixedSparseMatrix<REAL> fixedA;
    CudaSolver cudaSolver;
#else
	VectorXd b, vRes;
#endif

    OctreeTiltENodeGrid3Ptr tTiltE;
    OctreeOmniFaceCenteredGrid3Ptr tVelocity;
    OctreeOmniScalarGrid3Ptr tHeat;
    OctreeOmniScalarGrid3Ptr tDensity;

    // solvers
#ifndef CUDA_SOLVER
    PCGSolver<REAL> pcgSolver;
#endif
    PCG::SparseMatrix<REAL> A;
	OmniBoundarySolver3Ptr boundarySolver;
    OctreeCellCenteredCharGrid3Ptr marker0 = nullptr, marker1 = nullptr;
    CellCenteredScalarGrid3Ptr boundarySdf;
    AtomicCounter aCounter;

protected:
	REAL _buoyancySmokeDensityFactor = -0.000625;
	REAL _buoyancyTemperatureFactor = 5.0; // for large-scale demo, the value is 0.5
	REAL _smokeDecayFactor = 0.999;
	REAL _temperatureDecayFactor = 0.999;
};

}   // end of namespace Omni