#include "constants.h"
// #include "util.h"

using namespace Omni;

REAL Constants::OMNI_DUAL_GRAPH_MATRIX2[2][2][2][2][2][4] = {0}; // x, y, xx, yy
REAL Constants::OMNI_DUAL_GRAPH_MATRIX3[2][2][2][2][3][7] = {0};
// REAL Constants::OCTREE_OMNI_DUAL_GRAPH_MATRIX3[2][2][2][2][2][2][2][3][7] = {0};
Eigen::Matrix2d Constants::ROTATE45 = Eigen::Matrix2d::Identity();
Eigen::Matrix2d Constants::TILT_PRINCIPLE_AXIS2 = Eigen::Matrix2d::Identity();
Matrix37f Constants::OCTREE_OMNI_DUAL_GRAPH_EIGEN_MATRIX3[2][2][2][2][2][2][2] = {Matrix37f::Zero()};
Matrix37f Constants::OMNI_DUAL_GRAPH_EIGEN_MATRIX3[2][2][2][2] = {Matrix37f::Zero()};
Matrix34d Constants::OMNI_DUAL_GRAPH_TILT_EIGEN_MATRIX3 = Matrix34d::Zero();
Eigen::Matrix3d Constants::TILT_PRINCIPLE_AXIS[8] = {Eigen::Matrix3d::Identity()};

void Constants::init() {
	memset(Constants::OMNI_DUAL_GRAPH_MATRIX2, 0, sizeof(REAL) * 2 * 2 * 2 * 2 * 2 * 4);
    memset(Constants::OMNI_DUAL_GRAPH_MATRIX3, 0, sizeof(REAL) * 2 * 2 * 2 * 2 * 3 * 7);
    // memset(Constants::OCTREE_OMNI_DUAL_GRAPH_MATRIX3, 0, sizeof(REAL) * 2 * 2 * 2 * 2 * 2 * 2 * 2 * 3 * 7);

    TILT_PRINCIPLE_AXIS2.col(0) = Vector2f(SQRT2/2.0, -SQRT2/2.0);
    TILT_PRINCIPLE_AXIS2.col(1) = Vector2f(SQRT2/2.0, SQRT2/2.0);
    Constants::ROTATE45 << cos(PI/4), -sin(PI/4), sin(PI/4), cos(PI/4);

    OMNI_DUAL_GRAPH_TILT_EIGEN_MATRIX3 << INV_SQRT3, -INV_SQRT3, -INV_SQRT3, INV_SQRT3,
                                          INV_SQRT3, INV_SQRT3, -INV_SQRT3, -INV_SQRT3,
                                          INV_SQRT3, INV_SQRT3, INV_SQRT3, INV_SQRT3;
    OMNI_DUAL_GRAPH_TILT_EIGEN_MATRIX3 *= 3.0/4.0;

    // OMNI_DUAL_GRAPH_MATRIX3 initialization
    for(int n0 = 0; n0 <= 1; ++n0) {
    for(int n1 = 0; n1 <= 1; ++n1) {
    for(int n2 = 0; n2 <= 1; ++n2) {
    for(int n3 = 0; n3 <= 1; ++n3) {
        Eigen::MatrixXd m(n0+n1+n2+n3+3, 3);
        m(0, 0) = 1;m(0, 1) = 0;m(0, 2) = 0;
        m(1, 0) = 0;m(1, 1) = 1;m(1, 2) = 0;
        m(2, 0) = 0;m(2, 1) = 0;m(2, 2) = 1;
        int count = 3;
        if(n0) {
            m(count, 0) = INV_SQRT3;m(count, 1) = INV_SQRT3;m(count, 2) = INV_SQRT3;
            count++;
        }
        if(n1) {
            m(count, 0) = -INV_SQRT3;m(count, 1) = INV_SQRT3;m(count, 2) = INV_SQRT3;
            count++;
        }
        if(n2) {
            m(count, 0) = -INV_SQRT3;m(count, 1) = -INV_SQRT3;m(count, 2) = INV_SQRT3;
            count++;
        }
        if(n3) {
            m(count, 0) = INV_SQRT3;m(count, 1) = -INV_SQRT3;m(count, 2) = INV_SQRT3;
            count++;
        }

        Eigen::Matrix3Xd ans = (m.transpose() * m).inverse() * m.transpose();
        for(int k = 0; k < 3; ++k) {
            Constants::OMNI_DUAL_GRAPH_MATRIX3[n0][n1][n2][n3][k][0] = ans(k,0);   // x
            Constants::OMNI_DUAL_GRAPH_MATRIX3[n0][n1][n2][n3][k][0] = ans(k,0);
            Constants::OMNI_DUAL_GRAPH_MATRIX3[n0][n1][n2][n3][k][0] = ans(k,0);
            Constants::OMNI_DUAL_GRAPH_MATRIX3[n0][n1][n2][n3][k][1] = ans(k,1);   // y
            Constants::OMNI_DUAL_GRAPH_MATRIX3[n0][n1][n2][n3][k][1] = ans(k,1);
            Constants::OMNI_DUAL_GRAPH_MATRIX3[n0][n1][n2][n3][k][1] = ans(k,1);
            Constants::OMNI_DUAL_GRAPH_MATRIX3[n0][n1][n2][n3][k][2] = ans(k,2);   // z
            Constants::OMNI_DUAL_GRAPH_MATRIX3[n0][n1][n2][n3][k][2] = ans(k,2);
            Constants::OMNI_DUAL_GRAPH_MATRIX3[n0][n1][n2][n3][k][2] = ans(k,2);
            
            Constants::OMNI_DUAL_GRAPH_EIGEN_MATRIX3[n0][n1][n2][n3](k, 0) = ans(k,0);   // x
            Constants::OMNI_DUAL_GRAPH_EIGEN_MATRIX3[n0][n1][n2][n3](k, 0) = ans(k,0);
            Constants::OMNI_DUAL_GRAPH_EIGEN_MATRIX3[n0][n1][n2][n3](k, 0) = ans(k,0);
            Constants::OMNI_DUAL_GRAPH_EIGEN_MATRIX3[n0][n1][n2][n3](k, 1) = ans(k,1);   // y
            Constants::OMNI_DUAL_GRAPH_EIGEN_MATRIX3[n0][n1][n2][n3](k, 1) = ans(k,1);
            Constants::OMNI_DUAL_GRAPH_EIGEN_MATRIX3[n0][n1][n2][n3](k, 1) = ans(k,1);
            Constants::OMNI_DUAL_GRAPH_EIGEN_MATRIX3[n0][n1][n2][n3](k, 2) = ans(k,2);   // z
            Constants::OMNI_DUAL_GRAPH_EIGEN_MATRIX3[n0][n1][n2][n3](k, 2) = ans(k,2);
            Constants::OMNI_DUAL_GRAPH_EIGEN_MATRIX3[n0][n1][n2][n3](k, 2) = ans(k,2);
            count = 3;

            if(n0) {
                Constants::OMNI_DUAL_GRAPH_MATRIX3[n0][n1][n2][n3][k][3] = ans(k,count);
                Constants::OMNI_DUAL_GRAPH_MATRIX3[n0][n1][n2][n3][k][3] = ans(k,count);
                Constants::OMNI_DUAL_GRAPH_MATRIX3[n0][n1][n2][n3][k][3] = ans(k,count);
                
                Constants::OMNI_DUAL_GRAPH_EIGEN_MATRIX3[n0][n1][n2][n3](k, 3) = ans(k,count);
                Constants::OMNI_DUAL_GRAPH_EIGEN_MATRIX3[n0][n1][n2][n3](k, 3) = ans(k,count);
                Constants::OMNI_DUAL_GRAPH_EIGEN_MATRIX3[n0][n1][n2][n3](k, 3) = ans(k,count);
                count++;
            }
            if(n1) {
                Constants::OMNI_DUAL_GRAPH_MATRIX3[n0][n1][n2][n3][k][4] = ans(k,count);
                Constants::OMNI_DUAL_GRAPH_MATRIX3[n0][n1][n2][n3][k][4] = ans(k,count);
                Constants::OMNI_DUAL_GRAPH_MATRIX3[n0][n1][n2][n3][k][4] = ans(k,count);
                
                Constants::OMNI_DUAL_GRAPH_EIGEN_MATRIX3[n0][n1][n2][n3](k, 4) = ans(k,count);
                Constants::OMNI_DUAL_GRAPH_EIGEN_MATRIX3[n0][n1][n2][n3](k, 4) = ans(k,count);
                Constants::OMNI_DUAL_GRAPH_EIGEN_MATRIX3[n0][n1][n2][n3](k, 4) = ans(k,count);
                count++;
            }
            if(n2) {
                Constants::OMNI_DUAL_GRAPH_MATRIX3[n0][n1][n2][n3][k][5] = ans(k,count);
                Constants::OMNI_DUAL_GRAPH_MATRIX3[n0][n1][n2][n3][k][5] = ans(k,count);
                Constants::OMNI_DUAL_GRAPH_MATRIX3[n0][n1][n2][n3][k][5] = ans(k,count);
                
                Constants::OMNI_DUAL_GRAPH_EIGEN_MATRIX3[n0][n1][n2][n3](k, 5) = ans(k,count);
                Constants::OMNI_DUAL_GRAPH_EIGEN_MATRIX3[n0][n1][n2][n3](k, 5) = ans(k,count);
                Constants::OMNI_DUAL_GRAPH_EIGEN_MATRIX3[n0][n1][n2][n3](k, 5) = ans(k,count);
                count++;
            }
            if(n3) {
                Constants::OMNI_DUAL_GRAPH_MATRIX3[n0][n1][n2][n3][k][6] = ans(k,count);
                Constants::OMNI_DUAL_GRAPH_MATRIX3[n0][n1][n2][n3][k][6] = ans(k,count);
                Constants::OMNI_DUAL_GRAPH_MATRIX3[n0][n1][n2][n3][k][6] = ans(k,count);
                
                Constants::OMNI_DUAL_GRAPH_EIGEN_MATRIX3[n0][n1][n2][n3](k, 6) = ans(k,count);
                Constants::OMNI_DUAL_GRAPH_EIGEN_MATRIX3[n0][n1][n2][n3](k, 6) = ans(k,count);
                Constants::OMNI_DUAL_GRAPH_EIGEN_MATRIX3[n0][n1][n2][n3](k, 6) = ans(k,count);
            }
        }
    }}}}

    // OMNI_DUAL_GRAPH_MATRIX2 initialization
    for(int x = 0; x <= 1; ++x) {
    for(int y = 0; y <= 1; ++y) {
    for(int n0 = 0; n0 <= 1; ++n0) {
    for(int n1 = 0; n1 <= 1; ++n1) {
        if(x+y+n0+n1 <= 1)
            continue;
        Eigen::MatrixXd m(n0+n1+x+y, 2);
        int count = 0;
        if(x) {
            m(count, 0) = 1; m(count, 1) = 0;
            count++;
        }
        if(y) {
            m(count, 0) = 0; m(count, 1) = 1;
            count++;
        }
        if(n0) {
            m(count, 0) = INV_SQRT2; m(count, 1) = INV_SQRT2;
            count++;
        }
        if(n1) {
            m(count, 0) = INV_SQRT2; m(count, 1) = -INV_SQRT2;
            count++;
        }

        Eigen::Matrix2Xd ans = (m.transpose() * m).inverse() * m.transpose();
        for(int k = 0; k < 2; ++k) {
            count = 0;
            if(x) {
                Constants::OMNI_DUAL_GRAPH_MATRIX2[x][y][n0][n1][k][0] = ans(k, count);
                count++;
            }
            if(y) {
                Constants::OMNI_DUAL_GRAPH_MATRIX2[x][y][n0][n1][k][1] = ans(k, count);
                count++;
            }
            if(n0) {
                Constants::OMNI_DUAL_GRAPH_MATRIX2[x][y][n0][n1][k][2] = ans(k, count);
                count++;
            }
            if(n1) {
                Constants::OMNI_DUAL_GRAPH_MATRIX2[x][y][n0][n1][k][3] = ans(k, count);
                count++;
            }
        }
    }}}}

    // OCTREE_OMNI_DUAL_GRAPH_MATRIX3 initialization
    for(int x = 0; x <= 1; ++x) {
    for(int y = 0; y <= 1; ++y) {
    for(int z = 0; z <= 1; ++z) {
    for(int n0 = 0; n0 <= 1; ++n0) {
    for(int n1 = 0; n1 <= 1; ++n1) {
    for(int n2 = 0; n2 <= 1; ++n2) {
    for(int n3 = 0; n3 <= 1; ++n3) {
        if(x+y+z+n0+n1+n2+n3 <= 1)  continue;
        Eigen::MatrixXd m(x+y+z+n0+n1+n2+n3, 3);
        int count = 0;
        if(x) {
            m(count, 0) = 1;m(count, 1) = 0;m(count, 2) = 0;
            count++;
        }
        if(y) {
            m(count, 0) = 0;m(count, 1) = 1;m(count, 2) = 0;
            count++;
        }
        if(z) {
            m(count, 0) = 0;m(count, 1) = 0;m(count, 2) = 1;
            count++;
        }
        if(n0) {
            m(count, 0) = INV_SQRT3;m(count, 1) = INV_SQRT3;m(count, 2) = INV_SQRT3;
            count++;
        }
        if(n1) {
            m(count, 0) = -INV_SQRT3;m(count, 1) = INV_SQRT3;m(count, 2) = INV_SQRT3;
            count++;
        }
        if(n2) {
            m(count, 0) = -INV_SQRT3;m(count, 1) = -INV_SQRT3;m(count, 2) = INV_SQRT3;
            count++;
        }
        if(n3) {
            m(count, 0) = INV_SQRT3;m(count, 1) = -INV_SQRT3;m(count, 2) = INV_SQRT3;
            count++;
        }

        Eigen::Matrix3Xd ans = (m.transpose() * m).inverse() * m.transpose();
        Constants::OCTREE_OMNI_DUAL_GRAPH_EIGEN_MATRIX3[x][y][z][n0][n1][n2][n3] = Matrix37f::Zero();
        for(int k = 0; k < 3; ++k) {
            count = 0;
            if(x) {
                Constants::OCTREE_OMNI_DUAL_GRAPH_EIGEN_MATRIX3[x][y][z][n0][n1][n2][n3](k, 0) = ans(k,count);
                Constants::OCTREE_OMNI_DUAL_GRAPH_EIGEN_MATRIX3[x][y][z][n0][n1][n2][n3](k, 0) = ans(k,count);
                Constants::OCTREE_OMNI_DUAL_GRAPH_EIGEN_MATRIX3[x][y][z][n0][n1][n2][n3](k, 0) = ans(k,count);
                count++;
            }
            if(y) {
                Constants::OCTREE_OMNI_DUAL_GRAPH_EIGEN_MATRIX3[x][y][z][n0][n1][n2][n3](k, 1) = ans(k,count);
                Constants::OCTREE_OMNI_DUAL_GRAPH_EIGEN_MATRIX3[x][y][z][n0][n1][n2][n3](k, 1) = ans(k,count);
                Constants::OCTREE_OMNI_DUAL_GRAPH_EIGEN_MATRIX3[x][y][z][n0][n1][n2][n3](k, 1) = ans(k,count);
                count++;
            }
            if(z) {
                Constants::OCTREE_OMNI_DUAL_GRAPH_EIGEN_MATRIX3[x][y][z][n0][n1][n2][n3](k, 2) = ans(k,count);
                Constants::OCTREE_OMNI_DUAL_GRAPH_EIGEN_MATRIX3[x][y][z][n0][n1][n2][n3](k, 2) = ans(k,count);
                Constants::OCTREE_OMNI_DUAL_GRAPH_EIGEN_MATRIX3[x][y][z][n0][n1][n2][n3](k, 2) = ans(k,count);
                count++;
            }
            if(n0) {
                Constants::OCTREE_OMNI_DUAL_GRAPH_EIGEN_MATRIX3[x][y][z][n0][n1][n2][n3](k, 3) = ans(k,count);
                Constants::OCTREE_OMNI_DUAL_GRAPH_EIGEN_MATRIX3[x][y][z][n0][n1][n2][n3](k, 3) = ans(k,count);
                Constants::OCTREE_OMNI_DUAL_GRAPH_EIGEN_MATRIX3[x][y][z][n0][n1][n2][n3](k, 3) = ans(k,count);
                count++;
            }
            if(n1) {
                Constants::OCTREE_OMNI_DUAL_GRAPH_EIGEN_MATRIX3[x][y][z][n0][n1][n2][n3](k, 4) = ans(k,count);
                Constants::OCTREE_OMNI_DUAL_GRAPH_EIGEN_MATRIX3[x][y][z][n0][n1][n2][n3](k, 4) = ans(k,count);
                Constants::OCTREE_OMNI_DUAL_GRAPH_EIGEN_MATRIX3[x][y][z][n0][n1][n2][n3](k, 4) = ans(k,count);
                count++;
            }
            if(n2) {
                Constants::OCTREE_OMNI_DUAL_GRAPH_EIGEN_MATRIX3[x][y][z][n0][n1][n2][n3](k, 5) = ans(k,count);
                Constants::OCTREE_OMNI_DUAL_GRAPH_EIGEN_MATRIX3[x][y][z][n0][n1][n2][n3](k, 5) = ans(k,count);
                Constants::OCTREE_OMNI_DUAL_GRAPH_EIGEN_MATRIX3[x][y][z][n0][n1][n2][n3](k, 5) = ans(k,count);
                count++;
            }
            if(n3) {
                Constants::OCTREE_OMNI_DUAL_GRAPH_EIGEN_MATRIX3[x][y][z][n0][n1][n2][n3](k, 6) = ans(k,count);
                Constants::OCTREE_OMNI_DUAL_GRAPH_EIGEN_MATRIX3[x][y][z][n0][n1][n2][n3](k, 6) = ans(k,count);
                Constants::OCTREE_OMNI_DUAL_GRAPH_EIGEN_MATRIX3[x][y][z][n0][n1][n2][n3](k, 6) = ans(k,count);
            }
        }
    }}}}}}}


    for(int i = 0; i < 8; ++i) {
        REAL x_sign = (i==0||i==3||i==4||i==7?1:-1);
        REAL y_sign = (i==0||i==1||i==4||i==5?1:-1);
        REAL z_sign = i < 4 ? 1 : -1;
        Vector3f axis = TILT3_UNIT_DIR[i];
        Vector3f axis_1 = Vector3f((i%2) * x_sign, !(i%2) * y_sign, 0) - Vector3f(!(i%2) * x_sign, (i%2) * z_sign, 0);
        Vector3f axis_2 = Vector3f(0, 0, z_sign) 
                - Vector3f(0.5*x_sign, 0.5*y_sign, 0);
        Constants::TILT_PRINCIPLE_AXIS[i] = Eigen::Matrix3d::Zero();
        Constants::TILT_PRINCIPLE_AXIS[i].col(0) = axis.normalized();
        Constants::TILT_PRINCIPLE_AXIS[i].col(1) = axis_1.normalized();
        Constants::TILT_PRINCIPLE_AXIS[i].col(2) = axis_2.normalized();
    }
}