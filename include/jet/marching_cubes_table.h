// Copyright (c) 2018 Doyub Kim
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.
//
// Marching Cubes Example Program
// by Cory Bloyd (corysama@yahoo.com)
//
// A simple, portable and complete implementation of the Marching Cubes
// and Marching Tetrahedrons algorithms in a single source file.
// There are many ways that this code could be made faster, but the
// intent is for the code to be easy to understand.
//
// For a description of the algorithm go to
// http://astronomy.swin.edu.au/pbourke/modelling/polygonise/
//
// This code is public domain.
//

#ifndef SRC_JET_MARCHING_CUBES_TABLE_H_
#define SRC_JET_MARCHING_CUBES_TABLE_H_

namespace jet {

// vertexOffset lists the positions, relative to vertex0, of each of the 8
// vertices of a cube
static const double vertexOffset[8][3] = {
    {0.0f, 0.0f, 0.0f}, {1.0f, 0.0f, 0.0f},
    {1.0f, 0.0f, 1.0f}, {0.0f, 0.0f, 1.0f},
    {0.0f, 1.0f, 0.0f}, {1.0f, 1.0f, 0.0f},
    {1.0f, 1.0f, 1.0f}, {0.0f, 1.0f, 1.0f}
};

// edgeConnection lists the index of the endpoint vertices for each of the 12
// edges of the cube
static const int edgeConnection[12][2] = {
    {0, 1}, {1, 2}, {3, 2}, {0, 3},
    {4, 5}, {5, 6}, {7, 6}, {4, 7},
    {0, 4}, {1, 5}, {2, 6}, {3, 7}
};

// edgeDirection lists the direction vector (vertex[i+1]-vertex[i])/2 for each
// edge in the cube
static const double edgeDirection[12][3] = {
    {1.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 1.0f},
    {1.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 1.0f},
    {1.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 1.0f},
    {1.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 1.0f},
    {0.0f, 1.0f, 0.0f}, {0.0f, 1.0f, 0.0f},
    {0.0f, 1.0f, 0.0f}, {0.0f, 1.0f, 0.0f}
};

// edgeCenter lists for calculating counterclockwise order
static const double edgeCenter[12][3] = {
    {0.5f, 0.0f, 0.0f}, {1.0f, 0.0f, 0.5f},
    {0.5f, 0.0f, 1.0f}, {0.0f, 0.0f, 0.5f},
    {0.5f, 1.0f, 0.0f}, {1.0f, 1.0f, 0.5f},
    {0.5f, 1.0f, 1.0f}, {0.0f, 1.0f, 0.5f},
    {0.0f, 0.5f, 0.0f}, {1.0f, 0.5f, 0.0f},
    {1.0f, 0.5f, 1.0f}, {0.0f, 0.5f, 1.0f}
};

static const int threeNeighborNodes[8][3] = {
    {1, 3, 4}, {0, 2, 5}, {1, 3, 6}, {0, 2, 7},
    {0, 5, 7}, {1, 4, 6}, {2, 5, 7}, {3, 4, 6}
};

static const int threeNeighborEdges[8][3] = {
    {0, 3, 8}, {0, 1, 9}, {1, 2, 10}, {2, 3, 11},
    {4, 7, 8}, {4, 5, 9}, {5, 6, 10}, {6, 7, 11}
};

// For any edge, if one vertex is inside of the surface and the other is outside
// of the surface then the edge intersects the surface
// For each of the 8 vertices of the cube can be two possible states : either
// inside or outside of the surface
// For any cube there are 2^8=256 possible sets of vertex states
// This table lists the edges intersected by the surface for all 256 possible
// vertex states
// There are 12 edges.  For each entry in the table, if edge #n is intersected,
// then bit #n is set to 1
static const int cubeEdgeFlags[256] = {
    0x000, 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
    0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
    0x190, 0x099, 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
    0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
    0x230, 0x339, 0x033, 0x13a, 0x636, 0x73f, 0x435, 0x53c,
    0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
    0x3a0, 0x2a9, 0x1a3, 0x0aa, 0x7a6, 0x6af, 0x5a5, 0x4ac,
    0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
    0x460, 0x569, 0x663, 0x76a, 0x066, 0x16f, 0x265, 0x36c,
    0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
    0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0x0ff, 0x3f5, 0x2fc,
    0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
    0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x055, 0x15c,
    0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
    0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0x0cc,
    0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
    0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
    0x0cc, 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
    0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
    0x15c, 0x055, 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
    0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
    0x2fc, 0x3f5, 0x0ff, 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
    0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
    0x36c, 0x265, 0x16f, 0x066, 0x76a, 0x663, 0x569, 0x460,
    0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
    0x4ac, 0x5a5, 0x6af, 0x7a6, 0x0aa, 0x1a3, 0x2a9, 0x3a0,
    0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
    0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x033, 0x339, 0x230,
    0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
    0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x099, 0x190,
    0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
    0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x000
};

// For each of the possible vertex states listed in CubeEdgeFlags
// there is a specific triangulation of the edge intersection points.
// TriangleConnectionTable lists all of them in the form of 0-5 edge triples
// with the list terminated by the invalid value -1.
// For example: TriangleConnectionTable[3] list the 2 triangles formed
// when corner[0] and corner[1] are inside of the surface, but the rest of the
// cube is not.
// three vertices of a triangle = counter-clock-wise order
static const int triangleConnectionTable3D[256][16] = {                 // EdgeFlag, NodeFlag   special case for MDC
    { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0x000,   0000 0000   X
    {  0,  8,  3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0x109,   0000 0001   O   {  0, 8,  8, 3,  3, 0
    {  0,  1,  9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0x203,   0000 0010   O   {  0, 1,  1, 9,  9, 0
    {  1,  8,  3,  9,  8,  1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0x30a,   0000 0011   O   {  1, 9,  9, 8,  8, 3,  3, 1
    {  1,  2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0x406,   0000 0100   O   {  1, 2,  2, 10, 10, 1
    {  0,  8,  3,  1,  2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0x50f,   0000 0101   X
    {  9,  2, 10,  0,  2,  9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0x605,   0000 0110   O   {  0, 2,  2, 10, 10, 9,  9, 0
    {  2,  8,  3,  2, 10,  8, 10,  9,  8, -1, -1, -1, -1, -1, -1, -1 }, // 0x70c,   0000 0111   O   {  3, 2,  2, 10, 10, 9,  9, 8,  8, 3
    {  3, 11,  2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0x80c,   0000 1000   O   {  3, 11, 11, 2,  2, 3
    {  0, 11,  2,  8, 11,  0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0x905,   0000 1001   O   {  2, 0,  0, 8,  8, 11, 11, 2
    {  1,  9,  0,  2,  3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0xa0f,   0000 1010   X
    {  1, 11,  2,  1,  9, 11,  9,  8, 11, -1, -1, -1, -1, -1, -1, -1 }, // 0xb06,   0000 1011   O   {  2, 1,  1, 9,  9, 8,  8, 11, 11, 2
    {  3, 10,  1, 11, 10,  3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0xc0a,   0000 1100   O   {  1, 3,  3, 11, 11, 10, 10, 1
    {  0, 10,  1,  0,  8, 10,  8, 11, 10, -1, -1, -1, -1, -1, -1, -1 }, // 0xd03,   0000 1101   O   {  1, 0,  0, 8,  8, 11, 11, 10, 10, 1
    {  3,  9,  0,  3, 11,  9, 11, 10,  9, -1, -1, -1, -1, -1, -1, -1 }, // 0xe09,   0000 1110   O   {  0, 3,  3, 11, 11, 10, 10, 9,  9, 0
    {  9,  8, 10, 10,  8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0xf00,   0000 1111   O   {  9, 8,  8, 11, 11, 10, 10, 9

    {  4,  7,  8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0x190,   0001 0000   O   {  4, 7,  7, 8,  8, 4
    {  4,  3,  0,  7,  3,  4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0x099,   0001 0001   O   {  3, 0,  0, 4,  4, 7,  7, 3
    {  0,  1,  9,  8,  4,  7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0x393,   0001 0010   X
    {  4,  1,  9,  4,  7,  1,  7,  3,  1, -1, -1, -1, -1, -1, -1, -1 }, // 0x29a,   0001 0011   O   {  3, 1,  1, 9,  9, 4,  4, 7,  7, 3
    {  1,  2, 10,  8,  4,  7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0x596,   0001 0100   X
    {  3,  4,  7,  3,  0,  4,  1,  2, 10, -1, -1, -1, -1, -1, -1, -1 }, // 0x49f,   0001 0101   X
    {  9,  2, 10,  9,  0,  2,  8,  4,  7, -1, -1, -1, -1, -1, -1, -1 }, // 0x795,   0001 0110   X
    {  2, 10,  9,  2,  9,  7,  2,  7,  3,  7,  9,  4, -1, -1, -1, -1 }, // 0x69c,   0001 0111   O   {  3, 2,  2, 10, 10, 9,  9, 4,  4, 7,  7, 3
    {  8,  4,  7,  3, 11,  2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0x99c,   0001 1000   X
    { 11,  4,  7, 11,  2,  4,  2,  0,  4, -1, -1, -1, -1, -1, -1, -1 }, // 0x895,   0001 1001   O   {  2, 0,  0, 4,  4, 7,  7, 11, 11, 2
    {  9,  0,  1,  8,  4,  7,  2,  3, 11, -1, -1, -1, -1, -1, -1, -1 }, // 0xb9f,   0001 1010   X
    {  4,  7, 11,  9,  4, 11,  9, 11,  2,  9,  2,  1, -1, -1, -1, -1 }, // 0xa96,   0001 1011   O   {  2, 1,  1, 9,  9, 4,  4, 7,  7, 11, 11, 2
    {  3, 10,  1,  3, 11, 10,  7,  8,  4, -1, -1, -1, -1, -1, -1, -1 }, // 0xd9a,   0001 1100   X
    {  1, 11, 10,  1,  4, 11,  1,  0,  4,  7, 11,  4, -1, -1, -1, -1 }, // 0xc93,   0001 1101   O   {  1, 0,  0, 4,  4, 7,  7, 11, 11, 10, 10, 1
    {  4,  7,  8,  9,  0, 11,  9, 11, 10, 11,  0,  3, -1, -1, -1, -1 }, // 0xf99,   0001 1110   X
    {  4,  7, 11,  4, 11,  9,  9, 11, 10, -1, -1, -1, -1, -1, -1, -1 }, // 0xe90,   0001 1111   O   {  4, 7,  7, 11, 11, 10, 10, 9,  9, 4,

    {  9,  5,  4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0x230,   0010 0000   O   {  9, 5,  5, 4,  4, 9
    {  9,  5,  4,  0,  8,  3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0x339,   0010 0001   X
    {  0,  5,  4,  1,  5,  0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0x033,   0010 0010   O   {  0, 1,  1, 5,  5, 4,  4, 0
    {  8,  5,  4,  8,  3,  5,  3,  1,  5, -1, -1, -1, -1, -1, -1, -1 }, // 0x13a,   0010 0011   O   {  3, 1,  1, 5,  5, 4,  4, 8,  8, 3
    {  1,  2, 10,  9,  5,  4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0x636,   0010 0100   X
    {  3,  0,  8,  1,  2, 10,  4,  9,  5, -1, -1, -1, -1, -1, -1, -1 }, // 0x73f,   0010 0101   X
    {  5,  2, 10,  5,  4,  2,  4,  0,  2, -1, -1, -1, -1, -1, -1, -1 }, // 0x435,   0010 0110   O   {  0, 2,  2, 10, 10, 5,  5, 4,  4, 0
    {  2, 10,  5,  3,  2,  5,  3,  5,  4,  3,  4,  8, -1, -1, -1, -1 }, // 0x53c,   0010 0111   O   {  3, 2,  2, 10, 10, 5,  5, 4,  4, 8,  8, 3
    {  9,  5,  4,  2,  3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0xa3c,   0010 1000   X
    {  0, 11,  2,  0,  8, 11,  4,  9,  5, -1, -1, -1, -1, -1, -1, -1 }, // 0xb35,   0010 1001   X
    {  0,  5,  4,  0,  1,  5,  2,  3, 11, -1, -1, -1, -1, -1, -1, -1 }, // 0x83f,   0010 1010   X
    {  2,  1,  5,  2,  5,  8,  2,  8, 11,  4,  8,  5, -1, -1, -1, -1 }, // 0x936,   0010 1011   O   {  2, 1,  1, 5,  5, 4,  4, 8,  8, 11, 11, 2
    { 10,  3, 11, 10,  1,  3,  9,  5,  4, -1, -1, -1, -1, -1, -1, -1 }, // 0xe3a,   0010 1100   X
    {  4,  9,  5,  0,  8,  1,  8, 10,  1,  8, 11, 10, -1, -1, -1, -1 }, // 0xf33,   0010 1101   X
    {  5,  4,  0,  5,  0, 11,  5, 11, 10, 11,  0,  3, -1, -1, -1, -1 }, // 0xc39,   0010 1110   O   {  0, 3,  3, 11, 11, 10, 10, 5,  5, 4,  4, 0
    {  5,  4,  8,  5,  8, 10, 10,  8, 11, -1, -1, -1, -1, -1, -1, -1 }, // 0xd30,   0010 1111   O   {  5, 4,  4, 8,  8, 11, 11, 10, 10, 5

    {  9,  7,  8,  5,  7,  9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0x3a0,   0011 0000   O   {  5, 7,  7, 8,  8, 9,  9, 5
    {  9,  3,  0,  9,  5,  3,  5,  7,  3, -1, -1, -1, -1, -1, -1, -1 }, // 0x2a9,   0011 0001   O   {  3, 0,  0, 9,  9, 5,  5, 7,  7, 3
    {  0,  7,  8,  0,  1,  7,  1,  5,  7, -1, -1, -1, -1, -1, -1, -1 }, // 0x1a3,   0011 0010   O   {  0, 1,  1, 5,  5, 7,  7, 8,  8, 0
    {  1,  5,  3,  3,  5,  7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0x0aa,   0011 0011   O   {  3, 1,  1, 5,  5, 7,  7, 3
    {  9,  7,  8,  9,  5,  7, 10,  1,  2, -1, -1, -1, -1, -1, -1, -1 }, // 0x7a6,   0011 0100   X
    { 10,  1,  2,  9,  5,  0,  5,  3,  0,  5,  7,  3, -1, -1, -1, -1 }, // 0x6af,   0011 0101   X
    {  8,  0,  2,  8,  2,  5,  8,  5,  7, 10,  5,  2, -1, -1, -1, -1 }, // 0x5a5,   0011 0110   O   {  0, 2,  2, 10, 10, 5,  5, 7,  7, 8,  8, 0
    {  2, 10,  5,  2,  5,  3,  3,  5,  7, -1, -1, -1, -1, -1, -1, -1 }, // 0x4ac,   0011 0111   O   {  3, 2,  2, 10, 10, 5,  5, 7,  7, 3
    {  7,  9,  5,  7,  8,  9,  3, 11,  2, -1, -1, -1, -1, -1, -1, -1 }, // 0xbac,   0011 1000   X
    {  9,  5,  7,  9,  7,  2,  9,  2,  0,  2,  7, 11, -1, -1, -1, -1 }, // 0xaa5,   0011 1001   O   {  2, 0,  0, 9,  9, 5,  5, 7,  7, 11, 11, 2
    {  2,  3, 11,  0,  1,  8,  1,  7,  8,  1,  5,  7, -1, -1, -1, -1 }, // 0x9af,   0011 1010   X
    { 11,  2,  1, 11,  1,  7,  7,  1,  5, -1, -1, -1, -1, -1, -1, -1 }, // 0x8a6,   0011 1011   O   {  2, 1,  1, 5,  5, 7,  7, 11, 11, 2
    {  9,  5,  8,  8,  5,  7, 10,  1,  3, 10,  3, 11, -1, -1, -1, -1 }, // 0xfaa,   0011 1100   X
    {  5,  7,  0,  5,  0,  9,  7, 11,  0,  1,  0, 10, 11, 10,  0, -1 }, // 0xea3,   0011 1101   X
    { 11, 10,  0, 11,  0,  3, 10,  5,  0,  8,  0,  7,  5,  7,  0, -1 }, // 0xda9,   0011 1110   X
    { 11, 10,  5,  7, 11,  5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0xca0,   0011 1111   O   {  5, 7,  7, 11, 11, 10, 10, 5

    { 10,  6,  5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0x460,   0100 0000   O   { 10, 6,  6, 5,  5, 10
    {  0,  8,  3,  5, 10,  6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0x569,   0100 0001   X
    {  9,  0,  1,  5, 10,  6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0x663,   0100 0010   X
    {  1,  8,  3,  1,  9,  8,  5, 10,  6, -1, -1, -1, -1, -1, -1, -1 }, // 0x76a,   0100 0011   X
    {  1,  6,  5,  2,  6,  1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0x066,   0100 0100   O   {  1, 2,  2, 6,  6, 5,  5, 1
    {  1,  6,  5,  1,  2,  6,  3,  0,  8, -1, -1, -1, -1, -1, -1, -1 }, // 0x16f,   0100 0101   X
    {  9,  6,  5,  9,  0,  6,  0,  2,  6, -1, -1, -1, -1, -1, -1, -1 }, // 0x265,   0100 0110   O   {  0, 2,  2, 6,  6, 5,  5, 9,  9, 0
    {  5,  9,  8,  5,  8,  2,  5,  2,  6,  3,  2,  8, -1, -1, -1, -1 }, // 0x36c,   0100 0111   O   {  3, 2,  2, 6,  6, 5,  5, 9,  9, 8,  8, 3
    {  2,  3, 11, 10,  6,  5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0xc6c,   0100 1000   X
    { 11,  0,  8, 11,  2,  0, 10,  6,  5, -1, -1, -1, -1, -1, -1, -1 }, // 0xd65,   0100 1001   X
    {  0,  1,  9,  2,  3, 11,  5, 10,  6, -1, -1, -1, -1, -1, -1, -1 }, // 0xe6f,   0100 1010   X
    {  5, 10,  6,  1,  9,  2,  9, 11,  2,  9,  8, 11, -1, -1, -1, -1 }, // 0xf66,   0100 1011   X
    {  6,  3, 11,  6,  5,  3,  5,  1,  3, -1, -1, -1, -1, -1, -1, -1 }, // 0x86a,   0100 1100   O   {  5, 1,  1, 3,  3, 11, 11, 6,  6, 5
    {  0,  8, 11,  0, 11,  5,  0,  5,  1,  5, 11,  6, -1, -1, -1, -1 }, // 0x963,   0100 1101   O   {  1, 0,  0, 8,  8, 11, 11, 6,  6, 5,  5, 1
    {  3, 11,  6,  0,  3,  6,  0,  6,  5,  0,  5,  9, -1, -1, -1, -1 }, // 0xa69,   0100 1110   O   {  0, 3,  3, 11, 11, 6,  6, 5,  5, 9,  9, 0
    {  6,  5,  9,  6,  9, 11, 11,  9,  8, -1, -1, -1, -1, -1, -1, -1 }, // 0xb60,   0100 1111   O   {  6, 5,  5, 9,  9, 8,  8, 11, 11, 6

    {  5, 10,  6,  4,  7,  8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0x5f0,   0101 0000   X
    {  4,  3,  0,  4,  7,  3,  6,  5, 10, -1, -1, -1, -1, -1, -1, -1 }, // 0x4f9,   0101 0001   X
    {  1,  9,  0,  5, 10,  6,  8,  4,  7, -1, -1, -1, -1, -1, -1, -1 }, // 0x7f3,   0101 0010   X
    { 10,  6,  5,  1,  9,  7,  1,  7,  3,  7,  9,  4, -1, -1, -1, -1 }, // 0x6fa,   0101 0011   X
    {  6,  1,  2,  6,  5,  1,  4,  7,  8, -1, -1, -1, -1, -1, -1, -1 }, // 0x1f6,   0101 0100   X
    {  1,  2,  5,  5,  2,  6,  3,  0,  4,  3,  4,  7, -1, -1, -1, -1 }, // 0x0ff,   0101 0101   X
    {  8,  4,  7,  9,  0,  5,  0,  6,  5,  0,  2,  6, -1, -1, -1, -1 }, // 0x3f5,   0101 0110   X
    {  7,  3,  9,  7,  9,  4,  3,  2,  9,  5,  9,  6,  2,  6,  9, -1 }, // 0xxfc,   0101 0111   X
    {  3, 11,  2,  7,  8,  4, 10,  6,  5, -1, -1, -1, -1, -1, -1, -1 }, // 0xdfc,   0101 1000   X
    {  5, 10,  6,  4,  7,  2,  4,  2,  0,  2,  7, 11, -1, -1, -1, -1 }, // 0xcf5,   0101 1001   X
    {  0,  1,  9,  4,  7,  8,  2,  3, 11,  5, 10,  6, -1, -1, -1, -1 }, // 0xfff,   0101 1010   X
    {  9,  2,  1,  9, 11,  2,  9,  4, 11,  7, 11,  4,  5, 10,  6, -1 }, // 0xef6,   0101 1011   X
    {  8,  4,  7,  3, 11,  5,  3,  5,  1,  5, 11,  6, -1, -1, -1, -1 }, // 0x9fa,   0101 1100   X
    {  5,  1, 11,  5, 11,  6,  1,  0, 11,  7, 11,  4,  0,  4, 11, -1 }, // 0x8f3,   0101 1101   X
    {  0,  5,  9,  0,  6,  5,  0,  3,  6, 11,  6,  3,  8,  4,  7, -1 }, // 0xbf9,   0101 1110   X
    {  6,  5,  9,  6,  9, 11,  4,  7,  9,  7, 11,  9, -1, -1, -1, -1 }, // 0xaf0,   0101 1111   O   {  4, 7,  7, 11, 11, 6,  6, 5,  5, 9,  9, 4

    { 10,  4,  9,  6,  4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0x650,   0110 0000   O   {  6, 4,  4, 9,  9, 10, 10, 6
    {  4, 10,  6,  4,  9, 10,  0,  8,  3, -1, -1, -1, -1, -1, -1, -1 }, // 0x759,   0110 0001   X
    { 10,  0,  1, 10,  6,  0,  6,  4,  0, -1, -1, -1, -1, -1, -1, -1 }, // 0x453,   0110 0010   O   {  0, 1,  1, 10, 10, 6,  6, 4,  4, 0
    {  8,  3,  1,  8,  1,  6,  8,  6,  4,  6,  1, 10, -1, -1, -1, -1 }, // 0x55a,   0110 0011   O   {  3, 1,  1, 10, 10, 6,  6, 4,  4, 8,  8, 3
    {  1,  4,  9,  1,  2,  4,  2,  6,  4, -1, -1, -1, -1, -1, -1, -1 }, // 0x256,   0110 0100   O   {  1, 2,  2, 6,  6, 4,  4, 9,  9, 1
    {  3,  0,  8,  1,  2,  9,  2,  4,  9,  2,  6,  4, -1, -1, -1, -1 }, // 0x35f,   0110 0101   X
    {  0,  2,  4,  4,  2,  6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0x055,   0110 0110   O   {  0, 2,  2, 6,  6, 4,  4, 0
    {  8,  3,  2,  8,  2,  4,  4,  2,  6, -1, -1, -1, -1, -1, -1, -1 }, // 0x15c,   0110 0111   O   {  3, 2,  2, 6,  6, 4,  4, 8,  8, 3
    { 10,  4,  9, 10,  6,  4, 11,  2,  3, -1, -1, -1, -1, -1, -1, -1 }, // 0xe5c,   0110 1000   X
    {  0,  8,  2,  2,  8, 11,  4,  9, 10,  4, 10,  6, -1, -1, -1, -1 }, // 0xf55,   0110 1001   X
    {  3, 11,  2,  0,  1,  6,  0,  6,  4,  6,  1, 10, -1, -1, -1, -1 }, // 0xc5f,   0110 1010   X
    {  6 , 4,  1,  6,  1, 10,  4,  8,  1,  2,  1, 11,  8, 11,  1, -1 }, // 0xd56,   0110 1011   O   {  2, 1,  1, 10, 10, 6,  6, 4,  4, 8,  8, 11, 11, 2
    {  9,  6,  4,  9,  3,  6,  9,  1,  3, 11,  6,  3, -1, -1, -1, -1 }, // 0xa5a,   0110 1100   O   {  1, 3,  3, 11, 11, 6,  6, 4,  4, 9,  9, 1
    {  8, 11,  1,  8,  1,  0, 11,  6,  1,  9,  1,  4,  6,  4,  1, -1 }, // 0xb53,   0110 1101   X
    {  3, 11,  6,  3,  6,  0,  0,  6,  4, -1, -1, -1, -1, -1, -1, -1 }, // 0x859,   0110 1110   O   {  0, 3,  3, 11, 11, 6,  6, 4,  4, 0
    {  6,  4,  8, 11,  6,  8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0x950,   0110 1111   O   {  6, 4,  4, 8,  8, 11, 11, 6

    {  7, 10,  6,  7,  8, 10,  8,  9, 10, -1, -1, -1, -1, -1, -1, -1 }, // 0x7c0,   0111 0000   O   {  6, 7,  7, 8,  8, 9,  9, 10, 10, 6
    {  0,  7,  3,  0, 10,  7,  0,  9, 10,  6,  7, 10, -1, -1, -1, -1 }, // 0x6c9,   0111 0001   O   {  0, 10, 10, 6,  6, 7,  7, 3,  3, 0
    { 10,  6,  7,  1, 10,  7,  1,  7,  8,  1,  8,  0, -1, -1, -1, -1 }, // 0x5c3,   0111 0010   O   {  0, 1,  1, 10, 10, 6,  6, 7,  7, 8,  8, 0
    { 10,  6,  7, 10,  7,  1,  1,  7,  3, -1, -1, -1, -1, -1, -1, -1 }, // 0x4ca,   0111 0011   O   {  3, 1,  1, 10, 10, 6,  6, 7,  7, 3
    {  1,  2,  6,  1,  6,  8,  1,  8,  9,  8,  6,  7, -1, -1, -1, -1 }, // 0x3c6,   0111 0100   O   {  1, 2,  2, 6,  6, 7,  7, 8,  8, 9,  9, 1
    {  2,  6,  9,  2,  9,  1,  6,  7,  9,  0,  9,  3,  7,  3,  9, -1 }, // 0x2cf,   0111 0101   X
    {  7,  8,  0,  7,  0,  6,  6,  0,  2, -1, -1, -1, -1, -1, -1, -1 }, // 0x1c5,   0111 0110   O   {  0, 2,  2, 6,  6, 7,  7, 8,  8, 0
    {  7,  3,  2,  6,  7,  2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0x0cc,   0111 0111   O   {  3, 2,  2, 6,  6, 7,  7, 3
    {  2,  3, 11, 10,  6,  8, 10,  8,  9,  8,  6,  7, -1, -1, -1, -1 }, // 0xfcc,   0111 1000   X
    {  2,  0,  7,  2,  7, 11,  0,  9,  7,  6,  7, 10,  9, 10,  7, -1 }, // 0xec5,   0111 1001   X
    {  1,  8,  0,  1,  7,  8,  1, 10,  7,  6,  7, 10,  2,  3, 11, -1 }, // 0xdcf,   0111 1010   X
    { 11,  2,  1, 11,  1,  7, 10,  6,  1,  6,  7,  1, -1, -1, -1, -1 }, // 0xcc6,   0111 1011   O   {  2, 1,  1, 10, 10, 6,  6, 7,  7, 11, 11, 2
    {  8,  9,  6,  8,  6,  7,  9,  1,  6, 11,  6,  3,  1,  3 , 6, -1 }, // 0xbca,   0111 1100   X
    {  0,  9,  1, 11,  6,  7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0xac3,   0111 1101   X
    {  7,  8,  0,  7,  0,  6,  3, 11,  0, 11,  6,  0, -1, -1, -1, -1 }, // 0x9c9,   0111 1110   O   {  0, 3,  3, 11, 11, 6,  6, 7,  7, 8,  8, 0
    {  7, 11,  6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0x8c0,   0111 1111   O   {  6, 7,  7, 11, 11, 6

    {  7,  6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0x8c0,   1000 0000   O   {  7, 6,  6, 11, 11, 7
    {  3,  0,  8, 11,  7,  6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0x9c9,   1000 0001   X
    {  0,  1,  9, 11,  7,  6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0xac3,   1000 0010   X
    {  8,  1,  9,  8,  3,  1, 11,  7,  6, -1, -1, -1, -1, -1, -1, -1 }, // 0xbca,   1000 0011   X
    { 10,  1,  2,  6, 11,  7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0xcc6,   1000 0100   X
    {  1,  2, 10,  3,  0,  8,  6, 11,  7, -1, -1, -1, -1, -1, -1, -1 }, // 0xdcf,   1000 0101   X
    {  2,  9,  0,  2, 10,  9,  6, 11,  7, -1, -1, -1, -1, -1, -1, -1 }, // 0xec5,   1000 0110   X
    {  6, 11,  7,  2, 10,  3, 10,  8,  3, 10,  9,  8, -1, -1, -1, -1 }, // 0xfcc,   1000 0111   X
    {  7,  2,  3,  6,  2,  7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0x0cc,   1000 1000   O   {  2, 3,  3, 7,  7, 6,  6, 2
    {  7,  0,  8,  7,  6,  0,  6,  2,  0, -1, -1, -1, -1, -1, -1, -1 }, // 0x1c5,   1000 1001   O   {  2, 0,  0, 8,  8, 7,  7, 6,  6, 2
    {  2,  7,  6,  2,  3,  7,  0,  1,  9, -1, -1, -1, -1, -1, -1, -1 }, // 0x2cf,   1000 1010   X
    {  1,  6,  2,  1,  8,  6,  1,  9,  8,  8,  7,  6, -1, -1, -1, -1 }, // 0x3c6,   1000 1011   O   {  2, 1,  1, 9,  9, 8,  8, 7,  7, 6,  6, 2
    { 10,  7,  6, 10,  1,  7,  1,  3,  7, -1, -1, -1, -1, -1, -1, -1 }, // 0x4ca,   1000 1100   O   {  1, 3,  3, 7,  7, 6,  6, 10, 10, 1
    { 10,  7,  6,  1,  7, 10,  1,  8,  7,  1,  0,  8, -1, -1, -1, -1 }, // 0x5c3,   1000 1101   O   {  1, 0,  0, 8,  8, 7,  7, 6,  6, 10, 10, 1
    {  0,  3,  7,  0,  7, 10,  0, 10,  9,  6, 10,  7, -1, -1, -1, -1 }, // 0x6c9,   1000 1110   O   {  0, 3,  3, 7,  7, 6,  6, 10, 10, 9,  9, 0
    {  7,  6, 10,  7, 10,  8,  8, 10,  9, -1, -1, -1, -1, -1, -1, -1 }, // 0x7c0,   1000 1111   O   {  7, 6,  6, 10, 10, 9,  9, 8,  8, 7

    {  6,  8,  4, 11,  8,  6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0x950,   1001 0000   O   {  4, 6,  6, 11, 11, 8,  8, 4
    {  3,  6, 11,  3,  0,  6,  0,  4,  6, -1, -1, -1, -1, -1, -1, -1 }, // 0x859,   1001 0001   O   {  3, 0,  0, 4,  4, 6,  6, 11, 11, 3
    {  8,  6, 11,  8,  4,  6,  9,  0,  1, -1, -1, -1, -1, -1, -1, -1 }, // 0xb53,   1001 0010   X
    {  9,  4,  6,  9,  6,  3,  9,  3,  1, 11,  3,  6, -1, -1, -1, -1 }, // 0xa5a,   1001 0011   O   {  1, 9,  9, 4,  4, 6,  6, 11, 11, 3,  3, 1
    {  6,  8,  4,  6, 11,  8,  2, 10,  1, -1, -1, -1, -1, -1, -1, -1 }, // 0xd56,   1001 0100   X
    {  1,  2, 10,  3,  0, 11,  0,  6, 11,  0,  4,  6, -1, -1, -1, -1 }, // 0xc5f,   1001 0101   X
    {  4, 11,  8,  4,  6, 11,  0,  2,  9,  2, 10,  9, -1, -1, -1, -1 }, // 0xf55,   1001 0110   X
    { 10,  9,  3, 10,  3,  2,  9,  4,  3, 11,  3,  6,  4,  6,  3, -1 }, // 0xe5c,   1001 0111   X
    {  8,  2,  3,  8,  4,  2,  4,  6,  2, -1, -1, -1, -1, -1, -1, -1 }, // 0x15c,   1001 1000   O   {  2, 3,  3, 8,  8, 4,  4, 6,  6, 2
    {  0,  4,  2,  4,  6,  2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0x055,   1001 1001   O   {  2, 0,  0, 4,  4, 6,  6, 2
    {  1,  9,  0,  2,  3,  4,  2,  4,  6,  4,  3,  8, -1, -1, -1, -1 }, // 0x35f,   1001 1010   X
    {  1,  9,  4,  1,  4,  2,  2,  4,  6, -1, -1, -1, -1, -1, -1, -1 }, // 0x256,   1001 1011   O   {  2, 1,  1, 9,  9, 4,  4, 6,  6, 2
    {  8,  1,  3,  8,  6,  1,  8,  4,  6,  6, 10,  1, -1, -1, -1, -1 }, // 0x55a,   1001 1100   O   {  1, 3,  3, 8,  8, 4,  4, 6,  6, 10, 10, 1
    { 10,  1,  0, 10,  0,  6,  6,  0,  4, -1, -1, -1, -1, -1, -1, -1 }, // 0x453,   1001 1101   O   {  1, 0,  0, 4,  4, 6,  6, 10, 10, 1
    {  4,  6,  3,  4,  3,  8,  6, 10,  3,  0,  3,  9, 10,  9,  3, -1 }, // 0x759,   1001 1110   X
    { 10,  9,  4,  6, 10,  4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0x650,   1001 1111   O   {  4, 6,  6, 10, 10, 9,  9, 4

    {  4,  9,  5,  7,  6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0xaf0,   1010 0000   X
    {  0,  8,  3,  4,  9,  5, 11,  7,  6, -1, -1, -1, -1, -1, -1, -1 }, // 0xbf9,   1010 0001   X
    {  5,  0,  1,  5,  4,  0,  7,  6, 11, -1, -1, -1, -1, -1, -1, -1 }, // 0x8f3,   1010 0010   X
    { 11,  7,  6,  8,  3,  4,  3,  5,  4,  3,  1,  5, -1, -1, -1, -1 }, // 0x9fa,   1010 0011   X
    {  9,  5,  4, 10,  1,  2,  7,  6, 11, -1, -1, -1, -1, -1, -1, -1 }, // 0xef6,   1010 0100   X
    {  6, 11,  7,  1,  2, 10,  0,  8,  3,  4,  9,  5, -1, -1, -1, -1 }, // 0xfff,   1010 0101   X
    {  7,  6, 11,  5,  4, 10,  4,  2, 10,  4,  0,  2, -1, -1, -1, -1 }, // 0xcf5,   1010 0110   X
    {  3,  4,  8,  3,  5,  4,  3,  2,  5, 10,  5,  2, 11,  7,  6, -1 }, // 0xdfc,   1010 0111   X
    {  7,  2,  3,  7,  6,  2,  5,  4,  9, -1, -1, -1, -1, -1, -1, -1 }, // 0x2fc,   1010 1000   X
    {  9,  5,  4,  0,  8,  6,  0,  6,  2,  6,  8,  7, -1, -1, -1, -1 }, // 0x3f5,   1010 1001   X
    {  3,  6,  2,  3,  7,  6,  1,  5,  0,  5,  4,  0, -1, -1, -1, -1 }, // 0x0ff,   1010 1010   X
    {  6,  2,  8,  6,  8,  7,  2,  1,  8,  4,  8,  5,  1,  5,  8, -1 }, // 0x1f6,   1010 1011   X
    {  9,  5,  4, 10,  1,  6,  1,  7,  6,  1,  3,  7, -1, -1, -1, -1 }, // 0x6fa,   1010 1100   X
    {  1,  6, 10,  1,  7,  6,  1,  0,  7,  8,  7,  0,  9,  5,  4, -1 }, // 0x7f3,   1010 1101   X
    {  4,  0, 10,  4, 10,  5,  0,  3, 10,  6, 10,  7,  3,  7, 10, -1 }, // 0x4f9,   1010 1110   X
    {  7,  6, 10,  7, 10,  8,  5,  4, 10,  4,  8, 10, -1, -1, -1, -1 }, // 0x5f0,   1010 1111   O   {  5, 4,  4, 8,  8, 7,  7, 6,  6, 10, 10, 5

    {  6,  9,  5,  6, 11,  9, 11,  8,  9, -1, -1, -1, -1, -1, -1, -1 }, // 0xb60,   1011 0000   O   {  5, 6,  6, 11, 11, 8,  8, 9,  9, 5
    {  3,  6, 11,  0,  6,  3,  0,  5,  6,  0,  9,  5, -1, -1, -1, -1 }, // 0xa69,   1011 0001   O   {  5, 6,  6, 11, 11, 3,  3, 0,  0, 9,  9, 5
    {  0, 11,  8,  0,  5, 11,  0,  1,  5,  5,  6, 11, -1, -1, -1, -1 }, // 0x963,   1011 0010   O   {  0, 1,  1, 5,  5, 6,  6, 11, 11, 8,  8, 0
    {  6, 11,  3,  6,  3,  5,  5,  3 , 1, -1, -1, -1, -1, -1, -1, -1 }, // 0x86a,   1011 0011   O   {  3, 1,  1, 5,  5, 6,  6, 11, 11, 3
    {  1,  2, 10,  9,  5, 11,  9, 11,  8, 11,  5,  6, -1, -1, -1, -1 }, // 0xf66,   1011 0100   X
    {  0, 11,  3,  0,  6, 11,  0,  9,  6,  5,  6,  9,  1,  2, 10, -1 }, // 0xe6f,   1011 0101   X
    { 11,  8,  5, 11,  5,  6,  8,  0,  5, 10,  5,  2,  0,  2,  5, -1 }, // 0xd65,   1011 0110   X
    {  3,  6, 11,  6,  3,  5,  2, 10,  3, 10,  5,  3, -1, -1, -1, -1 }, //{  2, 11,  3,  5,  6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0xc6c,  1011 0111   X
    {  5,  8,  9,  5,  2,  8,  5,  6,  2,  3,  8,  2, -1, -1, -1, -1 }, // 0x36c,   1011 1000   O   {  2, 3,  3, 8,  8, 9,  9, 5,  5, 6,  6, 2
    {  9,  5,  6,  9,  6,  0,  0,  6,  2, -1, -1, -1, -1, -1, -1, -1 }, // 0x265,   1011 1001   O   {  2, 0,  0, 9,  9, 5,  5, 6,  6, 2
    {  1,  5,  8,  1,  8,  0,  5,  6,  8,  3,  8,  2,  6,  2,  8, -1 }, // 0x16f,   1011 1010   X
    {  1,  5,  6,  2,  1,  6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0x066,   1011 1011   O   {  2, 1,  1, 5,  5, 6,  6, 2
    {  1,  3,  6,  1,  6, 10,  3,  8,  6,  5,  6,  9,  8,  9,  6, -1 }, // 0x76a,   1011 1100   X
    { 10,  1,  0, 10,  0,  6,  9,  5,  0,  5,  6, 0 , -1, -1, -1, -1 }, // 0x663,   1011 1101   O   {  1, 0,  0, 9,  9, 5,  5, 6,  6, 10, 10, 1
    {  0,  3,  8,  5,  6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0x569,   1011 1110   X
    { 10,  5,  6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0x460,   1011 1111   O   { 10, 5,  5, 6,  6, 10

    { 11,  5, 10,  7,  5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0xca0,   1100 0000   O   {  7, 5,  5, 10, 10, 11, 11, 7
    { 11,  5, 10, 11,  7,  5,  8,  3,  0, -1, -1, -1, -1, -1, -1, -1 }, // 0xda9,   1100 0001   X
    {  5, 11,  7,  5, 10, 11,  1,  9,  0, -1, -1, -1, -1, -1, -1, -1 }, // 0xea3,   1100 0010   X
    { 10,  7,  5, 10, 11,  7,  9,  8,  1,  8,  3,  1, -1, -1, -1, -1 }, // 0xfaa,   1100 0011   X
    { 11,  1,  2, 11,  7,  1,  7,  5,  1, -1, -1, -1, -1, -1, -1, -1 }, // 0x8a6,   1100 0100   O   {  1, 2,  2, 11, 11, 7,  7, 5,  5, 1
    {  0,  8,  3,  1,  2,  7,  1,  7,  5,  7,  2, 11, -1, -1, -1, -1 }, // 0x9af,   1100 0101   X
    {  9,  7,  5,  9,  2,  7,  9,  0,  2,  2, 11,  7, -1, -1, -1, -1 }, // 0xaa5,   1100 0110   O   {  0, 2,  2, 11, 11, 7,  7, 5,  5, 9,  9, 0
    {  7,  5,  2,  7,  2, 11,  5,  9,  2,  3,  2,  8,  9,  8,  2, -1 }, // 0xbac,   1100 0111   X
    {  2,  5, 10,  2,  3,  5,  3,  7,  5, -1, -1, -1, -1, -1, -1, -1 }, // 0x4ac,   1100 1000   O   {  2, 3,  3, 7,  7, 5,  5, 10, 10, 2
    {  8,  2,  0,  8,  5,  2,  8,  7,  5, 10,  2,  5, -1, -1, -1, -1 }, // 0x5a5,   1100 1001   O   {  2, 0,  0, 8,  8, 7,  7, 5,  5, 10, 10, 2
    {  9,  0,  1,  5, 10,  3,  5,  3,  7,  3, 10,  2, -1, -1, -1, -1 }, // 0x6af,   1100 1010   X
    {  9,  8,  2,  9,  2,  1,  8,  7,  2, 10,  2,  5,  7,  5,  2, -1 }, // 0x7a6,   1100 1011   X
    {  1,  3,  5,  3,  7,  5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0x0aa,   1100 1100   O   {  1, 3,  3, 7,  7, 5,  5, 1
    {  0,  8,  7,  0,  7,  1,  1,  7,  5, -1, -1, -1, -1, -1, -1, -1 }, // 0x1a3,   1100 1101   O   {  1, 0,  0, 8,  8, 7,  7, 5,  5, 1
    {  9,  0,  3,  9,  3,  5,  5,  3,  7, -1, -1, -1, -1, -1, -1, -1 }, // 0x2a9,   1100 1110   O   {  0, 3,  3, 7,  7, 5,  5, 9,  9, 0
    {  9,  8,  7,  5,  9,  7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0x3a0,   1100 1111   O   {  7, 5,  5, 9,  9, 8,  8, 7

    {  5,  8,  4,  5, 10,  8, 10, 11,  8, -1, -1, -1, -1, -1, -1, -1 }, // 0xd30,   1101 0000   O   {  4, 5,  5, 10, 10, 11, 11, 8,  8, 4
    {  5,  0,  4,  5, 11,  0,  5, 10, 11, 11,  3,  0, -1, -1, -1, -1 }, // 0xc39,   1101 0001   O   {  3, 0,  0, 4,  4, 5,  5, 10, 10, 11, 11, 3,
    {  0,  1,  9,  8,  4, 10,  8, 10, 11, 10,  4,  5, -1, -1, -1, -1 }, // 0xf33,   1101 0010   X
    { 10, 11,  4, 10,  4,  5, 11,  3,  4,  9,  4,  1,  3,  1,  4, -1 }, // 0xe3a,   1101 0011   X
    {  2,  5,  1,  2,  8,  5,  2, 11,  8,  4,  5,  8, -1, -1, -1, -1 }, // 0x936,   1101 0100   O   {  1, 2,  2, 11, 11, 8,  8, 4,  4, 5,  5, 1
    {  0,  4, 11,  0, 11,  3,  4,  5, 11,  2, 11,  1,  5,  1, 11, -1 }, // 0x83f,   1101 0101   X
    {  0,  2,  5,  0,  5,  9,  2, 11,  5,  4,  5,  8, 11,  8,  5, -1 }, // 0xb35,   1101 0110   X
    {  9,  4,  5,  2, 11,  3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0xa3c,   1101 0111   X
    {  2,  5, 10,  3,  5,  2,  3,  4,  5,  3,  8,  4, -1, -1, -1, -1 }, // 0x53c,   1101 1000   O   {  2, 3,  3, 8,  8, 4,  4, 5,  5, 10, 10, 2
    {  5, 10,  2,  5,  2,  4,  4,  2,  0, -1, -1, -1, -1, -1, -1, -1 }, // 0x435,   1101 1001   O   {  2, 0,  0, 4,  4, 5,  5, 10, 10, 2
    {  3, 10,  2,  3,  5, 10,  3,  8,  5 , 4,  5,  8,  0,  1,  9, -1 }, // 0x73f,   1101 1010   X
    {  5, 10,  2,  5,  2 , 4,  1,  9,  2,  9,  4,  2, -1, -1, -1, -1 }, // 0x636,   1101 1011   O   {  1, 9,  9, 4,  4, 5,  5, 10, 10, 2,  2, 1
    {  8,  4,  5,  8,  5,  3,  3,  5,  1, -1, -1, -1, -1, -1, -1, -1 }, // 0x13a,   1101 1100   O   {  1, 3,  3, 8,  8, 4,  4, 5,  5, 1
    {  0,  4,  5,  1,  0,  5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0x033,   1101 1101   O   {  1, 0,  0, 4,  4, 5,  5, 1
    {  8,  4,  5,  8,  5,  3,  9,  0,  5,  0,  3 , 5, -1, -1, -1, -1 }, // 0x339,   1101 1110   O   {  0, 3,  3, 8,  8, 4,  4, 5,  5, 9,  9, 0
    {  9,  4,  5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0x230,   1101 1111   O   {  9, 4,  4, 5,  5, 9

    {  4, 11,  7,  4,  9, 11,  9, 10, 11, -1, -1, -1, -1, -1, -1, -1 }, // 0xe90,   1110 0000   O   {  7, 4,  4, 9,  9, 10, 10, 11, 11, 7
    {  0,  8,  3,  4,  9,  7,  9, 11,  7,  9, 10, 11, -1, -1, -1, -1 }, // 0xf99,   1110 0001   X
    {  1, 10, 11,  1, 11,  4,  1,  4,  0,  7,  4, 11, -1, -1, -1, -1 }, // 0xc93,   1110 0010   O   {  0, 1,  1, 10, 10, 11, 11, 7,  7, 4,  4, 0
    {  3,  1,  4,  3,  4,  8,  1, 10,  4,  7,  4, 11, 10, 11,  4, -1 }, // 0xd9a,   1110 0011   X
    {  4, 11,  7,  9, 11,  4,  9,  2, 11,  9,  1,  2, -1, -1, -1, -1 }, // 0xa96,   1110 0100   O   {  1, 2,  2, 11, 11, 7,  7, 4,  4, 9,  9, 1
    {  9,  7,  4,  9, 11,  7,  9,  1, 11,  2, 11,  1,  0,  8,  3, -1 }, // 0xb9f,   1110 0101   X
    { 11,  7,  4, 11,  4,  2,  2,  4,  0, -1, -1, -1, -1, -1, -1, -1 }, // 0x895,   1110 0110   O   {  0, 2,  2, 11, 11, 7,  7, 4,  4, 0
    { 11,  7,  4, 11,  4,  2,  8,  3,  4,  3,  2,  4, -1, -1, -1, -1 }, // 0x99c,   1110 0111   O   {  3, 2,  2, 11, 11, 7,  7, 4,  4, 8,  8, 3
    {  2,  9, 10,  2,  7,  9,  2,  3,  7,  7,  4,  9, -1, -1, -1, -1 }, // 0x69c,   1110 1000   O   {  2, 3,  3, 7,  7, 4,  4, 9,  9, 10, 10, 2
    {  9, 10,  7,  9,  7,  4, 10,  2,  7,  8,  7,  0,  2,  0,  7, -1 }, // 0x795,   1110 1001   X
    {  3,  7, 10,  3, 10,  2,  7,  4, 10,  1, 10,  0,  4,  0, 10, -1 }, // 0x49f,   1110 1010   X
    {  1, 10,  2,  8,  7,  4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0x596,   1110 1011   X
    {  4,  9,  1 , 4,  1,  7,  7,  1,  3, -1, -1, -1, -1, -1, -1, -1 }, // 0x29a,   1110 1100   O   {  1, 3,  3, 7,  7, 4,  4, 9,  9, 1
    {  4,  9,  1,  4,  1,  7,  0,  8,  1,  8,  7,  1, -1, -1, -1, -1 }, // 0x393,   1110 1101   O   {  1, 0,  0, 8,  8, 7,  7, 4,  4, 9,  9, 1
    {  4,  0,  3,  7,  4,  3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0x099,   1110 1110   O   {  0, 3,  3, 7,  7, 4,  4, 0
    {  4,  8,  7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0x190,   1110 1111   O   {  4, 8,  8, 7,  7, 4

    {  9, 10,  8, 10, 11,  8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0xf00,   1111 0000   O   {  8, 9,  9, 10, 10, 11, 11, 8
    {  3,  0,  9,  3,  9, 11, 11,  9, 10, -1, -1, -1, -1, -1, -1, -1 }, // 0xe09,   1111 0001   O   {  3, 0,  0, 9,  9, 10, 10, 11, 11, 3
    {  0,  1, 10,  0, 10,  8,  8, 10, 11, -1, -1, -1, -1, -1, -1, -1 }, // 0xd03,   1111 0010   O   {  0, 1,  1, 10, 10, 11, 11, 8,  8, 0
    {  3,  1, 10, 11,  3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0xc0a,   1111 0011   O   {  3, 1,  1, 10, 10, 11, 11, 3
    {  1,  2, 11,  1, 11,  9,  9, 11,  8, -1, -1, -1, -1, -1, -1, -1 }, // 0xb06,   1111 0100   O   {  1, 2,  2, 11, 11, 8,  8, 9,  9, 1
    {  3,  0,  9,  3,  9, 11,  1,  2,  9,  2, 11,  9, -1, -1, -1, -1 }, // 0xa0f,   1111 0101   O   {  0, 9,  9, 1,  1, 2,  2, 11, 11, 3,  3, 0
    {  0,  2, 11,  8,  0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0x905,   1111 0110   O   {  0, 2,  2, 11, 11, 8,  8, 0
    {  3,  2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0x80c,   1111 0111   O   {  3, 2,  2, 11, 11, 3
    {  2,  3,  8,  2,  8, 10, 10,  8,  9, -1, -1, -1, -1, -1, -1, -1 }, // 0x70c,   1111 1000   O   {  2, 3,  3, 8,  8, 9,  9, 10, 10, 2
    {  9, 10,  2,  0,  9,  2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0x605,   1111 1001   O   {  2, 0,  0, 9,  9, 10, 10, 2
    {  2,  3,  8,  2,  8, 10,  0,  1,  8,  1, 10,  8, -1, -1, -1, -1 }, // 0x50f,   1111 1010   O   {  0, 1,  1, 10, 10, 2,  2, 3,  3, 8,  8, 0
    {  1, 10,  2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0x406,   1111 1011   O   {  1, 10, 10, 2,  2, 1
    {  1,  3,  8,  9,  1,  8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0x30a,   1111 1100   O   {  1, 3,  3, 8,  8, 9,  9, 1
    {  0,  9,  1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0x203,   1111 1101   O   {  0, 9,  9, 1,  1, 0
    {  0,  3,  8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }, // 0x109,   1111 1110   O   {  0, 3,  3, 8,  8, 0
    { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }  // 0x000,   1111 1111   X
};

}  // namespace jet

#endif  // SRC_JET_MARCHING_CUBES_TABLE_H_
