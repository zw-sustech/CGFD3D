#include <iostream>
#include <vector>
#include <cfloat>
#include <cmath>
#include "media_geometry3d.hpp"
//using namespace std;

int dcmp(float x) {
	if (fabs(x) < EPS) 
		return 0;
	else
		return x < 0?-1:1;
}

Vector3 operator + (const Vector3 &A, const Vector3 &B) {
	return Vector3(A.x+B.x, A.y+B.y, A.z+B.z);
}

Vector3 operator - (const Vector3 &A, const Vector3 &B) {
	return Vector3(A.x-B.x, A.y-B.y, A.z-B.z);
}

Vector3 operator * (const Vector3 &A, float p) {
	return Vector3(A.x*p, A.y*p, A.z*p);
}

Vector3 operator / (const Vector3 &A, float p) {
	return Vector3(A.x/p, A.y/p, A.z/p);
}

float Dot(const Vector3 &A, const Vector3 &B) {
	return A.x*B.x + A.y*B.y + A.z*B.z;
}


float Length(const Vector3 &A) {
	return sqrt(Dot(A, A));
}

float DistanceP2P(const Point3 &A, const Point3 &B) {
	return Length(B-A);	
}

float Angle(const Vector3 &A, const Vector3 &B) {
	return acos(Dot(A, B)/Length(A)/Length(B));
}

float DistanceToPlane(const Point3 &p, const Point3 &p0, const Vector3 &n) {
	return fabs(Dot(p-p0, n));
}


Point3 GetPlaneProjection(const Point3 &p, const Point3 &p0, const Vector3 &n) {
	return p - n*Dot(p-p0, n);
}

Point3 LinePlaneIntersection(
	const Point3 &p1, const Point3 &p2, const Point3 &p0, 
    const Vector3 &n)
{
	Vector3 v = p2 - p1;
	float t = (Dot(n, p0-p1) / Dot(n, p2-p1)); // denominator = 0?
	return p1 + v*t;
}

Vector3 Cross(const Vector3 &A, const Vector3 &B) {
	return Vector3(A.y*B.z - A.z*B.y, A.z*B.x - A.x*B.z, A.x*B.y - A.y*B.x);
}

float Area2(const Point3 &A, const Point3 &B, const Point3 &C) {
	return Length(Cross(B-A, C-A));
}


bool isPointInTri(const Point3 &P, const Point3 &P0, const Point3 &P1, const Point3 &P2) {
	float area1 = Area2(P, P0, P1);
	float area2 = Area2(P, P1, P2);
	float area3 = Area2(P, P2, P0);
	return dcmp(area1+area2+area3 - Area2(P0, P1, P2)) == 0;
}
/*
float Volume6(Point3 A, Point3 B, Point3 C, Point3 D) {
	return Dot(D-A, Cross(B-A, C-A));
}
*/

/* 
 * Check whether the point is in the polyhedron.
 * Note: The hexahedron must be convex!
 */
bool isPointInPolyhedron(const Point3 &p, const std::vector<Face> &fs) {
	for (Face const &f:fs) {
		// if the point on the triangle, return true 
		if ( isPointInTri(p, f.v[0], f.v[1], f.v[2]) ) {
			return true;
		}
		Vector3 p2f = f.v[0]-p;
		Vector3 A = f.normal();
		float sign = Dot(p2f, f.normal());
		sign /= Length(p2f); // normalization, for stability

		constexpr float bound = -1e-15;
		if(sign < 0.0) return false;
	}
	return true;
}

/*
 * Input: vx, vy, vz are the EIGHT vertices of the hexahedron 
 *
 *    ↑ +z       4----6
 *    |         /|   /|
 *             / 0--/-2
 *            5----7 /
 *            |/   |/
 *            1----3
 *
 */
bool isPointInHexahedron(float px, float py, float pz,
                         float *vx, float *vy, float *vz) 
{

	Point3 P(px, py, pz); 

	/* 
	 * Just for cgfd3D, in which the grid mesh maybe not a hexahedron,
	 */
	std::vector<Face> hexa{
		Face{ {Point3{vx[0], vy[0], vz[0]}, Point3{vx[4], vy[4], vz[4]}, Point3{vx[6], vy[6], vz[6]}} }, 
//		Face{ {Point3{vx[6], vy[6], vz[6]}, Point3{vx[2], vy[2], vz[2]}, Point3{vx[0], vy[0], vz[0]}} }, // back

		Face{ {Point3{vx[7], vy[7], vz[7]}, Point3{vx[5], vy[5], vz[5]}, Point3{vx[1], vy[1], vz[1]}} },
//		Face{ {Point3{vx[1], vy[1], vz[1]}, Point3{vx[3], vy[3], vz[3]}, Point3{vx[7], vy[7], vz[7]}} }, // front

		Face{ {Point3{vx[5], vy[5], vz[5]}, Point3{vx[4], vy[4], vz[4]}, Point3{vx[0], vy[0], vz[0]}} }, 
//		Face{ {Point3{vx[0], vy[0], vz[0]}, Point3{vx[1], vy[1], vz[1]}, Point3{vx[5], vy[5], vz[5]}} }, // left

		Face{ {Point3{vx[2], vy[2], vz[2]}, Point3{vx[6], vy[6], vz[6]}, Point3{vx[7], vy[7], vz[7]}} }, 
//		Face{ {Point3{vx[7], vy[7], vz[7]}, Point3{vx[3], vy[3], vz[3]}, Point3{vx[2], vy[2], vz[2]}} }, // right

		Face{ {Point3{vx[4], vy[4], vz[4]}, Point3{vx[5], vy[5], vz[5]}, Point3{vx[7], vy[7], vz[7]}} }, 
//		Face{ {Point3{vx[7], vy[7], vz[7]}, Point3{vx[6], vy[6], vz[6]}, Point3{vx[4], vy[4], vz[4]}} }, // top

//		Face{ {Point3{vx[0], vy[0], vz[0]}, Point3{vx[2], vy[2], vz[2]}, Point3{vx[3], vy[3], vz[3]}} }, 
		Face{ {Point3{vx[3], vy[3], vz[3]}, Point3{vx[1], vy[1], vz[1]}, Point3{vx[0], vy[0], vz[0]}} } // bottom
	};

  return isPointInPolyhedron(P, hexa);
}


/*
 * Input: vx, vy, vz are the EIGHT vertices of the hexahedron 
 *
 *    ↑ +z       4----6
 *    |         /|   /|
 *             / 0--/-2
 *            5----7 /
 *            |/   |/
 *            1----3
 *
 *
 * A more strict function, 
 *  if isPointInHexahedron cannot find the cell, use it. 
 */
bool isPointInHexahedron_strict(float px, float py, float pz,
                               float *vx, float *vy, float *vz) 
{

	Point3 P(px, py, pz); 

	Point3 v0(vx[0], vy[0], vz[0]);
	Point3 v1(vx[1], vy[1], vz[1]);
	Point3 v2(vx[2], vy[2], vz[2]);
	Point3 v3(vx[3], vy[3], vz[3]);
	Point3 v4(vx[4], vy[4], vz[4]);
	Point3 v5(vx[5], vy[5], vz[5]);
	Point3 v6(vx[6], vy[6], vz[6]);
	Point3 v7(vx[7], vy[7], vz[7]);

	/* 
	 * One hexahedron cut into five tetrahedron.
	 * Since the four points may not coplanar, cut twice, ten tetrahedron.
	 */
	std::vector< std::vector<Face> > tetrahs(10, std::vector<Face>(4));

	// 0, 1, 2, 4
	tetrahs[0] = { 
		Face{{v0, v4, v2}}, Face{{v2, v4, v1}}, Face{{v1, v4, v0}}, Face{{v0, v2, v1}}
	};
	// 2, 6, 4, 7 
	tetrahs[1] = { 
		Face{{v4, v6, v2}}, Face{{v2, v6, v7}}, Face{{v7, v6, v4}}, Face{{v4, v2, v7}}
	};
	// 3, 1, 7, 2
	tetrahs[2] = {
		Face{{v2, v7, v3}}, Face{{v3, v7, v1}}, Face{{v1, v2, v3}}, Face{{v1, v7, v2}}
	};
	// 1, 7, 4, 5
	tetrahs[3] = {
		Face{{v7, v4, v5}}, Face{{v1, v7, v5}}, Face{{v5, v4, v1}}, Face{{v1, v4, v7}}
	};
	// 1, 4, 2, 7
	tetrahs[4] = {
		Face{{v2, v4, v7}}, Face{{v7, v4, v1}}, Face{{v2, v7, v1}}, Face{{v4, v2, v1}}
	};


	// 4, 5, 6, 0
	tetrahs[5] = {
		Face{{v4, v5, v6}}, Face{{v4, v6, v0}}, Face{{v4, v0, v5}}, Face{{v6, v5, v0}}
	};
	// 0, 2, 6, 3
	tetrahs[6] = {
		Face{{v6, v2, v0}}, Face{{v6, v3, v2}}, Face{{v0, v2, v3}}, Face{{v3, v6, v0}}
	};
	// 5, 6, 3, 7
	tetrahs[7] = {
		Face{{v6, v5, v7}}, Face{{v3, v6, v7}}, Face{{v7, v5, v3}}, Face{{v5, v6, v3}}
	};
	// 5, 0, 1, 3
	tetrahs[8] = {
		Face{{v5, v0, v1}}, Face{{v3, v5, v1}}, Face{{v0, v3, v1}}, Face{{v0, v5, v3}}
	};
	// 0, 5, 6, 3
	tetrahs[9] = {
		Face{{v5, v6, v0}}, Face{{v5, v0, v3}}, Face{{v0, v6, v3}}, Face{{v3, v6, v5}}
	};

	for (auto &tetrah:tetrahs) {
		if (isPointInPolyhedron(P, tetrah)) {
				return true;
		}
	}

  return false;
}