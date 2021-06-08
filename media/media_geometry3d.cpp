#include <iostream>
#include <vector>
#include <cfloat>
#include <cmath>
#include "media_geometry3d.hpp"
//using namespace std;

int dcmp(double x) {
	if (fabs(x) < EPS) 
		return 0;
	else
		return x < 0?-1:1;
}

Vector3 operator + (Vector3 A, Vector3 B) {
	return Vector3(A.x+B.x, A.y+B.y, A.z+B.z);
}

Vector3 operator - (Vector3 A, Vector3 B) {
	return Vector3(A.x-B.x, A.y-B.y, A.z-B.z);
}

Vector3 operator * (Vector3 A, double p) {
	return Vector3(A.x*p, A.y*p, A.z*p);
}

Vector3 operator / (Vector3 A, double p) {
	return Vector3(A.x/p, A.y/p, A.z/p);
}

double Dot(Vector3 A, Vector3 B) {
	return A.x*B.x + A.y*B.y + A.z*B.z;
}


double Length(Vector3 A) {
	return sqrt(Dot(A, A));
}

double DistanceP2P(Point3 A, Point3 B) {
	return Length(B-A);	
}


double Angle(Vector3 A, Vector3 B) {
	return acos(Dot(A, B)/Length(A)/Length(B));
}

double DistanceToPlane(const Point3 &p, const Point3 &p0, const Vector3 &n) {
	return fabs(Dot(p-p0, n));
}


Point3 GetPlaneProjection(const Point3 &p, const Point3 &p0, const Vector3 &n) {
	return p - n*Dot(p-p0, n);
}

Point3 LinePlaneIntersection(Point3 p1, Point3 p2, Point3 p0, Vector3 n){
	Vector3 v = p2 - p1;
	double t = (Dot(n, p0-p1) / Dot(n, p2-p1)); // denominator = 0?
	return p1 + v*t;
}

Vector3 Cross(Vector3 A, Vector3 B) {
	return Vector3(A.y*B.z - A.z*B.y, A.z*B.x - A.x*B.z, A.x*B.y - A.y*B.x);
}

double Area2 (Point3 A, Point3 B, Point3 C) {
	return Length(Cross(B-A, C-A));
}

bool PointInTri(Point3 P, Point3 P0, Point3 P1, Point3 P2) {
	double area1 = Area2(P, P0, P1);
	double area2 = Area2(P, P1, P2);
	double area3 = Area2(P, P2, P0);
	return dcmp(area1+area2+area3 - Area2(P0, P1, P2)) == 0;
}

double Volume6(Point3 A, Point3 B, Point3 C, Point3 D) {
	return Dot(D-A, Cross(B-A, C-A));
}


// check whether the point is in the polyhedron
bool isPointInPolyhedron(const Point3 &p, const std::vector<Face> &fs) {
	for (Face const &f:fs) {
		Vector3 p2f = f.v[0]-p;
//		Vector3 A = f.normal();
		double sign = Dot(p2f, f.normal());

		constexpr double bound = -1e-15;
		if(sign < bound) return false;
	}
	return true;
}

/*
 * Input: vx, vy, vz are the EIGHT vertexes of the hexahedron 
 *
 *    ↑ +z       4----6
 *    |         /|   /|
 *             / 0--/-2
 *            5----7 /
 *            |/   |/
 *            1----3
 *
 * Note: The hexahedron must be convex!
 *
 */
bool isPointInHexahedron(float px, float py, float pz,
                         float *vx, float *vy, float *vz) 
{
	Point3 P(px, py, pz); 
	std::vector<Face> hexa{
		Face{ {Point3{vx[0], vy[0], vz[0]}, Point3{vx[4], vy[4], vz[4]}, Point3{vx[6], vy[6], vz[6]}} }, // back
		Face{ {Point3{vx[1], vy[1], vz[1]}, Point3{vx[3], vy[3], vz[3]}, Point3{vx[7], vy[7], vz[7]}} }, // front
		Face{ {Point3{vx[5], vy[5], vz[5]}, Point3{vx[4], vy[4], vz[4]}, Point3{vx[0], vy[0], vz[0]}} }, // left
		Face{ {Point3{vx[6], vy[6], vz[6]}, Point3{vx[7], vy[7], vz[7]}, Point3{vx[3], vy[3], vz[3]}} }, // right
		Face{ {Point3{vx[4], vy[4], vz[4]}, Point3{vx[5], vy[5], vz[5]}, Point3{vx[7], vy[7], vz[7]}} }, // top
		Face{ {Point3{vx[0], vy[0], vz[0]}, Point3{vx[2], vy[2], vz[2]}, Point3{vx[3], vy[3], vz[3]}} }, // bottom
	};



    return isPointInPolyhedron(P, hexa);
}