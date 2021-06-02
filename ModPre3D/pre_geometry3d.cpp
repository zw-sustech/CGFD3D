#include <iostream>
#include <vector>
#include <cfloat>
#include <cmath>
#include "pre_geometry3d.hpp"
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

Point3 GetplaneProjection(const Point3 &p, const Point3 &p0, const Vector3 &n) {
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





