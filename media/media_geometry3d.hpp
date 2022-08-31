#ifndef _MEDIA_GEOMETRY3D_H_
#define _MEDIA_GEOMETRY3D_H_

#include <iostream>
#include <vector>
#include <cfloat>
#include <cassert>

#define EPS 1e-10

struct Point3{
	float x, y, z;
	Point3(float x = 0, float y = 0, float z = 0):x(x), y(y), z(z) {}
};

typedef Point3 Vector3;

Vector3 operator + (const Vector3 &A, const Vector3 &B);
Vector3 operator - (const Vector3 &A, const Vector3 &B);
Vector3 operator * (const Vector3 &A, float p);
Vector3 operator / (const Vector3 &A, float p); 

float Dot(const Vector3 &A, const Vector3 &B);
float Length(const Vector3 &A);
float Angle(const Vector3 &A, const Vector3 &B);
float DistanceP2P(const Point3 &A, const Point3 &B);
float DistanceToPlane(const Point3 &p, const Point3 &p0, const Vector3 &n);
Point3 GetplaneProjection(const Point3 &p, const Point3 &p0, const Vector3 &n);
Point3 LinePlaneIntersection(
    const Point3 &p1, const Point3 &p2, const Point3 &p0, 
    const Vector3 &n);
Vector3 Cross(const Vector3 &A, const Vector3 &B);
float Area2 (const Point3 &A, const Point3 &B, const Point3 &C);
bool isPointInTri(const Point3 &P, const Point3 &P0, const Point3 &P1, const Point3 &P2);
//float Volume6(Point3 A, Point3 B, Point3 C, Point3 D);

/* 
 * If the number of points is 3, calculate the normal directly;
 * if the number of points > 3, just use first two points.
 */
struct Face {
    std::vector<Point3> v;
    Vector3 normal() const { 
        assert(v.size() > 2);
        Vector3 n = Cross(v[1]-v[0], v[2]-v[0]);
        /* Normalized the normal vector */
        float d = Length(n);
        n = n/d;
        return n;
    }
    Face() = default;
     ~Face() = default;
};

struct Mesh3 {
    Point3 v[8]; // vetex: counterclockwise, and from z = 0 plane to z = +z plane
    Mesh3 (Point3 A, Point3 B, Point3 C, Point3 D, Point3 E, Point3 F, Point3 G, Point3 H) {
        this->v[0] = A;
        this->v[1] = B;
        this->v[2] = C;
        this->v[3] = D;
        this->v[4] = E;
        this->v[5] = F;
        this->v[6] = G;
        this->v[7] = H;
    }
};

bool isPointInPolyhedron(const Point3 &p, const std::vector<Face> &fs);

#ifdef __cplusplus
extern "C" {
#endif
bool isPointInHexahedron(float px, float py, float pz,
                         float *vx, float *vy, float *vz);
bool isPointInHexahedron_strict(float px, float py, float pz,
                               float *vx, float *vy, float *vz);
#ifdef __cplusplus
}
#endif /* extern C */


#endif /*media_geometry*/
