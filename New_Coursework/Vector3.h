#pragma once
#include <math.h>
#include <iostream>

class Vector3;

class Vector3
{
public:
	Vector3(void);
	Vector3(const float _x, const float _y, const float _z);
	~Vector3(void);

	bool operator==(const Vector3 & rhs) const;
	Vector3 & operator+=(const Vector3 & rhs);
	Vector3 & operator-=(const Vector3 & rhs);
	Vector3 & operator*=(const Vector3 & rhs);
	Vector3 & operator/=(const Vector3 & rhs);

	float Length();

	float LengthSqr();

	// Output
	friend std::ostream& operator<<(std::ostream& out, Vector3& object);

	Vector3 Normalize();
	void NormalizeSelf();
	union {
		struct { float x, y, z; };
		struct { float r, g, b; };
		struct { float u, v, w; };
	};
};

inline Vector3 operator * (const Vector3 &v,
						   const float &s)
{
	return Vector3(v.x * s, v.y * s, v.z * s);
}

inline Vector3 operator / (const Vector3 &v,
						   const float &s)
{
	return Vector3(v.x / s, v.y / s, v.z / s);
}

inline Vector3 operator + (const Vector3 &v,
						   const float &s)
{
	return Vector3(v.x + s, v.y + s, v.z + s);
}

inline Vector3 operator - (const Vector3 &v,
						   const float &s)
{
	return Vector3(v.x - s, v.y - s, v.z - s);
}

inline Vector3 operator - (const Vector3 &v,
						   const Vector3 &s)
{
	return Vector3(v.x - s.x, v.y - s.y, v.z - s.z);
}

inline Vector3 operator + (const Vector3 &v,
						   const Vector3 &s)
{
	return Vector3(v.x + s.x, v.y + s.y, v.z + s.z);
}

inline float Dot(const Vector3 &lhs, const Vector3 &rhs)
{
	return (lhs.x * rhs.x + lhs.y * rhs.y + lhs.z + rhs.z);
}

inline Vector3 Cross(const Vector3 &lhs, const Vector3 &rhs)
{
	return Vector3(lhs.y * rhs.z - lhs.z * rhs.y, lhs.z * rhs.x - lhs.x * rhs.z, lhs.x * rhs.y - lhs.y * rhs.x); 
}

inline Vector3 Normal(const Vector3& v1, const Vector3& v2, const Vector3& v3)
{
	Vector3 vec1 = v2 - v1;
	Vector3 vec2 = v3 - v1;
	return Cross(vec1, vec2).Normalize();
}
