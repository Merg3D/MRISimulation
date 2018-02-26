#pragma once

#include <string>
#include <iostream>

template<class T> class Vec4
{
public:
	Vec4() { x = y = z = w = 0.0; }
	Vec4(T a) { x = y = z = w = a; }
	Vec4(T p_x, T p_y, T p_z, T p_w) { x = p_x; y = p_y; z = p_z; w = p_w; }
	
	Vec4(Vec3<T> xyz, T p_w) { x = xyz.x; y = xyz.y; z = xyz.z; w = p_w; }
	Vec4(T p_x, Vec3<T> yzw) { x = p_x; y = yzw.x; z = yzw.y; w = yzw.z; }

	//General operators
	bool operator==(const Vec4 &r) const { return x == r.x && y == r.y && z == r.z && w == r.w; }
	bool operator!=(const Vec4 &r) const { return x != r.x || y != r.y || z != r.z || w != r.w; }

	//Vec4 operators
	Vec4 operator+(const Vec4 &r) const { return Vec4(x + r.x, y + r.y, z + r.z, w + r.w); }
	Vec4 operator-(const Vec4 &r) const { return Vec4(x - r.x, y - r.y, z - r.z, w - r.w); }
	Vec4 operator*(const Vec4 &r) const { return Vec4(x * r.x, y * r.y, z * r.z, w * r.w); }
	Vec4 operator/(const Vec4 &r) const { return Vec4(x / r.x, y / r.y, z / r.z, w / r.w); }

	Vec4& operator+=(const Vec4 &r) { x += r.x; y += r.y; z += r.z; w += r.w; return *this; }
	Vec4& operator-=(const Vec4 &r) { x -= r.x; y -= r.y; z -= r.z; w -= r.w; return *this; }
	Vec4& operator*=(const Vec4 &r) { x *= r.x; y *= r.y; z *= r.z; w *= r.w; return *this; }
	Vec4& operator/=(const Vec4 &r) { x /= r.x; y /= r.y; z /= r.z; w /= r.w; return *this; }

	//Scalar operators
	Vec4 operator+(const T &r) const { return Vec4(x + r, y + r, z + r, w + r); }
	Vec4 operator-(const T &r) const { return Vec4(x - r, y - r, z - r, w - r); }
	Vec4 operator*(const T &r) const { return Vec4(x * r, y * r, z * r, w * r); }
	Vec4 operator/(const T &r) const { return Vec4(x / r, y / r, z / r, w / r); }

	Vec4& operator+=(const T &r) { x += r; y += r; z += r; w += r; return *this; }
	Vec4& operator-=(const T &r) { x -= r; y -= r; z -= r; w -= r; return *this; }
	Vec4& operator*=(const T &r) { x *= r; y *= r; z *= r; w *= r; return *this; }
	Vec4& operator/=(const T &r) { x /= r; y /= r; z /= r; w /= r; return *this; }

	//Index operator
	T& operator[](const int &i) { return v[i]; }

	//Conversion
	std::string to_string() const
	{
		std::string s = "{ " + std::to_string(x) + ", " + std::to_string(y) + ", " + std::to_string(z) + ", " + std::to_string(w) + " }";
		return s;
	}

	//Get info
	T length() const { return sqrt(x * x + y * y + z * z + w * w); }
	T angle(const Vec4<T> &v) const
	{
		T angle = fabs(acos(dot(v) / (length() * v.length())));
		return angle;
	}
	T dot(const Vec4<T> &v) const
	{
		T result = x * v.x + y * v.y + z * v.z + w * v.w;
		return result;
	}

	//properties
	T& get_x() { return x; }
	T& get_y() { return y; }
	T& get_z() { return z; }
	T& get_w() { return w; }

	void set_x(const T &p_x) { x = p_x; }
	void set_y(const T &p_y) { y = p_y; }
	void set_z(const T &p_z) { z = p_z; }
	void set_w(const T &p_w) { w = p_w; }
	
	Vec3<T> get_xyz() const
	{
		return{ x, y, z };
	}

	//Set
	Vec4& normalize()
	{
		T inv = 1.0f / length();
		x *= inv;
		y *= inv;
		z *= inv;
		w *= inv;
		return *this;
	}

	//Components
	union
	{
		struct { T x, y, z, w; };
		struct { T r, g, b, a; };
		T v[4];
	};
};

typedef Vec4<double> vec4d;
typedef Vec4<float> vec4;
typedef Vec4<int> vec4i;