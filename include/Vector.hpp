#ifndef FDTD_VECTOR_H_
#define FDTD_VECTOR_H_

#include <cassert>
#include <array>
#include <iostream>

template <typename T>
class Vector3 {
    private:
    T data[3];

    public:
    Vector3() : data{0, 0, 0} { };
    Vector3(T d[3]) : data{d[0], d[1], d[2]} { };
    Vector3(T d1, T d2, T d3) : data{d1, d2, d3} { }
    Vector3(const std::array<T, 3>& d) : data{d[0], d[1], d[2]} { }
    Vector3(const std::array<T, 3>&& d) : data{d[0], d[1], d[2]} { }
    Vector3(const Vector3& v) : data{v[0], v[1], v[2]} { }
    Vector3(const Vector3&& v) : data{v[0], v[1], v[2]} { }

    const T* GetData() const {
        return data;
    }

    T GetData(int i) const {
        assert(i>=0 && i<3);
        return data[i];
    }

    //---------------------- operators
    Vector3& operator=(const Vector3& v) {
        data = v.GetData();
        return *this;
    }

    friend Vector3 operator+(const Vector3& v1, const Vector3& v2) {
        const T* v1_data = v1.GetData();
        const T* v2_data = v2.GetData();
        T v_data[3]{v1_data[0] + v2_data[0],
                    v1_data[1] + v2_data[1],
                    v1_data[2] + v2_data[2]};
        return Vector3(v_data);
    }

    friend Vector3 operator+(const Vector3& v1, const T d) {
        const T* v1_data = v1.GetData();
        T v_data[3]{v1_data[0] + d,
                    v1_data[1] + d,
                    v1_data[2] + d};
        return Vector3(v_data);
    }

    friend Vector3 operator+(const T d, const Vector3& v1) {
        const T* v1_data = v1.GetData();
        T v_data[3]{d + v1_data[0],
                    d + v1_data[1],
                    d + v1_data[2]};
        return Vector3(v_data);
    }

    Vector3& operator+=(const Vector3& v) {
        const T* v_data = v.GetData();
        data[0] += v_data[0];
        data[1] += v_data[1];
        data[2] += v_data[2];
        return *this;
    }

    Vector3& operator+=(const T d) {
        data[0] += d;
        data[1] += d;
        data[2] += d;
        return *this;
    }

    friend Vector3 operator-(const Vector3& v1, const Vector3& v2) {
        const T* v1_data = v1.GetData();
        const T* v2_data = v2.GetData();
        T v_data[3]{v1_data[0] - v2_data[0],
                    v1_data[1] - v2_data[1],
                    v1_data[2] - v2_data[2]};
        return Vector3(v_data);
    }

    friend Vector3 operator-(const Vector3& v1, const T d) {
        const T* v1_data = v1.GetData();
        T v_data[3]{v1_data[0] - d,
                    v1_data[1] - d,
                    v1_data[2] - d};
        return Vector3(v_data);
    }

    friend Vector3 operator-(const T d, const Vector3& v1) {
        const T* v1_data = v1.GetData();
        T v_data[3]{d - v1_data[0],
                    d - v1_data[1],
                    d - v1_data[2]};
        return Vector3(v_data);
    }

    Vector3& operator-=(const Vector3& v) {
        const T* v_data = v.GetData();
        data[0] -= v_data[0];
        data[1] -= v_data[1];
        data[2] -= v_data[2];
        return *this;
    }

    Vector3& operator-=(const T d) {
        data[0] -= d;
        data[1] -= d;
        data[2] -= d;
        return *this;
    }

    friend Vector3 operator*(const Vector3& v1, const Vector3& v2) {
        const T* v1_data = v1.GetData();
        const T* v2_data = v2.GetData();
        T v_data[3]{v1_data[0] * v2_data[0],
                    v1_data[1] * v2_data[1],
                    v1_data[2] * v2_data[2]};
        return Vector3(v_data);
    }

    friend Vector3 operator*(const Vector3& v1, const T d) {
        const T* v1_data = v1.GetData();
        T v_data[3]{v1_data[0] * d,
                    v1_data[1] * d,
                    v1_data[2] * d};
        return Vector3(v_data);
    }

    friend Vector3 operator*(const T d, const Vector3& v1) {
        const T* v1_data = v1.GetData();
        T v_data[3]{d * v1_data[0],
                    d * v1_data[1],
                    d * v1_data[2]};
        return Vector3(v_data);
    }

    Vector3& operator*=(const Vector3& v) {
        const T* v_data = v.GetData();
        data[0] *= v_data[0];
        data[1] *= v_data[1];
        data[2] *= v_data[2];
        return *this;
    }

    Vector3& operator*=(const T d) {
        data[0] *= d;
        data[1] *= d;
        data[2] *= d;
        return *this;
    }

    friend Vector3 operator/(const Vector3& v1, const Vector3& v2) {
        const T* v1_data = v1.GetData();
        const T* v2_data = v2.GetData();
        T v_data[3]{v1_data[0] / v2_data[0],
                    v1_data[1] / v2_data[1],
                    v1_data[2] / v2_data[2]};
        return Vector3(v_data);
    }

    friend Vector3 operator/(const Vector3& v1, const T d) {
        const T* v1_data = v1.GetData();
        T v_data[3]{v1_data[0] / d,
                    v1_data[1] / d,
                    v1_data[2] / d};
        return Vector3(v_data);
    }

    friend Vector3 operator/(const T d, const Vector3& v1) {
        const T* v1_data = v1.GetData();
        T v_data[3]{d / v1_data[0],
                    d / v1_data[1],
                    d / v1_data[2]};
        return Vector3(v_data);
    }

    Vector3& operator/=(const Vector3& v) {
        const T* v_data = v.GetData();
        data[0] /= v_data[0];
        data[1] /= v_data[1];
        data[2] /= v_data[2];
        return *this;
    }

    Vector3& operator/=(const T d) {
        data[0] /= d;
        data[1] /= d;
        data[2] /= d;
        return *this;
    }

    T operator[](std::size_t i) const {
        assert(i>=0 && i<3);
        return data[i];
    }

    friend std::ostream& operator<<(std::ostream& out, const Vector3& v)
    {
        return out << "[" << v[0] << ", " << v[1] << ", " << v[2] << "]" ;
    }

};


#endif  // FDTD_VECTOR_H_
