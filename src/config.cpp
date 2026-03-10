#include "config.h"

//complex & int
complex_type operator+(const complex_type& lhs, int rhs) {
    return lhs + static_cast<float_type>(rhs);
}

complex_type operator+(int lhs, const complex_type& rhs) {
    return static_cast<float_type>(lhs) + rhs;
}

complex_type operator-(const complex_type& lhs, int rhs) {
    return lhs - static_cast<float_type>(rhs);
}

complex_type operator-(int lhs, const complex_type& rhs) {
    return static_cast<float_type>(lhs) - rhs;
}

complex_type operator*(const complex_type& lhs, int rhs) {
    return lhs * static_cast<float_type>(rhs);
}

complex_type operator*(int lhs, const complex_type& rhs) {
    return static_cast<float_type>(lhs) * rhs;
}

complex_type operator/(const complex_type& lhs, int rhs) {
    return lhs / static_cast<float_type>(rhs);
}

complex_type operator/(int lhs, const complex_type& rhs) {
    return static_cast<float_type>(lhs) / rhs;
}


//complex & double
complex_type operator+(const complex_type& lhs, double rhs) {
    return lhs + static_cast<float_type>(rhs);
}

complex_type operator+(double lhs, const complex_type& rhs) {
    return static_cast<float_type>(lhs) + rhs;
}

complex_type operator-(const complex_type& lhs, double rhs) {
    return lhs - static_cast<float_type>(rhs);
}

complex_type operator-(double lhs, const complex_type& rhs) {
    return static_cast<float_type>(lhs) - rhs;
}

complex_type operator*(const complex_type& lhs, double rhs) {
    return lhs * static_cast<float_type>(rhs);
}

complex_type operator*(double lhs, const complex_type& rhs) {
    return static_cast<float_type>(lhs) * rhs;
}

complex_type operator/(const complex_type& lhs, double rhs) {
    return lhs / static_cast<float_type>(rhs);
}

complex_type operator/(double lhs, const complex_type& rhs) {
    return static_cast<float_type>(lhs) / rhs;
}


// float & double
float_type operator*(const float_type& lhs, double rhs) {
    return lhs * static_cast<float_type>(rhs);
}

float_type operator*(double lhs, const float_type& rhs) {
    return static_cast<float_type>(lhs) * rhs;
}

float_type operator+(const float_type& lhs, double rhs) {
    return lhs + static_cast<float_type>(rhs);
}

float_type operator+(double lhs, const float_type& rhs) {
    return static_cast<float_type>(lhs) + rhs;
}

float_type operator-(const float_type& lhs, double rhs) {
    return lhs - static_cast<float_type>(rhs);
}

float_type operator-(double lhs, const float_type& rhs) {
    return static_cast<float_type>(lhs) - rhs;
}

float_type operator/(const float_type& lhs, double rhs) {
    return lhs / static_cast<float_type>(rhs);
}

float_type operator/(double lhs, const float_type& rhs) {
    return static_cast<float_type>(lhs) / rhs;
}


//float & int
float_type operator+(const float_type& lhs, int rhs) {
    return lhs + static_cast<float_type>(rhs);
}

float_type operator+(int lhs, const float_type& rhs) {
    return static_cast<float_type>(lhs) + rhs;
}

float_type operator-(const float_type& lhs, int rhs) {
    return lhs - static_cast<float_type>(rhs);
}

float_type operator-(int lhs, const float_type& rhs) {
    return static_cast<float_type>(lhs) - rhs;
}

float_type operator*(const float_type& lhs, int rhs) {
    return lhs * static_cast<float_type>(rhs);
}

float_type operator*(int lhs, const float_type& rhs) {
    return static_cast<float_type>(lhs) * rhs;
}

float_type operator/(const float_type& lhs, int rhs) {
    return lhs / static_cast<float_type>(rhs);
}

float_type operator/(int lhs, const float_type& rhs) {
    return static_cast<float_type>(lhs) / rhs;
}
