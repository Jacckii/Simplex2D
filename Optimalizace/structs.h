#pragma once

struct Point {
    float x;
    float y;
   
    Point(float x_, float y_) : x(x_), y(y_) {}
    Point() : x(0), y(0) {};

    Point operator*(float value) {
        return Point(this->x * value, this->y * value);
    }
    Point operator+(Point b) {
        return Point(this->x + b.x, this->y + b.y);
    }
    Point operator-(Point b) {
        return Point(this->x - b.x, this->y - b.y);
    }
};  

struct Simplex2D {
    Point x1;
    Point x2;
    Point x3;

    Simplex2D(Point a_, Point b_, Point c_) : x1(a_), x2(b_), x3(c_) {};
    Simplex2D() : x1(0, 0), x2(0, 0), x3(0, 0) {};
};

struct SimplexRecord {
    Simplex2D simplex;
    float delta;
    int iteration_index;
};