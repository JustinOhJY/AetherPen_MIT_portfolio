#ifndef QUATERNIONS_H
#define QUATERNIONS_H

#include <Arduino.h>
#include <math.h>

struct quaternion {
    float qw;
    float qx;
    float qy;
    float qz;
    float dt;
    float q_norm;
    float angle;
    
    quaternion();
    
    quaternion(float gx, float gy, float gz, float time);
    
    quaternion(const quaternion &q);
    
    void normalize();
    quaternion operator*(const quaternion &q1);
    quaternion inv();
};

#endif
