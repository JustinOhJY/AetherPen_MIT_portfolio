#include "quaternions.h"

quaternion::quaternion() {
    qw = 1.0;
    qx = 0.0;
    qy = 0.0;
    qz = 0.0;
    angle = 0;
}

quaternion::quaternion(float gx, float gy, float gz, float time) {
    qx = gx;
    qy = gy;
    qz = gz;
    dt = time;
    float gyro_norm = sqrt(pow(gx, 2) + pow(gy, 2) + pow(gz, 2));
    q_norm = gyro_norm;
    angle = (gyro_norm * dt) / 2.0;
    qw = cos(angle);
    normalize();
}

quaternion::quaternion(const quaternion &q) {
    qw = q.qw;
    qx = q.qx;
    qy = q.qy;
    qz = q.qz;
    angle = q.angle;
}

void quaternion::normalize() {
    float factor = sin(acos(qw));
    float norm = sqrt(pow(qx, 2) + pow(qy, 2) + pow(qz, 2));
    if (norm >= 1e-12) {
        qx /= norm;
        qy /= norm;
        qz /= norm;
        qx *= factor;
        qy *= factor;
        qz *= factor;
    }
    else {
        qw = 1;
        qx = 0;
        qy = 0;
        qz = 0;
    }
}

quaternion quaternion::operator*(const quaternion &q1) {
    quaternion q_output;
    q_output.qw = qw * q1.qw - qx * q1.qx - qy * q1.qy - qz * q1.qz;
    q_output.qx = qw * q1.qx + qx * q1.qw + qy * q1.qz - qz * q1.qy;
    q_output.qy = qw * q1.qy - qx * q1.qz + qy * q1.qw + qz * q1.qx;
    q_output.qz = qw * q1.qz + qx * q1.qy - qy * q1.qx + qz * q1.qw;
    return q_output;
}

quaternion quaternion::inv() {
    quaternion inverse;
    inverse.qw = qw;
    inverse.qx = -qx;
    inverse.qy = -qy;
    inverse.qz = -qz;
    return inverse;
}
