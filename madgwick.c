void madgwick_filter(float accel[3], float gyro[3], float *q_a, float *q_b, float *q_c, float *q_d, float dt_s)
{
    float recipNorm;
    float s0, s1, s2, s3;
    float qDot0, qDot1, qDot2, qDot3;
    float _2q0, _2q1, _2q2, _2q3;
    float q0q0, q1q1, q2q2, q3q3;

    // Normalize the accelerometer measurement
    float ax = accel[0], ay = accel[1], az = accel[2];
    recipNorm = sqrtf(ax * ax + ay * ay + az * az);
    if (recipNorm < 1e-8f)
        return; // avoid division by zero
    recipNorm = 1.0f / recipNorm;
    ax *= recipNorm;
    ay *= recipNorm;
    az *= recipNorm;

    _2q0 = 2.0f * (*q_a);
    _2q1 = 2.0f * (*q_b);
    _2q2 = 2.0f * (*q_c);
    _2q3 = 2.0f * (*q_d);
    q0q0 = (*q_a) * (*q_a);
    q1q1 = (*q_b) * (*q_b);
    q2q2 = (*q_c) * (*q_c);
    q3q3 = (*q_d) * (*q_d);

    // Compute the objective function (error between estimated and measured gravity)
    // The estimated gravity direction from the quaternion is:
    // # Math: v = \begin{bmatrix} 2(q_1q_3 - q_0q_2) \\ 2(q_0q_1 + q_2q_3) \\ q_0q_0 - q_1q_1 - q_2q_2 + q_3q_3 \end{bmatrix}
    float f0 = _2q1 * (*q_d) - _2q0 * (*q_c) - ax;
    float f1 = _2q0 * (*q_b) + _2q2 * (*q_d) - ay;
    float f2 = q0q0 - q1q1 - q2q2 + q3q3 - az;

    // Compute the gradient (matrix multiplication of Jacobian^T and f)
    s0 = -_2q2 * f0 + _2q1 * f1;
    s1 = _2q3 * f0 + _2q0 * f1 - 4.0f * (*q_b) * f2;
    s2 = -_2q0 * f0 + _2q3 * f1 - 4.0f * (*q_c) * f2;
    s3 = _2q1 * f0 + _2q2 * f1;

    // Normalize
    recipNorm = sqrtf(s0 * s0 + s1 * s1 + s2 * s2 + s3 * s3);
    if (recipNorm < 1e-8f)
        return;
    recipNorm = 1.0f / recipNorm;
    s0 *= recipNorm;
    s1 *= recipNorm;
    s2 *= recipNorm;
    s3 *= recipNorm;

    // quaternion derivative from gyroscope
    float gx = gyro[0], gy = gyro[1], gz = gyro[2];
    qDot0 = 0.5f * (-(*q_b) * gx - (*q_c) * gy - (*q_d) * gz) - BETA * s0;
    qDot1 = 0.5f * ((*q_a) * gx + (*q_c) * gz - (*q_d) * gy) - BETA * s1;
    qDot2 = 0.5f * ((*q_a) * gy - (*q_b) * gz + (*q_d) * gx) - BETA * s2;
    qDot3 = 0.5f * ((*q_a) * gz + (*q_b) * gy - (*q_c) * gx) - BETA * s3;

    // Integrate to get new quaternion
    *q_a += qDot0 * dt_s;
    *q_b += qDot1 * dt_s;
    *q_c += qDot2 * dt_s;
    *q_d += qDot3 * dt_s;

    // Normalize
    recipNorm = sqrtf((*q_a) * (*q_a) + (*q_b) * (*q_b) + (*q_c) * (*q_c) + (*q_d) * (*q_d));
    if (recipNorm < 1e-8f)
        return;
    recipNorm = 1.0f / recipNorm;
    *q_a *= recipNorm;
    *q_b *= recipNorm;
    *q_c *= recipNorm;
    *q_d *= recipNorm;
}