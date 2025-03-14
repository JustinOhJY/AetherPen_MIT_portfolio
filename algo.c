static void quat_normalize(float q[4]) {
    float norm = sqrtf(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
    q[0] /= norm;
    q[1] /= norm;
    q[2] /= norm;
    q[3] /= norm;
}

void ekf_init(ekf_t* ekf, float sigma_wx, float sigma_wy, float sigma_wz, float sigma_ax, float sigma_ay, float sigma_az) {
    // State x = [q0 q1 q2 q3 bwx bwy bwz]
    ekf->x[0] = 1.0f;
    ekf->x[1] = 0.0f;
    ekf->x[2] = 0.0f;
    ekf->x[3] = 0.0f;


    // Gyro spectral noise covariance matrix
    ekf->sigma_w_data[0] = sigma_wx * sigma_wx;
    ekf->sigma_w_data[1] = 0.0f;
    ekf->sigma_w_data[2] = 0.0f;
    ekf->sigma_w_data[3] = 0.0f;
    ekf->sigma_w_data[4] = sigma_wy * sigma_wy;
    ekf->sigma_w_data[5] = 0.0f;
    ekf->sigma_w_data[6] = 0.0f;
    ekf->sigma_w_data[7] = 0.0f;
    ekf->sigma_w_data[8] = sigma_wz * sigma_wz;

    // Measurement noise covariance matrix
    ekf->R_data[0] = sigma_ax * sigma_ax;
    ekf->R_data[1] = 0.0f;
    ekf->R_data[2] = 0.0f;
    ekf->R_data[3] = 0.0f;
    ekf->R_data[4] = sigma_ay * sigma_ay;
    ekf->R_data[5] = 0.0f;
    ekf->R_data[6] = 0.0f;
    ekf->R_data[7] = 0.0f;
    ekf->R_data[8] = sigma_az * sigma_az;

    
}

void ekf_predict(ekf_t *ekf, float w_x, float w_y, float w_z, float dt) {
    // State prediction, xp: # Math: \hat{q}_{t}
        // # Math: \Omega_t = \begin{bmatrix} 0 & -\omega_x & -\omega_y & -\omega_z \\ \omega_x & 0 & \omega_z & -\omega_y \\ \omega_y & -\omega_z & 0 & \omega_x \\ \omega_z & \omega_y & -\omega_x & 0 \end{bmatrix}
        // # Math: \hat{q}_{t} = q_{t-1} + \frac{1}{2} \Omega_t q_{t-1} \cdot \Delta t
        // # Math: \begin{bmatrix} \hat{q}_w \\ \hat{q}_x \\ \hat{q}_y \\ \hat{q}_z \end{bmatrix} = \begin{bmatrix} q_w - \frac{\Delta t}{2}\omega_x q_x - \frac{\Delta t}{2}\omega_y q_y - \frac{\Delta t}{2}\omega_z q_z \\ q_x + \frac{\Delta t}{2}\omega_x q_w - \frac{\Delta t}{2}\omega_y q_z + \frac{\Delta t}{2}\omega_z q_y \\ q_y + \frac{\Delta t}{2}\omega_x q_z + \frac{\Delta t}{2}\omega_y q_w - \frac{\Delta t}{2}\omega_z q_x \\ q_z - \frac{\Delta t}{2}\omega_x q_y + \frac{\Delta t}{2}\omega_y q_x + \frac{\Delta t}{2}\omega_z q_w \end{bmatrix}
        
    ekf->xp[0] = ekf->x[0] - 0.5f * (w_x * ekf->x[1] + w_y * ekf->x[2] + w_z * ekf->x[3]) * dt;
    ekf->xp[1] = ekf->x[1] + 0.5f * (w_x * ekf->x[0] - w_y * ekf->x[3] + w_z * ekf->x[2]) * dt;
    ekf->xp[2] = ekf->x[2] + 0.5f * (w_x * ekf->x[3] + w_y * ekf->x[0] - w_z * ekf->x[1]) * dt;
    ekf->xp[3] = ekf->x[3] + 0.5f * (-w_x * ekf->x[2] + w_y * ekf->x[1] + w_z * ekf->x[0]) * dt;
    
    quat_normalize(ekf->xp);

    // State Transition Matrix, F
        // # Math: F(q_{t-1}, \omega_t) = \frac{\partial f(q_{t-1}, \omega_t)}{\partial q} = \begin{bmatrix} \frac{\partial f(q_{t-1}, \omega_t)}{\partial q_w} & \frac{\partial f(q_{t-1}, \omega_t)}{\partial q_x} & \frac{\partial f(q_{t-1}, \omega_t)}{\partial q_y} & \frac{\partial f(q_{t-1}, \omega_t)}{\partial q_z} \end{bmatrix} = \begin{bmatrix} 1 & -\frac{\Delta t}{2}\omega_x & -\frac{\Delta t}{2}\omega_y & -\frac{\Delta t}{2}\omega_z \\ \frac{\Delta t}{2}\omega_x & 1 & \frac{\Delta t}{2}\omega_z & -\frac{\Delta t}{2}\omega_y \\ \frac{\Delta t}{2}\omega_y & -\frac{\Delta t}{2}\omega_z & 1 & \frac{\Delta t}{2}\omega_x \\ \frac{\Delta t}{2}\omega_z & \frac{\Delta t}{2}\omega_y & -\frac{\Delta t}{2}\omega_x & 1 \end{bmatrix}
    arm_matrix_instance_f32 F;
    float32_t F_data[16] = {
        1.0f,            -0.5f * w_x * dt, -0.5f * w_y * dt, -0.5f * w_z * dt,
        0.5f * w_x * dt,  1.0f,             0.5f * w_z * dt, -0.5f * w_y * dt,
        0.5f * w_y * dt, -0.5f * w_z * dt,  1.0f,             0.5f * w_x * dt,
        0.5f * w_z * dt,  0.5f * w_y * dt, -0.5f * w_x * dt,  1.0f
    };
    arm_mat_init_f32(&F, 4, 4, F_data);
    
    float32_t Ft_data[16];
    arm_matrix_instance_f32 Ft;
    arm_mat_init_f32(&Ft, 4, 4, Ft_data);
    arm_mat_trans_f32(&F, &Ft);

    // Process Noise Covariance, Q: # Math: Q_t = W_t \Sigma_\omega W_t^T
        // # Math: W_t = \frac{\partial f(q_{t-1}, \omega_t)}{\partial \omega} = \begin{bmatrix} \frac{\partial f(q_{t-1}, \omega_t)}{\partial \omega_x} & \frac{\partial f(q_{t-1}, \omega_t)}{\partial \omega_y} & \frac{\partial f(q_{t-1}, \omega_t)}{\partial \omega_z} \end{bmatrix} = \frac{\Delta t}{2} \begin{bmatrix} -q_x & -q_y & -q_z \\ q_w & -q_z & q_y \\ q_z & q_w & -q_x \\ -q_y & q_x & q_w \end{bmatrix}
        // # Math: \Sigma_\omega = \begin{bmatrix} \sigma_{\omega x}^2 & 0 & 0 \\ 0 & \sigma_{\omega y}^2 & 0 \\ 0 & 0 & \sigma_{\omega z}^2 \end{bmatrix}
    arm_matrix_instance_f32 W;
    float32_t W_data[12] = {
        -ekf->xp[1] * dt / 2, -ekf->xp[2] * dt / 2, -ekf->xp[3] * dt / 2,
         ekf->xp[0] * dt / 2, -ekf->xp[3] * dt / 2,  ekf->xp[2] * dt / 2,
         ekf->xp[3] * dt / 2,  ekf->xp[0] * dt / 2, -ekf->xp[1] * dt / 2,
        -ekf->xp[2] * dt / 2,  ekf->xp[1] * dt / 2,  ekf->xp[0] * dt / 2
    };
    arm_mat_init_f32(&W, 4, 3, W_data);

    float32_t Wt_data[12];
    arm_matrix_instance_f32 Wt;
    arm_mat_init_f32(&Wt, 3, 4, Wt_data);
    arm_mat_trans_f32(&W, &Wt);

    arm_matrix_instance_f32 sigma_w;
    arm_mat_init_f32(&sigma_w, 3, 3, ekf->sigma_w_data);
    
    float32_t Wsigma_data[16];
    arm_matrix_instance_f32 Wsigma;
    arm_mat_init_f32(&Wsigma, 4, 4, Wsigma_data);
    arm_mat_mult_f32(&W, &sigma_w, &Wsigma);

    float32_t Q_data[16];
    arm_matrix_instance_f32 Q;
    arm_mat_init_f32(&Q, 4, 4, Q_data);
    arm_mat_mult_f32(&Wsigma, &Wt, &Q);

    // Predicted Covariance Matrix, Pp: # Math: \hat{P}_t = F_t P_{t-1} F_t^T + Q_t
    arm_matrix_instance_f32 P;
    arm_mat_init_f32(&P, 4, 4, ekf->P_data);

    float32_t FP_data[16];
    arm_matrix_instance_f32 FP;
    arm_mat_init_f32(&FP, 4, 4, FP_data);
    arm_mat_mult_f32(&F, &P, &FP);

    float32_t FPFt_data[16];
    arm_matrix_instance_f32 FPFt;
    arm_mat_init_f32(&FPFt, 4, 4, FPFt_data);
    arm_mat_mult_f32(&FP, &Ft, &FPFt);

    arm_matrix_instance_f32 Pp;
    arm_mat_init_f32(&Pp, 4, 4, ekf->Pp_data);
    arm_mat_add_f32(&FPFt, &Q, &Pp);


    #if 0
    // Update the state and covariance matrix
    for (int i = 0; i < 4; i++) {
        ekf->x[i] = ekf->xp[i];
        for (int j = 0; j < 4; j++) {
            ekf->P_data[i * 4 + j] = ekf->Pp_data[i * 4 + j];
        }
    }
    #endif

}

void ekf_update(ekf_t *ekf, float a_x, float a_y, float a_z) {
    // Correction step
    // Corrected state, q_t: # Math: q_t = \hat{q}_t + K_t (z_t - h(q_t))
        // Measurement vector, z_t: # Math: z_t = \begin{bmatrix} a_x \\ a_y \\ a_z \end{bmatrix}
    // Measurement model, h(q_t): # Math: h(\hat{q}_t) = \hat{a}  = C(\hat{q})^T g
        // Rotation Matrix, C # Math: C(\hat{q}) = 2 \begin{bmatrix} g_x(\frac{1}{2} - q_y^2 - q_z^2) + g_y(q_w q_z + q_x q_y) + g_z(q_x q_z - q_w q_y) \\ g_x(q_x q_y - q_w q_z) + g_y(\frac{1}{2} - q_x^2 - q_z^2) + g_z(q_w q_x + q_y q_z) \\ g_x(q_w q_y + q_x q_z) + g_y(q_y q_z - q_w q_x) + g_z(\frac{1}{2} - q_x^2 - q_y^2) \end{bmatrix}
        // # Math: g = \begin{bmatrix} 0 \\ 0 \\ 1 \end{bmatrix}
        // # Math: \hat{a} = \begin{bmatrix} 2(q_x q_z - q_w q_y) \\ 2(q_w q_x + q_y q_z) \\ q_w^2 - q_x^2 - q_y^2 + q_z^2 \end{bmatrix}
        // For magnetometer TODO: # Math: h(\hat{q}_t) = \begin{bmatrix} \hat{a} \\ \hat{m} \end{bmatrix} = \begin{bmatrix} C(\hat{q})^T g \\ C(\hat{q})^T r \end{bmatrix}
            // r being earth's magnetic field as 3D vec # Math: 2 \begin{bmatrix} g_x(\frac{1}{2} - q_y^2 - q_z^2) + g_y(q_w q_z + q_x q_y) + g_z(q_x q_z - q_w q_y) \\ g_x(q_x q_y - q_w q_z) + g_y(\frac{1}{2} - q_x^2 - q_z^2) + g_z(q_w q_x + q_y q_z) \\ g_x(q_w q_y + q_x q_z) + g_y(q_y q_z - q_w q_x) + g_z(\frac{1}{2} - q_x^2 - q_y^2) \\ r_x(\frac{1}{2} - q_y^2 - q_z^2) + r_y(q_w q_z + q_x q_y) + r_z(q_x q_z - q_w q_y) \\ r_x(q_x q_y - q_w q_z) + r_y(\frac{1}{2} - q_x^2 - q_z^2) + r_z(q_w q_x + q_y q_z) \\ r_x(q_w q_y + q_x q_z) + r_y(q_y q_z - q_w q_x) + r_z(\frac{1}{2} - q_x^2 - q_y^2) \end{bmatrix}
    // Innovation, v_t: # Math: v_t = z_t - h(\hat{q}_t)
    // normalize acceleration
    float32_t a_norm = sqrtf(a_x * a_x + a_y * a_y + a_z * a_z);
    a_x /= a_norm;
    a_y /= a_norm;
    a_z /= a_norm;
    float32_t z_data[3] = {a_x, a_y, a_z};
    arm_matrix_instance_f32 z;
    arm_mat_init_f32(&z, 3, 1, z_data);
    
    // construct h(q_t)
    float32_t h_data[3] = {
        2.0f * (ekf->xp[1] * ekf->xp[3] - ekf->xp[0] * ekf->xp[2]),
        2.0f * (ekf->xp[0] * ekf->xp[1] + ekf->xp[2] * ekf->xp[3]),
        ekf->xp[0] * ekf->xp[0] - ekf->xp[1] * ekf->xp[1] - ekf->xp[2] * ekf->xp[2] + ekf->xp[3] * ekf->xp[3]
    };
    arm_matrix_instance_f32 h;
    arm_mat_init_f32(&h, 3, 1, h_data);

    float32_t v_data[3];
    arm_matrix_instance_f32 v;
    arm_mat_init_f32(&v, 3, 1, v_data);
    arm_mat_sub_f32(&z, &h, &v);

    // Jacobian of the measurement model, H: # Math: H(q) = \begin{bmatrix} -2q_y & 2q_z & -2q_w & 2q_x \\ 2q_x & 2q_w & 2q_z & 2q_y \\ 2q_w & -2q_x & -2q_y & 2q_z \end{bmatrix}
    float32_t H_data[12] = {
        -2.0f * ekf->xp[2],  2.0f * ekf->xp[3], -2.0f * ekf->xp[0],  2.0f * ekf->xp[1],
         2.0f * ekf->xp[1],  2.0f * ekf->xp[0],  2.0f * ekf->xp[3],  2.0f * ekf->xp[2],
         2.0f * ekf->xp[0], -2.0f * ekf->xp[1], -2.0f * ekf->xp[2],  2.0f * ekf->xp[3]
    };
    arm_matrix_instance_f32 H;
    arm_mat_init_f32(&H, 3, 4, H_data);

    // Measurement Noise Covariance, R: # Math: R_t = \begin{bmatrix} \sigma_{a x}^2 & 0 & 0 \\ 0 & \sigma_{a y}^2 & 0 \\ 0 & 0 & \sigma_{a z}^2 \end{bmatrix}
    arm_matrix_instance_f32 R;
    arm_mat_init_f32(&R, 3, 3, ekf->R_data);

    // Kalman Gain, K: # Math: K_t = \hat{P}_t H_t^T (H_t \hat{P}_t H_t^T + R_t)^{-1}
    arm_matrix_instance_f32 Pp;
    arm_mat_init_f32(&Pp, 4, 4, ekf->Pp_data);

    float32_t Ht_data[12];
    arm_matrix_instance_f32 Ht;
    arm_mat_init_f32(&Ht, 4, 3, Ht_data);
    arm_mat_trans_f32(&H, &Ht);

    float32_t PHt_data[12];
    arm_matrix_instance_f32 PHt;
    arm_mat_init_f32(&PHt, 4, 3, PHt_data);
    arm_mat_mult_f32(&Pp, &Ht, &PHt);

    // # Math: S_t = H_t \hat{P}_t H_t^T + R_t
    float32_t S_data[9];
    arm_matrix_instance_f32 S;
    arm_mat_init_f32(&S, 3, 3, S_data);
    arm_mat_mult_f32(&H, &PHt, &S);
    arm_mat_add_f32(&S, &R, &S);

    float32_t Si_data[9];   
    arm_matrix_instance_f32 Si;
    arm_mat_init_f32(&Si, 3, 3, Si_data);
    arm_mat_inverse_f32(&S, &Si);

    float32_t K_data[12];
    arm_matrix_instance_f32 K;
    arm_mat_init_f32(&K, 4, 3, K_data);
    arm_mat_mult_f32(&PHt, &Si, &K);

    // Update the state and covariance matrix
        // # Math: \hat{q}_t = \hat{q}_t + K_t v_t
        // # Math: \hat{P}_t = (I - K_t H_t) \hat{P}_t  

    float32_t Kv_data[4];
    arm_matrix_instance_f32 Kv;
    arm_mat_init_f32(&Kv, 4, 1, Kv_data);
    arm_mat_mult_f32(&K, &v, &Kv);

    for (int i = 0; i < 4; i++) {
        ekf->x[i] = ekf->xp[i] + Kv_data[i];
    }

    // Compute I - K*H
    float32_t I_KH_data[16];
    arm_matrix_instance_f32 I_KH;
    arm_mat_init_f32(&I_KH, 4, 4, I_KH_data);
    float32_t KH_data[16];
    arm_matrix_instance_f32 KH;
    arm_mat_init_f32(&KH, 4, 4, KH_data);
    arm_mat_mult_f32(&K, &H, &KH);

    // Create identity matrix I (4x4)
    float32_t I_data[16] = {
        1,0,0,0,
        0,1,0,0,
        0,0,1,0,
        0,0,0,1
    };
    arm_matrix_instance_f32 I;
    arm_mat_init_f32(&I, 4, 4, I_data);

    // Compute I - KH
    for (int i = 0; i < 16; i++) {
        I_KH_data[i] = I_data[i] - KH_data[i];
    }

    // Now compute new covariance: P = (I - KH)*Pp
    float32_t P_new_data[16];
    arm_matrix_instance_f32 P_new;
    arm_mat_init_f32(&P_new, 4, 4, P_new_data);
    arm_mat_mult_f32(&I_KH, &Pp, &P_new);

    // Copy P_new_data to ekf->P_data.
    for (int i = 0; i < 16; i++) {
        ekf->P_data[i] = P_new_data[i];
    }

    quat_normalize(ekf->x);

}

void ekf_process(ekf_t *ekf, float w_x, float w_y, float w_z, float a_x, float a_y, float a_z, float dt) {
    printf("%f %f %f %f\n", ekf->x[0], ekf->x[1], ekf->x[2], ekf->x[3]);
    #if 0
    float omega_norm = sqrtf(w_x * w_x + w_y * w_y + w_z * w_z);
    if (omega_norm < 0.005)
    {
        // No significant rotation; keep the quaternion as is
        return;
    }
    #endif
    ekf_predict(ekf, w_x, w_y, w_z, dt);
    ekf_update(ekf, a_x, a_y, a_z);

}
