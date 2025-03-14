#include <Adafruit_ICM20X.h>
#include <Adafruit_ICM20948.h>
#include <Adafruit_Sensor.h>
#include <Wire.h>
#include <math.h>
#include <quaternions.h>
#include <linalg_funcs.h>
#include <other_funcs.h>

// EKF state structure: quaternion state (x), predicted state (xp),
// covariance matrices (P_data, Pp_data), and noise covariances.
typedef struct {
  float x[4];          // Current quaternion state [q0, q1, q2, q3]
  float xp[4];         // Predicted state (quaternion)
  float P_data[16] = {0};    // Covariance matrix (4x4)
  float Pp_data[16];   // Predicted covariance (4x4)
  float sigma_w_data[9] = {0}; // Gyro noise covariance (3x3)
  float R_data[36] = {0};     // Measurement noise covariance (3x3)
} ekf_t;

ekf_t ekf;

Adafruit_ICM20948 icm;

const float gravity = 10.015;

unsigned long t1 = 0, t2 = 0;
int counter = 0;
const int limit = 10;  

float east[3] = {0};
float north[3] = {0};

bool adjusted = false;
bool adjust_start = true;
int stationary_counter = 0;

// offset corrected values
float a[3] = {0};
float g[3] = {0};
float m[3] = {0};

// offset values (for accelerometer and gyro)
float ao[3] = {0};
float go[3] = {0};
float mo[3] = {0};

//global quaternions
struct quaternion q_r;
struct quaternion Q;

bool is_stationary(int counter, float (&arr)[3]) {
  if (abs(get_norm(arr) - gravity) < 0.15) {
    if (counter == 10) {
      return true;
    }
    else {
      counter += 1;
      return false;
    }
  }
  else {
    counter = 0;
    return false;
  }
}

//------ Matrix Math Helper Functions ------//

// Multiply matrix A (m x n) by matrix B (n x p) and store the result in C (m x p)
void matMultiply(const float *A, const float *B, float *C, int m, int n, int p) {
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < p; j++) {
      C[i * p + j] = 0.0f;
      for (int k = 0; k < n; k++) {
        C[i * p + j] += A[i * n + k] * B[k * p + j];
      }
    }
  }
}

// Compute the transpose of matrix A (m x n) and store in At (n x m)
void matTranspose(const float *A, float *At, int m, int n) {
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      At[j * m + i] = A[i * n + j];
    }
  }
}

// Add two matrices A and B (both m x n) and store the result in C (m x n)
void matAdd(const float *A, const float *B, float *C, int m, int n) {
  int total = m * n;
  for (int i = 0; i < total; i++) {
    C[i] = A[i] + B[i];
  }
}

// Subtract matrix B from A (both m x n) and store the result in C (m x n)
void matSub(const float *A, const float *B, float *C, int m, int n) {
  int total = m * n;
  for (int i = 0; i < total; i++) {
    C[i] = A[i] - B[i];
  }
}

// Invert a 3x3 matrix A; returns 0 if successful, -1 if singular.
int matInverse3x3(const float *A, float *invA) {
  float a11 = A[0], a12 = A[1], a13 = A[2];
  float a21 = A[3], a22 = A[4], a23 = A[5];
  float a31 = A[6], a32 = A[7], a33 = A[8];
  
  float det = a11*(a22*a33 - a23*a32) - a12*(a21*a33 - a23*a31) + a13*(a21*a32 - a22*a31);
  if (fabs(det) < 1e-6) return -1;
  float invDet = 1.0f / det;
  
  invA[0] = (a22 * a33 - a23 * a32) * invDet;
  invA[1] = (a13 * a32 - a12 * a33) * invDet;
  invA[2] = (a12 * a23 - a13 * a22) * invDet;
  invA[3] = (a23 * a31 - a21 * a33) * invDet;
  invA[4] = (a11 * a33 - a13 * a31) * invDet;
  invA[5] = (a13 * a21 - a11 * a23) * invDet;
  invA[6] = (a21 * a32 - a22 * a31) * invDet;
  invA[7] = (a12 * a31 - a11 * a32) * invDet;
  invA[8] = (a11 * a22 - a12 * a21) * invDet;
  
  return 0;
}

int matInverse6x6(const float *A, float *Ainv)
{
    // Copy A into a working matrix M, and initialize I to identity
    float M[36], I[36];
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            M[i*6 + j] = A[i*6 + j];
            I[i*6 + j] = (i == j) ? 1.0f : 0.0f;
        }
    }

    // Forward elimination with partial pivoting
    for (int c = 0; c < 6; c++) {
        // 1) Find pivot row
        int pivot = c;
        float maxAbs = fabsf(M[c*6 + c]);
        for (int r = c + 1; r < 6; r++) {
            float val = fabsf(M[r*6 + c]);
            if (val > maxAbs) {
                pivot = r;
                maxAbs = val;
            }
        }

        // 2) Check for near-singular
        if (maxAbs < 1e-6f) {
            // Matrix is singular or close to it
            return -1;
        }

        // 3) Swap pivot row with current row c if needed
        if (pivot != c) {
            for (int k = 0; k < 6; k++) {
                float tmpM = M[c*6 + k];
                M[c*6 + k] = M[pivot*6 + k];
                M[pivot*6 + k] = tmpM;

                float tmpI = I[c*6 + k];
                I[c*6 + k] = I[pivot*6 + k];
                I[pivot*6 + k] = tmpI;
            }
        }

        // 4) Normalize pivot row (so pivot element = 1)
        float pivotVal = M[c*6 + c];
        for (int k = 0; k < 6; k++) {
            M[c*6 + k] /= pivotVal;
            I[c*6 + k] /= pivotVal;
        }

        // 5) Eliminate below pivot
        for (int r = c + 1; r < 6; r++) {
            float factor = M[r*6 + c];
            for (int k = 0; k < 6; k++) {
                M[r*6 + k] -= factor * M[c*6 + k];
                I[r*6 + k] -= factor * I[c*6 + k];
            }
        }
    }

    // Back substitution to get reduced row echelon form
    for (int c = 5; c >= 0; c--) {
        for (int r = c - 1; r >= 0; r--) {
            float factor = M[r*6 + c];
            for (int k = 0; k < 6; k++) {
                M[r*6 + k] -= factor * M[c*6 + k];
                I[r*6 + k] -= factor * I[c*6 + k];
            }
        }
    }

    // Now I[] is the inverse of A[]
    for (int i = 0; i < 36; i++) {
        Ainv[i] = I[i];
    }

    return 0; // success
}

//------ Quaternion and EKF Functions ------//

// Normalize a quaternion (4-element vector)
void quat_normalize(float q[4]) {
  float norm = sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
  if (norm > 0.0f) {
    q[0] /= norm;
    q[1] /= norm;
    q[2] /= norm;
    q[3] /= norm;
  }
}

// Initialize the EKF structure and noise covariance matrices
void ekf_init(ekf_t* ekf, float sigma_wx, float sigma_wy, float sigma_wz, 
              float sigma_ax, float sigma_ay, float sigma_az, float sigma_mx, float sigma_my, float sigma_mz) {
  // Set initial quaternion state (no rotation)
  ekf->x[0] = 1.0f;
  ekf->x[1] = 0.0f;
  ekf->x[2] = 0.0f;
  ekf->x[3] = 0.0f;

  // Initialize covariance P_data to identity (4x4)
  ekf->P_data[0]  = 1.0f;
  ekf->P_data[5]  = 1.0f;
  ekf->P_data[10] = 1.0f;
  ekf->P_data[15] = 1.0f;
  
  // Set gyro spectral noise covariance matrix (3x3)
  ekf->sigma_w_data[0] = sigma_wx * sigma_wx;
  ekf->sigma_w_data[4] = sigma_wy * sigma_wy;
  ekf->sigma_w_data[8] = sigma_wz * sigma_wz;
  
  // Set measurement noise covariance matrix (6x6)
  ekf->R_data[0] = sigma_ax * sigma_ax;
  ekf->R_data[7] = sigma_ay * sigma_ay;
  ekf->R_data[14] = sigma_az * sigma_az;
  ekf->R_data[21] = sigma_mx * sigma_mx;
  ekf->R_data[28] = sigma_my * sigma_my;
  ekf->R_data[35] = sigma_mz * sigma_mz;
}

// Predict step: update quaternion state and covariance based on gyro data.
void ekf_predict(ekf_t *ekf, float w_x, float w_y, float w_z, float dt) {
  // State prediction for quaternion (xp = x + 0.5 * Omega * x * dt)
  ekf->xp[0] = ekf->x[0] - 0.5f * (w_x * ekf->x[1] + w_y * ekf->x[2] + w_z * ekf->x[3]) * dt;
  ekf->xp[1] = ekf->x[1] + 0.5f * (w_x * ekf->x[0] - w_y * ekf->x[3] + w_z * ekf->x[2]) * dt;
  ekf->xp[2] = ekf->x[2] + 0.5f * (w_x * ekf->x[3] + w_y * ekf->x[0] - w_z * ekf->x[1]) * dt;
  ekf->xp[3] = ekf->x[3] + 0.5f * (-w_x * ekf->x[2] + w_y * ekf->x[1] + w_z * ekf->x[0]) * dt;
  quat_normalize(ekf->xp);

  // Build the state transition matrix F (4x4)
  float F[16] = {
      1.0f,            -0.5f * w_x * dt, -0.5f * w_y * dt, -0.5f * w_z * dt,
      0.5f * w_x * dt,  1.0f,             0.5f * w_z * dt, -0.5f * w_y * dt,
      0.5f * w_y * dt, -0.5f * w_z * dt,  1.0f,             0.5f * w_x * dt,
      0.5f * w_z * dt,  0.5f * w_y * dt, -0.5f * w_x * dt,  1.0f
  };
  float Ft[16];
  matTranspose(F, Ft, 4, 4);

  // Process noise covariance Q = W * sigma_w * Wt
  // Construct W (4x3)
  float W[12] = {
      -ekf->xp[1] * dt / 2, -ekf->xp[2] * dt / 2, -ekf->xp[3] * dt / 2,
       ekf->xp[0] * dt / 2, -ekf->xp[3] * dt / 2,  ekf->xp[2] * dt / 2,
       ekf->xp[3] * dt / 2,  ekf->xp[0] * dt / 2, -ekf->xp[1] * dt / 2,
      -ekf->xp[2] * dt / 2,  ekf->xp[1] * dt / 2,  ekf->xp[0] * dt / 2
  };
  // Compute intermediate: Wsigma = W * sigma_w (4x3 * 3x3 = 4x3)
  float Wsigma[12];
  matMultiply(W, ekf->sigma_w_data, Wsigma, 4, 3, 3);
  // Compute W transpose, Wt (3x4)
  float Wt[12];
  matTranspose(W, Wt, 4, 3);
  // Compute Q = Wsigma * Wt (4x3 * 3x4 = 4x4)
  float Q[16];
  matMultiply(Wsigma, Wt, Q, 4, 3, 4);

  // Predicted covariance: Pp = F * P * Ft + Q
  float FP[16];
  matMultiply(F, ekf->P_data, FP, 4, 4, 4);
  float FPFt[16];
  matMultiply(FP, Ft, FPFt, 4, 4, 4);
  float Pp[16];
  matAdd(FPFt, Q, Pp, 4, 4);
  for (int i = 0; i < 16; i++) {
    ekf->Pp_data[i] = Pp[i];
  }
}

// Update step: correct the predicted state with accelerometer AND mag data.
void ekf_update(ekf_t *ekf, float a_x, float a_y, float a_z, float m_x, float m_y, float m_z) {
  // Normalize accelerometer measurement (assumed to measure gravity)
  float a_norm = sqrt(a_x*a_x + a_y*a_y + a_z*a_z);
  if (a_norm > 0.0f) {
    a_x /= a_norm;
    a_y /= a_norm;
    a_z /= a_norm;
  }

  float m_norm = sqrt(m_x*m_x + m_y*m_y + m_z*m_z);
  if (m_norm > 0.0f) {
    m_x /= m_norm;
    m_y /= m_norm;
    m_z /= m_norm;
  }
  // Measurement vector z (6x1)
  float z[6] = { a_x, a_y, a_z, m_x, m_y, m_z};

  float r[3] = {0, cos(73.0f * M_PI/180.0f), -sin(73.0f * M_PI / 180.0f)};

  // Construct h(xp) from predicted quaternion (xp)
  float h[6] = {
    2.0f * (ekf->xp[1] * ekf->xp[3] - ekf->xp[0] * ekf->xp[2]),
    2.0f * (ekf->xp[0] * ekf->xp[1] + ekf->xp[2] * ekf->xp[3]),
    2.0f * (0.5f - ekf->xp[1] * ekf->xp[1] - ekf->xp[2] * ekf->xp[2]),
    2.0f * (r[0] * (0.5f - ekf->xp[2] * ekf->xp[2] - ekf->xp[3] * ekf->xp[3]) + r[1] * (ekf->xp[0] * ekf->xp[3] + ekf->xp[1] * ekf->xp[2]) + r[2] * (ekf->xp[1] * ekf->xp[3] - ekf->xp[0] * ekf->xp[2])),
    2.0f * (r[0] * (ekf->xp[1] * ekf->xp[2] - ekf->xp[0] * ekf->xp[3]) + r[1] * (0.5f - ekf->xp[1] * ekf->xp[1] - ekf->xp[3] * ekf->xp[3]) + r[2] * (ekf->xp[0] * ekf->xp[1] + ekf->xp[2] * ekf->xp[3])),
    2.0f * (r[0] * (ekf->xp[0] * ekf->xp[2] + ekf->xp[1] * ekf->xp[3]) + r[1] * (ekf->xp[2] * ekf->xp[3] - ekf->xp[0] * ekf->xp[1]) + r[2] * (0.5f - ekf->xp[1] * ekf->xp[1] - ekf->xp[2] * ekf->xp[2]))
    };
  // Innovation: v = z - h (6x1)
  float v[6];
  for (int i = 0; i < 6; i++) {
    v[i] = z[i] - h[i];
  }

  // Measurement Jacobian H (6x4)
  float H[24] = {
    -2.0f * ekf->xp[2],  2.0f * ekf->xp[3], -2.0f * ekf->xp[0],  2.0f * ekf->xp[1],
     2.0f * ekf->xp[1],  2.0f * ekf->xp[0],  2.0f * ekf->xp[3],  2.0f * ekf->xp[2],
     2.0f * ekf->xp[0], -2.0f * ekf->xp[1], -2.0f * ekf->xp[2],  2.0f * ekf->xp[3],

     2.0f * (r[0] * ekf->xp[0] + r[1] * ekf->xp[3] - r[2] * ekf->xp[2]),
     2.0f * (r[0] * ekf->xp[1] + r[1] * ekf->xp[2] + r[2] * ekf->xp[3]),
     2.0f * (-r[0] * ekf->xp[2] + r[1] * ekf->xp[1] - r[2] * ekf->xp[0]),
     2.0f * (-r[0] * ekf->xp[3] + r[1] * ekf->xp[0] + r[2] * ekf->xp[1]),

     2.0f * (-r[0] * ekf->xp[3] + r[1] * ekf->xp[0] + r[2] * ekf->xp[1]),
     2.0f * (r[0] * ekf->xp[2] - r[1] * ekf->xp[1] + r[2] * ekf->xp[0]),
     2.0f * (r[0] * ekf->xp[1] + r[1] * ekf->xp[2] + r[2] * ekf->xp[3]),
     2.0f * (-r[0] * ekf->xp[0] - r[1] * ekf->xp[3] + r[2] * ekf->xp[2]),

     2.0f * (r[0] * ekf->xp[2] - r[1] * ekf->xp[1] + r[2] * ekf->xp[0]),
     2.0f * (r[0] * ekf->xp[3] - r[1] * ekf->xp[0] - r[2] * ekf->xp[1]),
     2.0f * (r[0] * ekf->xp[0] + r[1] * ekf->xp[3] - r[2] * ekf->xp[2]), 
     2.0f * (r[0] * ekf->xp[1] + r[1] * ekf->xp[2] + r[2] * ekf->xp[3])
  };

  // Kalman Gain K = Pp * H^T * inv(H * Pp * H^T + R)
  // Compute H transpose (4x6)
  float Ht[24];
  matTranspose(H, Ht, 6, 4);
  // PHt = Pp * Ht (4x4 * 4x6 = 4x6)
  float PHt[24];
  matMultiply(ekf->Pp_data, Ht, PHt, 4, 4, 6);
  // S = H * PHt (6x4 * 4x6 = 6x6)
  float S[36];
  matMultiply(H, PHt, S, 6, 4, 6);
  // Add measurement noise: S = S + R
  for (int i = 0; i < 36; i++) {
    S[i] += ekf->R_data[i];
  }
  // Invert S (6x6) to get Si
  float Si[36];
  if (matInverse6x6(S, Si) != 0) {
    // Inversion failed (singular matrix); skip update.
    return;
  }
  // Compute Kalman gain K = PHt * Si (4x6 * 6x6 = 4x6)
  float K[24];
  matMultiply(PHt, Si, K, 4, 6, 6);

  // Update state: x = xp + K * v
  float Kv[4] = {0, 0, 0, 0};
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 6; j++) {
      Kv[i] += K[i*6 + j] * v[j];
    }
    ekf->x[i] = ekf->xp[i] + Kv[i];
  }

  // Update covariance: P = (I - K*H) * Pp
  float KH[16];
  matMultiply(K, H, KH, 4, 6, 4);
  float I_KH[16];
  for (int i = 0; i < 16; i++) {
    int row = i / 4;
    int col = i % 4;
    float I_val = (row == col) ? 1.0f : 0.0f;
    I_KH[i] = I_val - KH[i];
  }
  float P_new[16];
  matMultiply(I_KH, ekf->Pp_data, P_new, 4, 4, 4);
  for (int i = 0; i < 16; i++) {
    ekf->P_data[i] = P_new[i];
  }
  
  quat_normalize(ekf->x);
}

// Process the EKF: prediction and update steps.
// Sensor data (gyro and accelerometer) and timestep dt are provided.
void ekf_process(ekf_t *ekf, float w_x, float w_y, float w_z,
                 float a_x, float a_y, float a_z, float m_x, float m_y, float m_z, float dt) {
  // Print current quaternion state via Serial
  Serial.print(ekf->x[0]); Serial.print(" ");
  Serial.print(ekf->x[1]); Serial.print(" ");
  Serial.print(ekf->x[2]); Serial.print(" ");
  Serial.println(ekf->x[3]);
  
  // (Optional: ignore very small rotations)
  // float omega_norm = sqrt(w_x*w_x + w_y*w_y + w_z*w_z);
  // if (omega_norm < 0.005f) { return; }
  
  ekf_predict(ekf, w_x, w_y, w_z, dt);
  ekf_update(ekf, a_x, a_y, a_z, m_x, m_y, m_z);
}

// Compute a quaternion that aligns the global frame {Xg, Yg, Zg}
// to the sensor frame {Xs, Ys, Zs}, based on one reading of accel & mag.
void computeInitialOrientation(float ax, float ay, float az,
                               float mx, float my, float mz,
                               float q_out[4])
{
  // 1) Normalize accel and mag
  float a_norm = sqrtf(ax*ax + ay*ay + az*az);
  if (a_norm < 1e-6f) return; // bad data
  ax /= a_norm; ay /= a_norm; az /= a_norm;

  float m_norm = sqrtf(mx*mx + my*my + mz*mz);
  if (m_norm < 1e-6f) return; // bad data
  mx /= m_norm; my /= m_norm; mz /= m_norm;

  // 2) Usually we define 'down' = -accel if the sensor reads +g upward.
  //    If your sensor sees +9.81 for Z-up, then 'down' = -accel.
  float dx = -ax, dy = -ay, dz = -az;

  // 3) Compute 'east' = mag x down
  float ex = my*dz - mz*dy;
  float ey = mz*dx - mx*dz;
  float ez = mx*dy - my*dx;
  float e_norm = sqrtf(ex*ex + ey*ey + ez*ez);
  if (e_norm < 1e-6f) return; // collinear vectors?
  ex /= e_norm; ey /= e_norm; ez /= e_norm;

  // 4) Compute 'north' = down x east
  float nx = dy*ez - dz*ey;
  float ny = dz*ex - dx*ez;
  float nz = dx*ey - dy*ex;
  // (dx,dy,dz) is 'down'

  // 5) Build rotation matrix R_s->g (sensor to global)
  //    If you prefer global->sensor, just transpose it.
  //    Each column is a global axis in sensor coords or vice versa.
  float R[9];
  R[0] = ex; R[1] = nx; R[2] = dx;
  R[3] = ey; R[4] = ny; R[5] = dy;
  R[6] = ez; R[7] = nz; R[8] = dz;

  // 6) Convert that 3x3 rotation matrix to a quaternion
  //    This formula returns a quaternion that rotates from sensor frame to global frame.
  float trace = R[0] + R[4] + R[8];
  float qw, qx, qy, qz;
  if (trace > 0.0f) {
    float s = 0.5f / sqrtf(trace + 1.0f);
    qw = 0.25f / s;
    qx = (R[7] - R[5]) * s;
    qy = (R[2] - R[6]) * s;
    qz = (R[3] - R[1]) * s;
  } else {
    // fallback for negative trace
    if (R[0] > R[4] && R[0] > R[8]) {
      float s = 2.0f * sqrtf(1.0f + R[0] - R[4] - R[8]);
      qw = (R[7] - R[5]) / s;
      qx = 0.25f * s;
      qy = (R[1] + R[3]) / s;
      qz = (R[2] + R[6]) / s;
    } else if (R[4] > R[8]) {
      float s = 2.0f * sqrtf(1.0f + R[4] - R[0] - R[8]);
      qw = (R[2] - R[6]) / s;
      qx = (R[1] + R[3]) / s;
      qy = 0.25f * s;
      qz = (R[5] + R[7]) / s;
    } else {
      float s = 2.0f * sqrtf(1.0f + R[8] - R[0] - R[4]);
      qw = (R[3] - R[1]) / s;
      qx = (R[2] + R[6]) / s;
      qy = (R[5] + R[7]) / s;
      qz = 0.25f * s;
    }
  }

  // 7) Normalize quaternion
  float normQ = sqrtf(qw*qw + qx*qx + qy*qy + qz*qz);
  q_out[0] = qw / normQ;
  q_out[1] = qx / normQ;
  q_out[2] = qy / normQ;
  q_out[3] = qz / normQ;
}

void setup() {
  Wire.begin(33,32);
  Serial.begin(115200);
  while (!Serial) {
    delay(10);
  }

  if (!icm.begin_I2C()) {
    Serial.println("Failed to find ICM20948");
    while (1) {
      delay(10);
    }
  }
  Serial.println("ICM20948 Found!");
  
  icm.setAccelRange(ICM20948_ACCEL_RANGE_4_G);
  Serial.print("Accelerometer range set to: ");
  switch (icm.getAccelRange()) {
    case ICM20948_ACCEL_RANGE_2_G:
      Serial.println("+-2G");
      break;
    case ICM20948_ACCEL_RANGE_4_G:
      Serial.println("+-4G");
      break;
    case ICM20948_ACCEL_RANGE_8_G:
      Serial.println("+-8G");
      break;
    case ICM20948_ACCEL_RANGE_16_G:
      Serial.println("+-16G");
      break;
  }

  icm.setGyroRange(ICM20948_GYRO_RANGE_500_DPS);
  Serial.print("Gyro range set to: ");
  switch (icm.getGyroRange()) {
    case ICM20948_GYRO_RANGE_250_DPS:
      Serial.println("250 degrees/s");
      break;
    case ICM20948_GYRO_RANGE_500_DPS:
      Serial.println("500 degrees/s");
      break;
    case ICM20948_GYRO_RANGE_1000_DPS:
      Serial.println("1000 degrees/s");
      break;
    case ICM20948_GYRO_RANGE_2000_DPS:
      Serial.println("2000 degrees/s");
      break;
  }

  uint16_t accel_divisor = icm.getAccelRateDivisor();
  float accel_rate = 1125 / (1.0 + accel_divisor);
  Serial.print("Accelerometer data rate divisor set to: ");
  Serial.println(accel_divisor);
  Serial.print("Accelerometer data rate (Hz) is approximately: ");
  Serial.println(accel_rate);

  uint8_t gyro_divisor = icm.getGyroRateDivisor();
  float gyro_rate = 1100 / (1.0 + gyro_divisor);
  Serial.print("Gyro data rate divisor set to: ");
  Serial.println(gyro_divisor);
  Serial.print("Gyro data rate (Hz) is approximately: ");
  Serial.println(gyro_rate);

  Serial.print("Magnetometer data rate set to: ");
  switch (icm.getMagDataRate()) {
    case AK09916_MAG_DATARATE_SHUTDOWN:
      Serial.println("Shutdown");
      break;
    case AK09916_MAG_DATARATE_SINGLE:
      Serial.println("Single/One shot");
      break;
    case AK09916_MAG_DATARATE_10_HZ:
      Serial.println("10 Hz");
      break;
    case AK09916_MAG_DATARATE_20_HZ:
      Serial.println("20 Hz");
      break;
    case AK09916_MAG_DATARATE_50_HZ:
      Serial.println("50 Hz");
      break;
    case AK09916_MAG_DATARATE_100_HZ:
      Serial.println("100 Hz");
      break;
  }
  Serial.println();
    ekf_init(&ekf,
           /* gyro std dev */ 0.015f, 0.015f, 0.015f,
           /* accel std dev */ 0.01f, 0.01f, 0.01f,
           /* mag std dev */   0.05f,    0.05f,    0.05f);

  // 2) Grab one reading
  sensors_event_t accel, gyro, temp, mag;
  icm.getEvent(&accel, &gyro, &temp, &mag);

  float ax = accel.acceleration.x;
  float ay = accel.acceleration.y;
  float az = accel.acceleration.z;
  float mx = mag.magnetic.x;
  float my = mag.magnetic.y;
  float mz = mag.magnetic.z;

  // 3) Compute orientation from that reading
  float q_init[4];
  computeInitialOrientation(ax, ay, az, mx, my, mz, q_init);

  // 4) Override the ekf->x[] with that “start = identity” quaternion
  ekf.x[0] = q_init[0];
  ekf.x[1] = q_init[1];
  ekf.x[2] = q_init[2];
  ekf.x[3] = q_init[3];
}

void loop() {
  sensors_event_t accel;
  sensors_event_t gyro;
  sensors_event_t mag;
  sensors_event_t temp;
  icm.getEvent(&accel, &gyro, &temp, &mag);

    a[0] = accel.acceleration.x;
    a[1] = accel.acceleration.y;
    a[2] = accel.acceleration.z;
    g[0] = gyro.gyro.x;
    g[1] = gyro.gyro.y;
    g[2] = gyro.gyro.z;
    m[0] = mag.magnetic.x;
    m[1] = mag.magnetic.y;
    m[2] = mag.magnetic.z;

    t1 = t2;
    t2 = micros();
    float t_f = (t2 - t1) / 1000000.0;

    cross(m, a, east);
    cross(a, east, north);
    unit(east);
    unit(north);

    ekf_process(&ekf, g[0], g[1], g[2], a[0], a[1], a[2], m[0], m[1], m[2], t_f);

    //print_data(north);

    Serial.print(ekf.x[0], 4); Serial.print(",");
    Serial.print(ekf.x[1], 4); Serial.print(",");
    Serial.print(ekf.x[2], 4); Serial.print(",");
    Serial.println(ekf.x[3], 4);

}
