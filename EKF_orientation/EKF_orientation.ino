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
  float P_data[16];    // Covariance matrix (4x4)
  float Pp_data[16];   // Predicted covariance (4x4)
  float sigma_w_data[9]; // Gyro noise covariance (3x3)
  float R_data[9];     // Measurement noise covariance (3x3)
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
              float sigma_ax, float sigma_ay, float sigma_az) {
  // Set initial quaternion state (no rotation)
  ekf->x[0] = 1.0f;
  ekf->x[1] = 0.0f;
  ekf->x[2] = 0.0f;
  ekf->x[3] = 0.0f;

  // Initialize covariance P_data to identity (4x4)
  for (int i = 0; i < 16; i++) {
    ekf->P_data[i] = 0.0f;
  }
  ekf->P_data[0]  = 1.0f;
  ekf->P_data[5]  = 1.0f;
  ekf->P_data[10] = 1.0f;
  ekf->P_data[15] = 1.0f;
  
  // Set gyro spectral noise covariance matrix (3x3)
  ekf->sigma_w_data[0] = sigma_wx * sigma_wx;
  ekf->sigma_w_data[1] = 0.0f;
  ekf->sigma_w_data[2] = 0.0f;
  ekf->sigma_w_data[3] = 0.0f;
  ekf->sigma_w_data[4] = sigma_wy * sigma_wy;
  ekf->sigma_w_data[5] = 0.0f;
  ekf->sigma_w_data[6] = 0.0f;
  ekf->sigma_w_data[7] = 0.0f;
  ekf->sigma_w_data[8] = sigma_wz * sigma_wz;
  
  // Set measurement noise covariance matrix (3x3)
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

// Update step: correct the predicted state with accelerometer data.
void ekf_update(ekf_t *ekf, float a_x, float a_y, float a_z) {
  // Normalize accelerometer measurement (assumed to measure gravity)
  float a_norm = sqrt(a_x*a_x + a_y*a_y + a_z*a_z);
  if (a_norm > 0.0f) {
    a_x /= a_norm;
    a_y /= a_norm;
    a_z /= a_norm;
  }
  // Measurement vector z (3x1)
  float z[3] = { a_x, a_y, a_z };

  // Construct h(xp) from predicted quaternion (xp)
  float h[3] = {
    2.0f * (ekf->xp[1] * ekf->xp[3] - ekf->xp[0] * ekf->xp[2]),
    2.0f * (ekf->xp[0] * ekf->xp[1] + ekf->xp[2] * ekf->xp[3]),
    ekf->xp[0]*ekf->xp[0] - ekf->xp[1]*ekf->xp[1] - ekf->xp[2]*ekf->xp[2] + ekf->xp[3]*ekf->xp[3]
  };
  // Innovation: v = z - h (3x1)
  float v[3];
  for (int i = 0; i < 3; i++) {
    v[i] = z[i] - h[i];
  }

  // Measurement Jacobian H (3x4)
  float H[12] = {
    -2.0f * ekf->xp[2],  2.0f * ekf->xp[3], -2.0f * ekf->xp[0],  2.0f * ekf->xp[1],
     2.0f * ekf->xp[1],  2.0f * ekf->xp[0],  2.0f * ekf->xp[3],  2.0f * ekf->xp[2],
     2.0f * ekf->xp[0], -2.0f * ekf->xp[1], -2.0f * ekf->xp[2],  2.0f * ekf->xp[3]
  };

  // Kalman Gain K = Pp * H^T * inv(H * Pp * H^T + R)
  // Compute H transpose (4x3)
  float Ht[12];
  matTranspose(H, Ht, 3, 4);
  // PHt = Pp * Ht (4x4 * 4x3 = 4x3)
  float PHt[12];
  matMultiply(ekf->Pp_data, Ht, PHt, 4, 4, 3);
  // S = H * PHt (3x4 * 4x3 = 3x3)
  float S[9];
  matMultiply(H, PHt, S, 3, 4, 3);
  // Add measurement noise: S = S + R
  for (int i = 0; i < 9; i++) {
    S[i] += ekf->R_data[i];
  }
  // Invert S (3x3) to get Si
  float Si[9];
  if (matInverse3x3(S, Si) != 0) {
    // Inversion failed (singular matrix); skip update.
    return;
  }
  // Compute Kalman gain K = PHt * Si (4x3 * 3x3 = 4x3)
  float K[12];
  matMultiply(PHt, Si, K, 4, 3, 3);

  // Update state: x = xp + K * v
  float Kv[4] = {0, 0, 0, 0};
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 3; j++) {
      Kv[i] += K[i*3 + j] * v[j];
    }
    ekf->x[i] = ekf->xp[i] + Kv[i];
  }

  // Update covariance: P = (I - K*H) * Pp
  float KH[16];
  matMultiply(K, H, KH, 4, 3, 4);
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
                 float a_x, float a_y, float a_z, float dt) {
  // Print current quaternion state via Serial
  Serial.print(ekf->x[0]); Serial.print(" ");
  Serial.print(ekf->x[1]); Serial.print(" ");
  Serial.print(ekf->x[2]); Serial.print(" ");
  Serial.println(ekf->x[3]);
  
  // (Optional: ignore very small rotations)
  // float omega_norm = sqrt(w_x*w_x + w_y*w_y + w_z*w_z);
  // if (omega_norm < 0.005f) { return; }
  
  ekf_predict(ekf, w_x, w_y, w_z, dt);
  ekf_update(ekf, a_x, a_y, a_z);
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
  ekf_init(&ekf, 0.015f, 0.015f, 0.015f, 0.23f, 0.23f, 0.23f);
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
    m[0] = (mag.magnetic.x - 25.200001) * 0.724638;
    m[1] = (mag.magnetic.y - 12.674999) * 1.010101;
    m[2] = (mag.magnetic.z - 90.975006) * 1.587301;

    t1 = t2;
    t2 = micros();
    float t_f = (t2 - t1) / 1000000.0;

    cross(m, a, east);
    cross(a, east, north);
    unit(east);
    unit(north);

    ekf_process(&ekf, g[0], g[1], g[2], a[0], a[1], a[2], t_f);

    //print_data(north);

    Serial.print(ekf.x[0], 4); Serial.print(",");
    Serial.print(ekf.x[1], 4); Serial.print(",");
    Serial.print(ekf.x[2], 4); Serial.print(",");
    Serial.println(ekf.x[3], 4);

}
