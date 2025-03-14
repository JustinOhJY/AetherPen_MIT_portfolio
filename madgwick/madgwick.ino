#include <Adafruit_ICM20X.h>
#include <Adafruit_ICM20948.h>
#include <Adafruit_Sensor.h>
#include <Wire.h>
#include <math.h>
#include <quaternions.h>
#include <linalg_funcs.h>
#include <other_funcs.h>

// Adjust the filter gain constant as needed
#define BETA 0.05f

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

// offset corrected sensor values
float a[3] = {0};
float g[3] = {0};
float m[3] = {0};

// offset values (for accelerometer and gyro)
float ao[3] = {0};
float go[3] = {0};
float mo[3] = {0};

// Global quaternion variables
struct quaternion q_r;
struct quaternion Q;

// ---- Magnetometer Calibration Globals ----
const int mag_limit = 500;  // number of samples for mag calibration
int mag_counter = 0;
float mag_calib_min[3] = { 1e6, 1e6, 1e6 };
float mag_calib_max[3] = { -1e6, -1e6, -1e6 };
float mag_offset[3] = { 0, 0, 0 };
float mag_scale[3] = { 1, 1, 1 };
bool mag_calibrated = false;

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

//---------------------------------------------------------------------
// Madgwick filter implementation
//---------------------------------------------------------------------
void madgwick_filter(float accel[3], float gyro[3],
                     float *q_a, float *q_b, float *q_c, float *q_d,
                     float dt_s) {
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

  // Precompute repeated terms
  _2q0 = 2.0f * (*q_a);
  _2q1 = 2.0f * (*q_b);
  _2q2 = 2.0f * (*q_c);
  _2q3 = 2.0f * (*q_d);
  q0q0 = (*q_a) * (*q_a);
  q1q1 = (*q_b) * (*q_b);
  q2q2 = (*q_c) * (*q_c);
  q3q3 = (*q_d) * (*q_d);

  // Compute the objective function (error between estimated and measured gravity)
  float f0 = _2q1 * (*q_d) - _2q0 * (*q_c) - ax;
  float f1 = _2q0 * (*q_b) + _2q2 * (*q_d) - ay;
  float f2 = q0q0 - q1q1 - q2q2 + q3q3 - az;

  // Compute the gradient (Jacobian^T * f)
  s0 = -_2q2 * f0 + _2q1 * f1;
  s1 =  _2q3 * f0 + _2q0 * f1 - 4.0f * (*q_b) * f2;
  s2 = -_2q0 * f0 + _2q3 * f1 - 4.0f * (*q_c) * f2;
  s3 =  _2q1 * f0 + _2q2 * f1;

  // Normalize the gradient
  recipNorm = sqrtf(s0 * s0 + s1 * s1 + s2 * s2 + s3 * s3);
  if (recipNorm < 1e-8f)
    return;
  recipNorm = 1.0f / recipNorm;
  s0 *= recipNorm;
  s1 *= recipNorm;
  s2 *= recipNorm;
  s3 *= recipNorm;

  // Compute quaternion derivative from gyroscope measurements
  float gx = gyro[0], gy = gyro[1], gz = gyro[2];
  qDot0 = 0.5f * (-(*q_b) * gx - (*q_c) * gy - (*q_d) * gz) - BETA * s0;
  qDot1 = 0.5f * ( (*q_a) * gx + (*q_c) * gz - (*q_d) * gy) - BETA * s1;
  qDot2 = 0.5f * ( (*q_a) * gy - (*q_b) * gz + (*q_d) * gx) - BETA * s2;
  qDot3 = 0.5f * ( (*q_a) * gz + (*q_b) * gy - (*q_c) * gx) - BETA * s3;

  // Integrate to yield new quaternion values
  *q_a += qDot0 * dt_s;
  *q_b += qDot1 * dt_s;
  *q_c += qDot2 * dt_s;
  *q_d += qDot3 * dt_s;

  // Normalize the updated quaternion
  recipNorm = sqrtf((*q_a) * (*q_a) + (*q_b) * (*q_b) +
                    (*q_c) * (*q_c) + (*q_d) * (*q_d));
  if (recipNorm < 1e-8f)
    return;
  recipNorm = 1.0f / recipNorm;
  *q_a *= recipNorm;
  *q_b *= recipNorm;
  *q_c *= recipNorm;
  *q_d *= recipNorm;
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

  // Initialize global quaternion Q to the identity quaternion
  Q.qw = 1.0f;
  Q.qx = 0.0f;
  Q.qy = 0.0f;
  Q.qz = 0.0f;

  // Initialize timing variables
  t2 = micros();
}

void loop() {
  sensors_event_t accel;
  sensors_event_t gyro;
  sensors_event_t mag;
  sensors_event_t temp;
  icm.getEvent(&accel, &gyro, &temp, &mag);

  // Get accelerometer and gyro values (apply gyro offsets)
  a[0] = accel.acceleration.x;
  a[1] = accel.acceleration.y;
  a[2] = accel.acceleration.z;  // Invert if needed for your setup
  g[0] = gyro.gyro.x;
  g[1] = gyro.gyro.y;
  g[2] = gyro.gyro.z;


  // Compute the elapsed time in seconds for this loop iteration
  t1 = t2;
  t2 = micros();
  float dt = (t2 - t1) / 1000000.0f;

  // Update the orientation using the Madgwick filter
  madgwick_filter(a, g, &Q.qw, &Q.qx, &Q.qy, &Q.qz, dt);

  // (Optional) Update east and north vectors if needed
  cross(m, a, east);
  cross(a, east, north);
  unit(east);
  unit(north);

  // Print the updated quaternion (orientation)
  Serial.print(Q.qw, 4); Serial.print(", ");
  Serial.print(Q.qx, 4); Serial.print(", ");
  Serial.print(Q.qy, 4); Serial.print(", ");
  Serial.println(Q.qz, 4);
}
