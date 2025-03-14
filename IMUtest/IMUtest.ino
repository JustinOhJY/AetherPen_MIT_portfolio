#include <Adafruit_ICM20X.h>
#include <Adafruit_ICM20948.h>
#include <Adafruit_Sensor.h>
#include <Wire.h>
#include <math.h>
#include <quaternions.h>
#include <linalg_funcs.h>
#include <other_funcs.h>

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

// // ---- Magnetometer Calibration Globals ----
// const int mag_limit = 500;  // number of samples for mag calibration
// int mag_counter = 0;
// float mag_calib_min[3] = { 1e6, 1e6, 1e6 };
// float mag_calib_max[3] = { -1e6, -1e6, -1e6 };
// float mag_offset[3] = { 0, 0, 0 };
// float mag_scale[3] = { 1, 1, 1 };
// bool mag_calibrated = false;

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
}

void loop() {
  sensors_event_t accel;
  sensors_event_t gyro;
  sensors_event_t mag;
  sensors_event_t temp;
  icm.getEvent(&accel, &gyro, &temp, &mag);

  else {
    a[0] = accel.acceleration.x;
    a[1] = accel.acceleration.y;
    a[2] = -1.0 * accel.acceleration.z;
    g[0] = gyro.gyro.x;
    g[1] = gyro.gyro.y;
    g[2] = gyro.gyro.z;
    
    // // --- Magnetometer Calibration Section ---
    // if (!mag_calibrated) {
    //   // Accumulate calibration samples for the magnetometer:
    //   mag_counter++;
    //   if (mag.magnetic.x < mag_calib_min[0]) mag_calib_min[0] = mag.magnetic.x;
    //   if (mag.magnetic.y < mag_calib_min[1]) mag_calib_min[1] = mag.magnetic.y;
    //   if (mag.magnetic.z < mag_calib_min[2]) mag_calib_min[2] = mag.magnetic.z;
    //   if (mag.magnetic.x > mag_calib_max[0]) mag_calib_max[0] = mag.magnetic.x;
    //   if (mag.magnetic.y > mag_calib_max[1]) mag_calib_max[1] = mag.magnetic.y;
    //   if (mag.magnetic.z > mag_calib_max[2]) mag_calib_max[2] = mag.magnetic.z;
      
    //   if (mag_counter >= mag_limit) {
    //     // Compute offsets
    //     for (int i = 0; i < 3; i++) {
    //       mag_offset[i] = (mag_calib_max[i] + mag_calib_min[i]) / 2.0;
    //     }
    //     // Optionally compute scales:
    //     float scale[3];
    //     for (int i = 0; i < 3; i++) {
    //       scale[i] = (mag_calib_max[i] - mag_calib_min[i]) / 2.0;
    //     }
    //     float avg_scale = (scale[0] + scale[1] + scale[2]) / 3.0;
    //     for (int i = 0; i < 3; i++) {
    //       // scale factor to bring each axis to the average
    //       mag_scale[i] = avg_scale / scale[i];
    //     }
    //     mag_calibrated = true;
    //     Serial.println("Magnetometer calibration complete.");
    //   }
    // }
    
    // if (mag_calibrated) {
    //   m[0] = (mag.magnetic.x - mag_offset[0]) * mag_scale[0];
    //   m[1] = (mag.magnetic.y - mag_offset[1]) * mag_scale[1];
    //   m[2] = (mag.magnetic.z - mag_offset[2]) * mag_scale[2];

    //   Serial.println("calibration values");
    //   Serial.println(mag_offset[0], 6);
    //   Serial.println(mag_scale[0], 6);
    //   Serial.println(mag_offset[1], 6);
    //   Serial.println(mag_scale[1], 6);
    //   Serial.println(mag_offset[2], 6);
    //   Serial.println(mag_scale[2], 6);
    // } else {
    //   m[0] = mag.magnetic.x;
    //   m[1] = mag.magnetic.y;
    //   m[2] = mag.magnetic.z;
    // }

      m[0] = (mag.magnetic.x - 25.200001) * 0.724638;
      m[1] = (mag.magnetic.y - 12.674999) * 1.010101;
      m[2] = (mag.magnetic.z - 90.975006) * 1.587301;

    // //accel-gyro orientation calibration
    // if (is_stationary(staionary_counter, a) == true) {
    //   float accel_norm = a.get_norm();
    //   quaternion adjustment_q;
    //   adjustment_q.qw = cos(acos(a[2])/2.0);
    //   adjustment_q.qx = a[1];
    //   adjustment_q.qy = a[2];
    //   adjustment_q.qz = 0;
    //   adjustment_q.normalize();
    // }

    // start orientation calibration
    // if (adjusted == false) {
    //   quaternion q_adj;
    //   float start_theta = acos((-a[2]/get_norm(a)))/2;
    //   q_adj.qw = cos(start_theta);
    //   q_adj.qx = (a[1] / sqrt(a[0] * a[0] + a[1] * a[1])) * sin(start_theta);
    //   q_adj.qy = (-a[0] / sqrt(a[0] * a[0] + a[1] * a[1])) * sin(start_theta);
    //   q_adj.qz = 0.0;
    //   q_adj.normalize();
    //   Q = q_adj * Q;
    //   Q.normalize();
    //   adjusted = true;
    // }

    // Serial.print(adjustment_q.qw, 4); Serial.print(",");
    // Serial.print(adjustment_q.qx, 4); Serial.print(",");
    // Serial.print(adjustment_q.qy, 4); Serial.print(",");
    // Serial.println(adjustment_q.qz, 4);

    t1 = t2;
    t2 = micros();
    float t_f = (t2 - t1) / 1000000.0;

    cross(m, a, east);
    cross(a, east, north);
    unit(east);
    unit(north);

    quaternion q;

    
    
    quaternion q_global = Q * q * Q.inv();
    
    Q = q_global * Q;
    Q.normalize();

    //print_data(north);

    // Serial.print(Q.qw, 4); Serial.print(",");
    // Serial.print(Q.qx, 4); Serial.print(",");
    // Serial.print(Q.qy, 4); Serial.print(",");
    // Serial.println(Q.qz, 4);

  }
}
