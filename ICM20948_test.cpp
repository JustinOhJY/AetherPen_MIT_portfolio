#include <Adafruit_ICM20X.h>
#include <Adafruit_ICM20948.h>
#include <Adafruit_Sensor.h>
#include <Wire.h>
#include <math.h>

Adafruit_ICM20948 icm;

const float gravity = 9.80665;

unsigned long t1 = 0, t2 = 0;
int counter = 0;

const int limit = 10;

float east[3] = {0};
float north[3] = {0};

bool stationary = true;

//vector pointing down
float down[3] = {0};

//offset corrected values
float a[3] = {0};
float g[3] = {0};
float m[3] = {0};

//offset values
float ao[3] = {0};
float go[3] = {0};
float mo[3] = {0};

float get_norm(const float (&arr)[3]) {
  float sum = 0;
  for (int i = 0; i < 3; i++) {
    sum += arr[i] * arr[i];
  }
  return sqrt(sum);
}

void unit(float (&arr)[3]) {
  float norm = get_norm(arr);
  if (norm >= 1e-6) {
    for (int i = 0; i < 3; i++) {
      arr[i] /= norm;
    }
  }
}


void cross(const float (&arr1)[3], const float (&arr2)[3], float (&arr3)[3]) {
  arr3[0] = arr1[1] * arr2[2] - arr1[2] * arr2[1];
  arr3[1] = arr2[0] * arr1[2] - arr2[2] * arr1[0];
  arr3[2] = arr1[0] * arr2[1] - arr1[1] * arr2[0];
}

struct quaternion {
  float qw;
  float qx;
  float qy;
  float qz;
  float dt;
  float q_norm;
  float angle;
  // default constructor
  quaternion() {
    qw = 1.0;
    qx = 0.0;
    qy = 0.0;
    qz = 0.0;
    angle = 0;
  }

  //parametrized constructor
  quaternion(float gx, float gy, float gz, float time) {
    qx = gx;
    qy = gy;
    qz = gz;
    dt = time;
    float gyro_norm = sqrt(pow(gx, 2) + pow(gy, 2) + pow(gz, 2));
    q_norm = gyro_norm;
    angle = (gyro_norm * dt)/2.0;
    qw = cos(angle);
    normalize();
  }

  //copy constructor
  quaternion(const quaternion &q) {
    qw = q.qw;
    qx = q.qx;
    qy = q.qy;
    qz = q.qz;
    angle = q.angle;
  }
  
  void normalize() {
    float factor = sin(angle);
    float norm = sqrt(pow(qx, 2) + pow(qy, 2) + pow(qz, 2));
    if (norm >= 1e-6) {
      qx /= norm;
      qy /= norm;
      qz /= norm;
      qx *= factor;
      qy *= factor;
      qz *= factor;
    }
  }

  void rotate(struct quaternion qr) {
    quaternion result;
    result.qw = qw;
    result.qx = qx * (1 - 2 * qr.qy * qr.qy - 2 * qr.qz * qr.qz) + qy * (2 * qr.qx * qr.qy - 2 * qr.qw * qr.qz) + qz * (2 * qr.qx * qr.qz + 2 * qr.qw * qr.qy);
    result.qy = qx * (2 * qr.qx * qr.qy + 2 * qr.qw * qr.qz) + qy * (1 - 2 * qr.qx * qr.qx - 2 * qr.qz * qr.qz) + qz * (2 * qr.qy * qr.qz - 2 * qr.qw * qr.qx);
    result.qz = qx * (2 * qr.qx * qr.qz - 2 * qr.qw * qr.qy) + qy * (2 * qr.qy * qr.qz + 2 * qr.qw * qr.qx) + qz * (1 - 2 * qr.qx * qr.qx - 2 * qr.qy * qr.qy);
    *this = result;
  }

  quaternion operator*(const quaternion &q1) {
    quaternion q_output;
    q_output.qw = this->qw * q1.qw - this->qx * q1.qx - this->qy * q1.qy - this->qz * q1.qz;
    q_output.qx = this->qw * q1.qx + this->qx * q1.qw + this->qy * q1.qz - this->qz * q1.qy;
    q_output.qy = this->qw * q1.qy - this->qx * q1.qz + this->qy * q1.qw + this->qz * q1.qx;
    q_output.qz = this->qw * q1.qz + this->qx * q1.qy - this->qy * q1.qx + this->qz * q1.qw;
    return q_output;
  }
};

float get_norm() {
  
}

struct quaternion rotation_quaternion;

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

  //  icm.setAccelRateDivisor(4095);
  uint16_t accel_divisor = icm.getAccelRateDivisor();
  float accel_rate = 1125 / (1.0 + accel_divisor);

  Serial.print("Accelerometer data rate divisor set to: ");
  Serial.println(accel_divisor);
  Serial.print("Accelerometer data rate (Hz) is approximately: ");
  Serial.println(accel_rate);

  //  icm.setGyroRateDivisor(255);
  uint8_t gyro_divisor = icm.getGyroRateDivisor();
  float gyro_rate = 1100 / (1.0 + gyro_divisor);

  Serial.print("Gyro data rate divisor set to: ");
  Serial.println(gyro_divisor);
  Serial.print("Gyro data rate (Hz) is approximately: ");
  Serial.println(gyro_rate);

  // icm.setMagDataRate(AK09916_MAG_DATARATE_10_HZ);
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

  if (counter < limit) {
    ao[0] += accel.acceleration.x;
    ao[1] += accel.acceleration.y;
    ao[2] += accel.acceleration.z;
    go[0] += gyro.gyro.x;
    go[1] += gyro.gyro.y;
    go[2] += gyro.gyro.z;
    mo[0] += mag.magnetic.x;
    mo[1] += mag.magnetic.y;
    mo[2] += mag.magnetic.z;

    if (counter == limit - 1) {
      for (int i = 0; i < 3; i++) {
        ao[i] /= limit;
        go[i] /= limit;
        mo[i] /= limit;
      }
    }
    counter += 1;
  }
  else {
    a[0] = accel.acceleration.x - ao[0];
    a[1] = accel.acceleration.y - ao[1];
    a[2] = accel.acceleration.z - ao[2] + gravity;
    g[0] = gyro.gyro.x - go[0];
    g[1] = gyro.gyro.y - go[1];
    g[2] = gyro.gyro.z - go[2];
    m[0] = mag.magnetic.x;
    m[1] = mag.magnetic.y;
    m[2] = mag.magnetic.z;

    t1 = t2;
    t2 = micros();
    float t_f = (t2 - t1) / 1000000.0;

    if(abs(get_norm(a) - gravity) <= 10 - gravity) {
      cross(a, m, east);
      cross(east, a, north);
      unit(east);
      unit(north);
    }

    quaternion q(g[0], g[1], g[2], t_f);
    quaternion q_final = q;
    q_final.rotate(rotation_quaternion);
    rotation_quaternion = rotation_quaternion * q;

    Serial.print(q_final.qw, 4); Serial.print(",");
    Serial.print(q_final.qx, 4); Serial.print(",");
    Serial.print(q_final.qy, 4); Serial.print(",");
    Serial.println(q_final.qz, 4);

    // // prints accel data
    // for (int i = 0; i < 3; i++) {
    //   Serial.print(a[i]); Serial.print(",");
    // }
    // Serial.println();

    //  //prints gyro data
    // for (int i = 0; i < 3; i++) {
    //   Serial.print(g[i]); Serial.print(",");
    // }
    // Serial.println();

    //  //prints mag data
    // for (int i = 0; i < 3; i++) {
    //   Serial.print(m[i]); Serial.print(",");
    // }
    // Serial.println();

    // //prints north
    // for (int i = 0; i < 3; i++) {
    //   Serial.print(north[i], 6); Serial.print(",");
    // }
    // Serial.println();

    // //prints east
    // for (int i = 0; i < 3; i++) {
    //   Serial.print(east[i]); Serial.print(",");
    // }
    // Serial.println();
  }
}
