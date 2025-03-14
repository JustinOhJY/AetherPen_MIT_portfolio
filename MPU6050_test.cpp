#include "I2Cdev.h"
#include "MPU6050.h"
#include "math.h"

#if I2CDEV_IMPLEMENTATION == I2CDEV_ARDUINO_WIRE
#include "Wire.h"
#endif

MPU6050 accelgyro;

int16_t ax, ay, az;
int16_t gx, gy, gz;
float ax_offset = 0.0, ay_offset = 0.0, az_offset = -2048.0, gx_offset = 0.0, gy_offset = 0.0, gz_offset = 0.0;
int counter = -200;
const int limit = 25;
const float gravity = 9.80665;
float offset_normalizer_accel[3][limit];
float offset_normalizer_gyro[3][limit];
bool is_offset = false;
unsigned long t1;
unsigned long t2;

#define LED_PIN 13
bool blinkState = false;

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
    angle = (gyro_norm * dt * M_PI) / 360.0;
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
    if (norm > 0.001) {
      qx /= norm;
      qy /= norm;
      qz /= norm;
      qx *= factor;
      qy *= factor;
      qz *= factor;
    } 
    else {
      qw = 1.0;
      qx = 0.0;
      qy = 0.0;
      qz = 0.0;
    }
  }

  void rotate(struct quaternion qr) {
    struct quaternion result;
    result.qw = qw;
    result.qx = qx * (1 - 2 * qr.qy * qr.qy - 2 * qr.qz * qr.qz) + qy * (2 * qr.qx * qr.qy - 2 * qr.qw * qr.qz) + qz * (2 * qr.qx * qr.qz + 2 * qr.qw * qr.qy);
    result.qy = qx * (2 * qr.qx * qr.qy + 2 * qr.qw * qr.qz) + qy * (1 - 2 * qr.qx * qr.qx - 2 * qr.qz * qr.qz) + qz * (2 * qr.qy * qr.qz - 2 * qr.qw * qr.qx);
    result.qz = qx * (2 * qr.qx * qr.qz - 2 * qr.qw * qr.qy) + qy * (2 * qr.qy * qr.qz + 2 * qr.qw * qr.qx) + qz * (1 - 2 * qr.qx * qr.qx - 2 * qr.qy * qr.qy);
    *this = result;
  }

  quaternion operator*(const quaternion &q1) {
     struct quaternion q_output;
    q_output.qw = this->qw * q1.qw - this->qx * q1.qx - this->qy * q1.qy - this->qz * q1.qz;
    q_output.qx = this->qw * q1.qx + this->qx * q1.qw + this->qy * q1.qz - this->qz * q1.qy;
    q_output.qy = this->qw * q1.qy - this->qx * q1.qz + this->qy * q1.qw + this->qz * q1.qx;
    q_output.qz = this->qw * q1.qz + this->qx * q1.qy - this->qy * q1.qx + this->qz * q1.qw;
    return q_output;
  }
};

struct quaternion rotation_quaternion;

void setup() {
// Join I2C bus
#if I2CDEV_IMPLEMENTATION == I2CDEV_ARDUINO_WIRE
  Wire.begin();
#elif I2CDEV_IMPLEMENTATION == I2CDEV_BUILTIN_FASTWIRE
  Fastwire::setup(400, true);
#endif

  // Initialize serial communication
  Serial.begin(115200);

  // Reset the MPU6050
  accelgyro.reset();
  delay(100);

  // Initialize MPU6050 and set the range
  Serial.println("Initializing I2C devices...");
  accelgyro.initialize();
  accelgyro.setFullScaleAccelRange(MPU6050_ACCEL_FS_16);  // Set range to ±16g
  accelgyro.setFullScaleGyroRange(MPU6050_GYRO_FS_2000);

  // Verify connection
  Serial.println("Testing device connections...");
  Serial.println(accelgyro.testConnection() ? "MPU6050 connection successful" : "MPU6050 connection failed");

  // Directly read the accel config register (0x1C) to verify the initial range setting
  uint8_t accelConfig = 0;
  I2Cdev::readByte(0x68, 0x1C, &accelConfig);  // Read from register 0x1C
  Serial.print("Accel Config Register (0x1C) initial setting: ");
  Serial.println(accelConfig, BIN);  // Output in binary to inspect individual bits

  // Force the range setting to ±16g by writing 0x18 (binary: 00011000) to register 0x1C
  I2Cdev::writeByte(0x68, 0x1C, 0x18);  // Write ±16g configuration to register 0x1C
  delay(10);                            // Short delay to ensure the write completes

  // Read the register again to confirm the updated setting
  I2Cdev::readByte(0x68, 0x1C, &accelConfig);
  Serial.print("Accel Config Register (0x1C) after forcing ±16g: ");
  Serial.println(accelConfig, BIN);

  // Initialize the array values
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < limit - 1; j++) {
      offset_normalizer_accel[i][j] = 0;
      offset_normalizer_gyro[i][j] = 0;
    }
  }
  // Configure LED pin for output
  pinMode(LED_PIN, OUTPUT);
  accelgyro.setDLPFMode(6);  // Use 6 for highest filtering
}

void loop() {
  t1 = t2;
  accelgyro.getMotion6(&ax, &ay, &az, &gx, &gy, &gz);
  if (counter != 0) {
    t2 = micros();
  }

  if (counter < limit) {
    if (counter >= 0) {
      offset_normalizer_accel[0][counter] = ax;
      offset_normalizer_accel[1][counter] = ay;
      offset_normalizer_accel[2][counter] = az;
      offset_normalizer_gyro[0][counter] = gx;
      offset_normalizer_gyro[1][counter] = gy;
      offset_normalizer_gyro[2][counter] = gz;
    }
    if (counter == limit - 1) {
      for (int i = 0; i < limit; i++) {
        ax_offset += offset_normalizer_accel[0][i] / (1.0 * limit);
        ay_offset += offset_normalizer_accel[1][i] / (1.0 * limit);
        az_offset += offset_normalizer_accel[2][i] / (1.0 * limit);
        gx_offset += offset_normalizer_gyro[0][i] / (1.0 * limit);
        gy_offset += offset_normalizer_gyro[1][i] / (1.0 * limit);
        gz_offset += offset_normalizer_gyro[2][i] / (1.0 * limit);
      }
    }
    counter++;
  }

  else {
    float ax_scaled = (ax - ax_offset) * gravity / 2048.0;
    float ay_scaled = (ay - ay_offset) * gravity / 2048.0;
    float az_scaled = (az - az_offset) * gravity / 2048.0;
    float gx_scaled = (gx - gx_offset) / (131.0 / 8);
    float gy_scaled = (gy - gy_offset) / (131.0 / 8);
    float gz_scaled = (gz - gz_offset) / (131.0 / 8);
    float t_f = (t2 - t1) / 1000000.0;

    struct quaternion q(gx_scaled, gy_scaled, gz_scaled, t_f);
    struct quaternion q_final = q;
    q_final.rotate(rotation_quaternion);
    rotation_quaternion = rotation_quaternion * q;

    //gx pos is down, gz pos is right, gy pos is right.

    // Serial.print(ax_scaled, 4); Serial.print(",");
    // Serial.print(ay_scaled, 4); Serial.print(",");
    // Serial.print(az_scaled, 4); Serial.print(",");
    // Serial.print(gx_scaled, 4); Serial.print(",");
    // Serial.print(gy_scaled, 4); Serial.print(",");
    // Serial.println(gz_scaled, 4);
    // Serial.flush();

    // if (q.qw != 1) {
      // Serial.print(q.qw, 4); Serial.print(",");
      // Serial.print(q.qx, 4); Serial.print(",");
      // Serial.print(q.qy, 4); Serial.print(",");
      // Serial.println(q.qz, 4);
      // Serial.flush();
      // Serial.print(rotation_quaternion.qw, 4);
      // Serial.print(",");
      // Serial.print(rotation_quaternion.qx, 4);
      // Serial.print(",");
      // Serial.print(rotation_quaternion.qy, 4);
      // Serial.print(",");
      // Serial.println(rotation_quaternion.qz, 4);
      // Serial.flush();
    // }

    Serial.print(q_final.qw, 6); Serial.print(",");
    Serial.print(q_final.qx, 6); Serial.print(",");
    Serial.print(q_final.qy, 6); Serial.print(",");
    Serial.println(q_final.qz, 6);
    Serial.flush();

    // Serial.print(rotation_quaternion.qw, 4); Serial.print(",");
    // Serial.print(rotation_quaternion.qx, 4); Serial.print(",");
    // Serial.print(rotation_quaternion.qy, 4); Serial.print(",");
    // Serial.println(rotation_quaternion.qz, 4);
    // Serial.flush();

    // blink LED to indicate activity
    blinkState = !blinkState;
    digitalWrite(LED_PIN, blinkState);
  }
}
