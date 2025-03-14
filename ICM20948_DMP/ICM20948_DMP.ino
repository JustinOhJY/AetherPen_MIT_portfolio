/****************************************************************
 * Modified Example6_DMP_Quat9_Orientation.ino
 * ICM 20948 Arduino Library Demo
 * This example initializes the DMP and outputs quaternion data
 * as a comma-separated line: w,x,y,z
 *
 * Make sure you have enabled DMP support by uncommenting
 * #define ICM_20948_USE_DMP in ICM_20948_C.h.
 *
 * Based on original code by SparkFun / Owen Lyke.
 ***************************************************************/

//#define QUAT_ANIMATION // Comment this out so we don't output JSON

#include "ICM_20948.h" // SparkFun ICM-20948 library header

//#define USE_SPI       // Uncomment this to use SPI

#define SERIAL_PORT Serial

#define SPI_PORT SPI // Used only when USE_SPI is defined
#define CS_PIN 2     // CS pin used when SPI is defined

#define WIRE_PORT Wire // Used when SPI is not defined
// Set AD0 value per your hardware; on SparkFun boards typically 1.
#define AD0_VAL 1

#ifdef USE_SPI
ICM_20948_SPI myICM; // Use SPI
#else
ICM_20948_I2C myICM; // Use I2C
#endif

void setup()
{
  SERIAL_PORT.begin(115200); // Start serial
  while (!SERIAL_PORT); // Wait for Serial to initialize
  SERIAL_PORT.println(F("ICM-20948 DMP Quaternion CSV Output"));

  delay(100);

#ifdef USE_SPI
  SPI_PORT.begin();
#else
  WIRE_PORT.begin();
  WIRE_PORT.setClock(400000);
#endif

  bool initialized = false;
  while (!initialized)
  {
#ifdef USE_SPI
    myICM.begin(CS_PIN, SPI_PORT);
#else
    myICM.begin(WIRE_PORT, AD0_VAL);
#endif

    SERIAL_PORT.print(F("Sensor init returned: "));
    SERIAL_PORT.println(myICM.statusString());
    if (myICM.status != ICM_20948_Stat_Ok)
    {
      SERIAL_PORT.println(F("Initialization failed. Trying again..."));
      delay(500);
    }
    else
    {
      initialized = true;
    }
  }

  SERIAL_PORT.println(F("Device connected!"));

  // Initialize the DMP (load firmware and configure registers)
  bool success = (myICM.initializeDMP() == ICM_20948_Stat_Ok);
  // Enable the DMP orientation sensor (Quat9)
  success &= (myICM.enableDMPSensor(INV_ICM20948_SENSOR_ORIENTATION) == ICM_20948_Stat_Ok);

  // Set DMP output data rate as needed (here we use the maximum for Quat9)
  success &= (myICM.setDMPODRrate(DMP_ODR_Reg_Quat9, 0) == ICM_20948_Stat_Ok);

  // Enable FIFO, enable DMP, and reset both FIFO and DMP
  success &= (myICM.enableFIFO() == ICM_20948_Stat_Ok);
  success &= (myICM.enableDMP() == ICM_20948_Stat_Ok);
  success &= (myICM.resetDMP() == ICM_20948_Stat_Ok);
  success &= (myICM.resetFIFO() == ICM_20948_Stat_Ok);

  if (success)
  {
    SERIAL_PORT.println(F("DMP enabled and configured successfully."));
  }
  else
  {
    SERIAL_PORT.println(F("DMP configuration failed!"));
    SERIAL_PORT.println(F("Make sure you have uncommented #define ICM_20948_USE_DMP in ICM_20948_C.h."));
    while (1);
  }
}

void loop()
{
  // Create structure to hold DMP data from FIFO
  icm_20948_DMP_data_t data;
  myICM.readDMPdataFromFIFO(&data);

  // Check if we received a valid frame or more data is available
  if ((myICM.status == ICM_20948_Stat_Ok) || (myICM.status == ICM_20948_Stat_FIFOMoreDataAvail))
  {
    // Check that the header indicates Quat9 data is present
    if ((data.header & DMP_header_bitmap_Quat9) > 0)
    {
      // The DMP outputs a 9-axis quaternion where only Q1, Q2, Q3 are directly given.
      // We compute Q0 as: Q0 = sqrt( 1 - (Q1^2 + Q2^2 + Q3^2) )
      // The raw quaternion data is scaled by 2^30.
      double q1 = ((double)data.Quat9.Data.Q1) / 1073741824.0;
      double q2 = ((double)data.Quat9.Data.Q2) / 1073741824.0;
      double q3 = ((double)data.Quat9.Data.Q3) / 1073741824.0;
      double q0 = sqrt(1.0 - (q1*q1 + q2*q2 + q3*q3));
      
      // Output CSV: w,x,y,z
      Serial.print(q0, 3);
      Serial.print(",");
      Serial.print(q1, 3);
      Serial.print(",");
      Serial.print(q2, 3);
      Serial.print(",");
      Serial.print(q3, 3);
    }
  }
  // If no new data available, delay a little (but if more data remains, avoid delay)
  if (myICM.status != ICM_20948_Stat_FIFOMoreDataAvail)
  {
    delay(10);
  }
}
