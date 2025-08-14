/* ============================================================================
   Teensy 4.1 — Thrust Vector Control + Reaction Wheel Flight Computer
   --------------------------------------------------------------------------
   Features:
     - State machine: GROUND, BOOST, COAST, APOGEE, DESCENT, LANDED
     - Two-axis TVC using servos (pitch & yaw)
     - Reaction wheel for roll stability
     - MPU9250 IMU with Kalman + complementary filtering
     - BMP388 barometer with oversampling & IIR filtering
     - PID loops for precise control
     - Recovery LED on apogee or failsafe tilt
     - Data logging to SD card
     - Noise minimization for both sensors
   --------------------------------------------------------------------------
   Hardware:
     Teensy 4.1
     MPU9250 IMU
     BMP388 barometer
     2x Servos for TVC
     DC motor + driver for reaction wheel
     LEDs for recovery indication & logging
   ============================================================================
*/

#include <Wire.h>
#include <SPI.h>
#include <Servo.h>
#include <SD.h>
#include <Adafruit_BMP3XX.h>
#include "MPU9250.h"

// --------------------- PIN CONFIG ---------------------
#define SERVO_PITCH_PIN 3
#define SERVO_YAW_PIN   4
#define RW_PWM_PIN      11
#define RW_DIR_PIN      9
#define SD_CS_PIN       10
#define RECOVERY_LED    12
#define LOGGING_LED     8
#define ARM_SWITCH_PIN  22
#define HEARTBEAT_LED   13

// --------------------- OBJECTS ------------------------
MPU9250 imu(Wire, 0x68);
Adafruit_BMP3XX bmp;
Servo servoPitch, servoYaw;
File logFile;

// --------------------- FLIGHT STATES ------------------
enum FlightState { GROUND, BOOST, COAST, APOGEE, DESCENT, LANDED };
FlightState state = GROUND;

// --------------------- FILTERS ------------------------
struct Kalman1D {
  float q, r, x, p;
  Kalman1D(): q(0.001f), r(0.03f), x(0), p(1) {}
  float update(float z, float dt) {
    p += q * dt;
    float k = p / (p + r);
    x += k * (z - x);
    p *= (1 - k);
    return x;
  }
};
Kalman1D kalPitch, kalRoll, kalAlt;

const float ALPHA = 0.98f; // complementary filter

// --------------------- PID CONTROL --------------------
struct PID {
  float kp, ki, kd;
  float iTerm, lastErr;
  float outMin, outMax;
  PID(): kp(0), ki(0), kd(0), iTerm(0), lastErr(0), outMin(-100), outMax(100) {}
  PID(float P, float I, float D, float mn, float mx) {
    kp=P; ki=I; kd=D; outMin=mn; outMax=mx;
    iTerm=0; lastErr=0;
  }
  float update(float setpoint, float measured, float dt) {
    float err = setpoint - measured;
    iTerm += err * dt;
    float dTerm = (err - lastErr) / dt;
    lastErr = err;
    float out = kp*err + ki*iTerm + kd*dTerm;
    if(out > outMax) out = outMax;
    if(out < outMin) out = outMin;
    return out;
  }
  void reset() { iTerm = 0; lastErr = 0; }
};
PID pidPitch(1.15f, 0.025f, 0.045f, -30, 30);
PID pidYaw(1.15f, 0.025f, 0.045f, -30, 30);
PID pidRoll(2.2f, 0.012f, 0.03f, -255, 255);

// --------------------- VARIABLES ----------------------
float groundPressure = 101325.0f;
float maxAltitude = 0.0f;
bool recoveryTriggered = false;
unsigned long boostStartTime = 0;

// --------------------- UTILS --------------------------
float pressureToAlt(float p, float p0) {
  return 44330.0f * (1.0f - pow(p/p0, 0.190294957f));
}
void logData(String line) {
  Serial.println(line);
  digitalWrite(LOGGING_LED, HIGH);
  logFile.println(line);
  logFile.flush();
  digitalWrite(LOGGING_LED, LOW);
}
void triggerRecovery(const char *reason) {
  if (!recoveryTriggered) {
    recoveryTriggered = true;
    digitalWrite(RECOVERY_LED, HIGH);
    logData(String("RECOVERY_TRIGGERED: ") + reason);
  }
}

// --------------------- SENSOR SETUP -------------------
void setupSensors() {
  if (imu.begin() != INV_SUCCESS) {
    Serial.println("MPU9250 init failed");
  } else {
    imu.setAccelRange(MPU9250::ACCEL_RANGE_8G);
    imu.setGyroRange(MPU9250::GYRO_RANGE_1000DPS);
    imu.setDlpfBandwidth(MPU9250::DLPF_BANDWIDTH_20HZ);
  }
  if (!bmp.begin_I2C()) {
    Serial.println("BMP388 init failed");
  } else {
    bmp.setTemperatureOversampling(BMP3_OVERSAMPLING_8X);
    bmp.setPressureOversampling(BMP3_OVERSAMPLING_8X);
    bmp.setIIRFilterCoeff(BMP3_IIR_FILTER_COEFF_15);
    bmp.setOutputDataRate(BMP3_ODR_50_HZ);
  }
}

// --------------------- SETUP --------------------------
void setup() {
  Serial.begin(115200);
  pinMode(HEARTBEAT_LED, OUTPUT);
  pinMode(RECOVERY_LED, OUTPUT);
  pinMode(LOGGING_LED, OUTPUT);
  pinMode(ARM_SWITCH_PIN, INPUT_PULLDOWN);
  pinMode(RW_PWM_PIN, OUTPUT);
  pinMode(RW_DIR_PIN, OUTPUT);

  servoPitch.attach(SERVO_PITCH_PIN);
  servoYaw.attach(SERVO_YAW_PIN);

  analogWriteFrequency(RW_PWM_PIN, 20000);

  if (SD.begin(SD_CS_PIN)) {
    logFile = SD.open("flightlog.csv", FILE_WRITE);
    logFile.println("time_ms,state,pitch,roll,yaw,ax,ay,az,gx,gy,gz,alt");
  } else {
    Serial.println("SD init failed");
  }

  setupSensors();

  if (bmp.performReading()) {
    groundPressure = bmp.pressure;
  }
}

// --------------------- STATE HANDLERS -----------------
void handleGround(unsigned long now) {
  static bool armed = false;

  // Step 1: Arm system when switch is flipped ON
  if (digitalRead(ARM_SWITCH_PIN) && !armed) {
    armed = true;
    logData("SYSTEM ARMED - Waiting for launch detection");
  }

  // Step 2: Detect launch only if armed
  if (armed) {
    // Read current acceleration magnitude
    float ax = imu.getAccelX_mss();
    float ay = imu.getAccelY_mss();
    float az = imu.getAccelZ_mss();
    float accelMag = sqrt(ax*ax + ay*ay + az*az);

    // If acceleration > 2.0g sustained for 100 ms, start BOOST phase
    static unsigned long accelHighSince = 0;
    if (accelMag > 19.62f) { // 2.0g in m/s²
      if (accelHighSince == 0) accelHighSince = now;
      if (now - accelHighSince > 100) {
        state = BOOST;
        boostStartTime = now;
        logData("STATE: BOOST (Launch detected)");
      }
    } else {
      accelHighSince = 0; // reset timer if acceleration drops
    }
  }
}
void handleBoost(unsigned long now, float accel, float alt) {
  if (now - boostStartTime > 1500 && accel < 6.0f) {
    state = COAST;
    logData("STATE: COAST");
  }
}
void handleCoast(unsigned long now, float alt) {
  if (alt > maxAltitude) maxAltitude = alt;
  else if (alt < maxAltitude - 0.5f) {
    state = APOGEE;
    triggerRecovery("Apogee detected");
    logData("STATE: APOGEE");
  }
}
void handleApogee() {
  state = DESCENT;
  logData("STATE: DESCENT");
}
void handleDescent(float alt) {
  if (alt < 2.0f) {
    state = LANDED;
    logData("STATE: LANDED");
  }
}
void handleLanded() {
  servoPitch.write(90);
  servoYaw.write(90);
  analogWrite(RW_PWM_PIN, 0);
}

// --------------------- CONTROL ------------------------
void controlTVC(float pitch, float roll, float yaw, float yawRate, float dt) {
  float pitchCmd = pidPitch.update(0.0f, pitch, dt);
  float yawCmd   = pidYaw.update(0.0f, roll, dt);
  servoPitch.write(constrain(90 + pitchCmd, 85, 95));
  servoYaw.write(constrain(90 + yawCmd, 85, 95));

  float rwOut = pidRoll.update(0.0f, yawRate, dt);
  int pwmVal = constrain(abs((int)rwOut), 0, 255);
  digitalWrite(RW_DIR_PIN, rwOut >= 0 ? HIGH : LOW);
  analogWrite(RW_PWM_PIN, pwmVal);
}

// --------------------- MAIN LOOP ----------------------
unsigned long lastTime = 0;
void loop() {
  unsigned long now = millis();
  float dt = (now - lastTime) / 1000.0f;
  lastTime = now;
  if (dt <= 0) dt = 0.001f;

  // Heartbeat
  digitalWrite(HEARTBEAT_LED, (now/500)%2);

  // IMU
  imu.readSensor();
  float ax = imu.getAccelX_mss(), ay = imu.getAccelY_mss(), az = imu.getAccelZ_mss();
  float gx = imu.getGyroX_rads() * 180/PI;
  float gy = imu.getGyroY_rads() * 180/PI;
  float gz = imu.getGyroZ_rads() * 180/PI;

  // Angles
  float pitchAcc = atan2(-ax, sqrt(ay*ay + az*az)) * 180/PI;
  float rollAcc  = atan2(ay, az) * 180/PI;
  static float pitchGyro=0, rollGyro=0, yawGyro=0;
  pitchGyro += gy * dt;
  rollGyro  += gx * dt;
  yawGyro   += gz * dt;
  float pitch = ALPHA*pitchGyro + (1-ALPHA)*pitchAcc;
  float roll  = ALPHA*rollGyro  + (1-ALPHA)*rollAcc;
  pitch = kalPitch.update(pitch, dt);
  roll  = kalRoll.update(roll, dt);

  // Barometer
  static float altitude = 0.0f;
  if (bmp.performReading()) {
    float altRaw = pressureToAlt(bmp.pressure, groundPressure);
    altitude = kalAlt.update(altRaw, dt);
  }

  float accelMag = sqrt(ax*ax + ay*ay + az*az);

  // State machine
  switch(state) {
    case GROUND:  handleGround(now); break;
    case BOOST:   handleBoost(now, accelMag, altitude); break;
    case COAST:   handleCoast(now, altitude); break;
    case APOGEE:  handleApogee(); break;
    case DESCENT: handleDescent(altitude); break;
    case LANDED:  handleLanded(); break;
  }

  // Control logic
  if (state == BOOST || state == COAST) controlTVC(pitch, roll, yawGyro, gz, dt);
  if (abs(pitch) > 30.0f) triggerRecovery("Tilt > 30 deg");

  // Log
  String logLine = String(now) + "," + String(state) + "," +
                   String(pitch,2) + "," + String(roll,2) + "," + String(yawGyro,2) + "," +
                   String(ax,2) + "," + String(ay,2) + "," + String(az,2) + "," +
                   String(gx,2) + "," + String(gy,2) + "," + String(gz,2) + "," +
                   String(altitude,2);
  logData(logLine);

  delay(2);
}
