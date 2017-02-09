#include <Servo.h>

const int MOTOR_PIN = 11;
const int SERVO_PIN= 12;

// Motor limits
// TODO are these the real limits?
//const int MOTOR_MAX = 120;
//const int MOTOR_MIN = 40;
const int MOTOR_NEUTRAL = 90;
// Optional: smaller values for testing safety
//const int MOTOR_MAX = 100; 
//const int MOTOR_MIN = 75;

// Steering limits
// TODO seems to me that the permissible steering range is really about [58,
// 120] judging from the sound of the servo pushing beyond a mechanical limit
// outside that range. The offset may be 2 or 3 deg and the d_theta_max is then
// ~31.
//const int D_THETA_MAX = 30;
//const int THETA_CENTER = 90;
//const int THETA_MAX = THETA_CENTER + D_THETA_MAX;
//const int THETA_MIN = THETA_CENTER - D_THETA_MAX;

// Interfaces to motor and steering actuators
Servo motor;
Servo steering;

// Initialize an instance of the Car class as car
//Car car;


void setup() {
  initActuators();
  armActuators();
}

void loop() {
  servoSweep();
}

/* Functions */
void initActuators() {
  motor.attach(MOTOR_PIN);
  steering.attach(SERVO_PIN);
}

void armActuators() {
  motor.write(MOTOR_NEUTRAL);
//  steering.write(THETA_CENTER);
  delay(1000);
}

void servoSweep(){
    for (pos = 0; pos <= 180; pos += 1) { // goes from 0 degrees to 180 degrees
    // in steps of 1 degree
    steering.write(pos);              // tell servo to go to position in variable 'pos'
    delay(15);                       // waits 15ms for the servo to reach the position
  }
  for (pos = 180; pos >= 0; pos -= 1) { // goes from 180 degrees to 0 degrees
    steering.write(pos);              // tell servo to go to position in variable 'pos'
    delay(15);                       // waits 15ms for the servo to reach the position
  }
}

