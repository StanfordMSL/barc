// basic hall effect sensor test 


volatile byte half_revolutions;
unsigned int rpm;
unsigned long timeold;

const int HALL_SENSOR_PIN = 2;

void setup()
{
 Serial.begin(9600);
 attachInterrupt(digitalPinToInterrupt(HALL_SENSOR_PIN), magnet_detect, RISING);//Initialize the intterrupt pin (Arduino digital pin 2)
 half_revolutions = 0;
 rpm = 0;
 timeold = 0;
 Serial.println("hall effect sensor test");
}
void loop()//Measure RPM
{
// if (half_revolutions >= 20) { 
//   rpm = 30*1000/(millis() - timeold)*half_revolutions;
//   timeold = millis();
//   half_revolutions = 0;
   //Serial.println(rpm,DEC);
// }
}
void magnet_detect()//ISR
{
 half_revolutions++;
 Serial.println((int) half_revolutions, DEC);
}
