#include <iostream>
#include <cmath>

using namespace std;

float Kp = 100.0;
float Ki = 2.0;
float Kd = 0.5;

class simple_torque_servo_motor_model{
	float u;
	float c_theta;
	float c_theta_dot;
	float c_theta_dotdot;
	float J;
	float B;
	simple_torque_servo_motor_model();
		
	public:
		simple_torque_servo_motor_model(float _J, float _B):J(_J), B(_B){
			cout << "Model intialized" << endl;
		}
		
		~simple_torque_servo_motor_model(){
			cout << "Model destroyed" << endl;
		}
		
		void run(float init_theta, float f_theta, float t, float step){
			c_theta = init_theta;
			c_theta_dot = 0;
			c_theta_dotdot = 0;
			int count = 100; // print resolution
			float cumm_error = 0;
			float error_dot = 0;
			float prev_error = 0;
			for (float i = 1; i <= t; i += step){
				// angle error in radians
				float error = (f_theta - c_theta);
				float deno = J;
				float a1 = -B; 
								
				// PID controller
				u = Kp * error + Ki * cumm_error + Kd * error_dot;

				// d^2(theta)/d^2(t) = -(B / J) * (d(theta) / d(t)) + ( u / J) where u = torque 
				c_theta_dotdot = (float(a1) / float(deno)) * c_theta_dot + (float(1) / float(deno)) * u;

				float prev_theta = c_theta;				
				float prev_theta_dot = c_theta_dot;
				
				// calculating updated theta
				c_theta = c_theta_dotdot  * pow(step , 2) + c_theta;
				
				// adjusting for radians
				while ((c_theta > (2 * M_PI)) || (c_theta < 0)){
					if (c_theta > (2 * M_PI)){
						c_theta = c_theta - (2 * M_PI);
					}
					if (c_theta < 0){
						c_theta = (2 * M_PI) + c_theta;
					}
				}
				
				// calculating theta_dot for next iteration
				c_theta_dot = float(c_theta - prev_theta) / float(step);
				
				// itegral error and derivative error
				cumm_error += error;
				error_dot = float(error - prev_error) / float(step);
				prev_error = error;
				count--;
				if (count == 0){
					count = 100;
					cout << "State X: " << c_theta << " , " << c_theta_dot << " and error: " << error << " at time t: " << i << endl;
				}
			}
		} 
};


int main(){
	simple_torque_servo_motor_model STSMM(0.02, 0.03);
	STSMM.run(float(0), float(1.0), float(10), float(0.001));

	return 0;
}
