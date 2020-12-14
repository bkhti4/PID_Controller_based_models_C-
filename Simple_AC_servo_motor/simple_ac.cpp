#include <iostream>
#include <cmath>

using namespace std;

float Kp = 100;
float Ki = 10;
float Kd = 40;

class simple_ac_servo_motor_model{
	float u;
	float c_x;
	float c_theta;
	float c_theta_dot;
	float c_theta_dotdot;
	float N;
	float J_m;
	float J_g;
	float J_l;
	float K1;
	float K2;
	float c;
	simple_ac_servo_motor_model();
		
	public:
		simple_ac_servo_motor_model(float _N, float _J_m, float _J_g, float _J_l, float _K1, float _K2, float _c):N(_N)
			,J_m(_J_m),J_g(_J_g),J_l(_J_l),K1(_K1),K2(_K2),c(_c){
			cout << "Model intialized" << endl;
		}
		
		~simple_ac_servo_motor_model(){
			cout << "Model destroyed" << endl;
		}
		
		void run(float init_theta, float f_theta, float t, float step){
			c_theta = init_theta;
			c_theta_dot = 0;
			c_theta_dotdot = 0;
			int count = 1000; // print resolution
			float cumm_error = 0;
			float error_dot = 0;
			float prev_error = 0;
			for (float i = 1; i <= t; i += step){
				// angle error in radians
				float error = (f_theta - c_theta);
				
				// d^2(theta)/d^2(t) = (1 / (J_m + J_g + (J_l/N^2)) ) [ -K2 * d(theta)/d(t) + K1 * u] where u = Stator Voltage
				float deno = (J_m + J_g + (float(J_l) / float(N * N)));
				float a1 = -K2; 
				float a2 = K1;
				// PID controller				
				u = Kp * error + Ki * cumm_error + Kd * error_dot;
				c_theta_dotdot = (float(a1) / float(deno)) * c_theta_dot + (float(a2) / float(deno)) * u;

				float prev_theta = c_theta;				
				float prev_theta_dot = c_theta_dot;
				
				// updating theta
				c_theta = c_theta_dotdot  * pow(step , 2) + c_theta;
				
				// to bound rotation but f_theta should be adjusted
				//while ((c_theta > (2 * M_PI)) || (c_theta < 0)){
				//	if (c_theta > (2 * M_PI)){
				//		c_theta = c_theta - (2 * M_PI);
				//	}
				//	if (c_theta < 0){
				//		c_theta = (2 * M_PI) + c_theta;
				//	}
				//}
				
				// updating theta_dot for next iteration
				c_theta_dot = float(c_theta - prev_theta) / float(step);
				
				// integral and derivative of error
				cumm_error += error;
				error_dot = float(error - prev_error) / float(step);
				prev_error = error;
				count--;
				if (count == 0){
					count = 1000;
					cout << "State X: " << c_theta << " , " << c_theta_dot << " and error: " << error << " at time t: " << i << endl;
				}
			}
		} 
};


int main(){
	simple_ac_servo_motor_model SASMM(5.6666, 0.00012, 0.01, 0.00206, 0.00843, 0.00703, 0.00); 
	SASMM.run(float(0), float(6.0), float(1000), float(0.001));

	return 0;
}
