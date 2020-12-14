#include <iostream>
#include <cmath>

using namespace std;

float Kp1 = 200.0;
float Ki1 = 22.37;
float Kd1 = 23.79;

float Kp2 = 184.07;
float Ki2 = 50.0;
float Kd2 = 7.96;

float g = 9.81;

class two_link_planar_arm_model{
	float u1, u2;
	float c_x1, c_x2;
	float c_theta1, c_theta2;
	float c_v1, c_v2;
	float c_theta1_dot, c_theta2_dot;
	float c_theta1_dotdot, c_theta2_dotdot;
	float m1;
	float m2;
	float length1;
	float length2;
	two_link_planar_arm_model();
		
	public:
		two_link_planar_arm_model(float _m1, float _m2, float _length1, float _length2):m1(_m1),m2(_m2),length1(_length1),length2(_length2){
			cout << "Model intialized" << endl;
		}
		
		~two_link_planar_arm_model(){
			cout << "Model destroyed" << endl;
		}
		
		void run(float init_theta1, float init_theta2, float f_theta1, float f_theta2, float t, float step){
			// initializer
			c_x1 = c_x2 = 0;
			c_theta1 = init_theta1;
			c_theta2 = init_theta2;
			c_v1 = c_theta1_dot = 0;
			c_a1 = c_theta1_dotdot = 0;
			c_v2 = c_theta2_dot = 0;
			c_a2 = c_theta2_dotdot = 0;
			int count = 100; // print resolution

			float cumm_error1 = 0;
			float error1_dot = 0;
			float prev_error1 = 0;

			float cumm_error2 = 0;
			float error2_dot = 0;
			float prev_error2 = 0;
			
			for (float i = 1; i <= t; i += step){
				float error1 = (f_theta1 - c_theta1);
				float error2 = (f_theta2 - c_theta2);
				
				
				// B Matrix and its inverse
				float b11 = ((m1 + m2) * pow(length1, 2)) + (m2 * pow(length2, 2)) + (2 * m2 * length1 * length2 * cos(c_theta2)); 
				float b12 = (m2 * pow(length2, 2)) + (2 * m2 * length1 * length2 * cos(c_theta2));
				float b21 = (m2 * pow(length2, 2)) + (2 * m2 * length1 * length2 * cos(c_theta2));
				float b22 = m2 * pow(length2, 2);
				float b_det = (b11 * b22) - (b12 * b21);
				float inv_b11 = b22;
				float inv_b12 = -b12;
				float inv_b21 = -b21;
				float inv_b22 = b11;
				
				// C Matrix well its vector
				float c11 = - m2 * length1 * length2 * sin(c_theta2) * ((2 * c_theta1_dot * c_theta2_dot) + pow(c_theta2_dot,2)); 
				float c21 = - m2 * length1 * length2 * sin(c_theta2) * (c_theta1_dot * c_theta2_dot);				

				// G vector
				float g11 = -((m1 + m2) * g * length1 * sin(c_theta1)) - (m2 * g * length2 * sin(c_theta1 + c_theta2));
				float g21 = -(m2 * g * length2 * sin(c_theta1 + c_theta2));

				// PID controller where u = torque				
				u1 = Kp1 * error1 + Ki1 * cumm_error1 + Kd1 * error1_dot;
				u2 = Kp2 * error2 + Ki2 * cumm_error2 + Kd2 * error2_dot;
				
				// The dynamic model
				c_theta1_dotdot = (float((inv_b11 * (-c11-g11)) + (inv_b12 * (-c21-g21))) / float(b_det)) + u1;
				c_theta2_dotdot = (float((inv_b21 * (-c11-g11)) + (inv_b22 * (-c21-g21))) / float(b_det)) + u2;
				
				float p_theta1 = c_theta1;
				float p_theta2 = c_theta2;
				
				// Updating thetas for next iteration
				c_theta1 = c_theta1_dotdot  * pow(step , 2) + c_theta1;
				c_theta2 = c_theta2_dotdot  * pow(step , 2) + c_theta2;
				
				// Bounding angles usually [-pi, pi]
				if (c_theta1 >= M_PI/4){
					c_theta1 = M_PI/4;
				}
				if (c_theta1 <= 0){
					c_theta1 = 0;
				}

				if (c_theta2 >= M_PI/4){
					c_theta2 = M_PI/4;
				}
				if (c_theta2 <= 0){
					c_theta2 = 0;
				}
						
				// Updating theta_dots for next iteraion				
				c_theta1_dot = float(c_theta1 - p_theta1)/float(step);
				c_theta2_dot = float(c_theta2 - p_theta2)/float(step);

				// integratl and derivative of errors
				cumm_error1 += error1;
				error1_dot = float(error1 - prev_error1) / float(step);
				prev_error1 = error1;

				cumm_error2 += error2;
				error2_dot = float(error2 - prev_error2) / float(step);
				prev_error2 = error2;

				count--;
				if (count == 0){
					count = 100;
					// print theta in degrees (understandable)
					cout << "State X: " << c_theta1 * (float(180)/float(M_PI))  << ", " << c_theta2 * (float(180)/float(M_PI));
					cout << " and error: " << error1 << " , " << error2 << " at time t: " << i << endl;
				}
			}
		} 
};


int main(){
	two_link_planar_arm_model TLPAM(1.0, 1.0, 1.0, 1.0);
	TLPAM.run(float(0.1), float(0.1), float(0.5), float(0.5), float(100), float(0.01));

	return 0;
}
