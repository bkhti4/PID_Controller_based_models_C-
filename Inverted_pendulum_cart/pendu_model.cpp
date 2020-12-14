#include <iostream>
#include <cmath>

using namespace std;

float Kp = 100.0;
float Ki = 1.0;
float Kd = 20.0;
float g = 9.81;

class inverted_pendulum_cart_model{
	float u;
	float c_x;
	float c_theta;
	float c_v;
	float c_theta_dot;
	float c_a;
	float c_theta_dotdot;
	float m_pend;
	float M_cart;
	float length;
	float I;
	float f_b;
	inverted_pendulum_cart_model();
		
	public:
		inverted_pendulum_cart_model(float _m_pend, float _M_cart, float _length, float _I, float _f_b):m_pend(_m_pend),M_cart(_M_cart),length(_length),I(_I)
		,f_b(_f_b){
			cout << "Model intialized" << endl;
		}
		
		~inverted_pendulum_cart_model(){
			cout << "Model destroyed" << endl;
		}
		
		void run(float init_x, float init_theta, float f_theta, float t, float step){
			c_x= init_x;
			c_theta = init_theta;
			c_v = c_theta_dot = 0;
			c_a = c_theta_dotdot = 0;
			int count = 100;
			float cumm_error = 0;
			float error_dot = 0;
			float prev_error = 0;
			for (float i = 1; i <= t; i += step){
				float error = (f_theta - c_theta);
				float deno = I * (M_cart + m_pend) + (M_cart * m_pend * length * length);
				float a1 = -(I + (m_pend *length * length)) * f_b; 
				float a2 = (m_pend * m_pend * length * length * g);				
				float a3 = (I + (m_pend *length * length));
				float a4 = -(m_pend *length * f_b);
				float a5 = (M_cart + m_pend) * (m_pend * g * length);
				float a6 = (m_pend * length);
								
				u = Kp * error + Ki * cumm_error + Kd * error_dot;
				c_a =  (float(a1) / float(deno)) * c_v + (float(a2) / float(deno)) * c_theta + (float(a3) / float(deno)) * u;
				c_theta_dotdot = (float(a4) / float(deno)) * c_v + (float(a5) / float(deno)) * c_theta + (float(a6) / float(deno)) * u;
				
				float p_x = c_x;
				float p_theta = c_theta;
				c_x = c_a  * pow(step , 2) + c_x;
				c_theta = c_theta_dotdot  * pow(step , 2) + c_theta;
				if (c_theta >= M_PI/4){
					c_theta = M_PI/4;
				}
				if (c_theta <= -M_PI/4){
					c_theta = -M_PI/4;
				}				
				c_v = float(c_x - p_x)/float(step);
				c_theta_dot = float(c_theta - p_theta)/float(step);
				cumm_error += error;
				error_dot = float(error - prev_error) / float(step);
				prev_error = error;
				count--;
				if (count == 0){
					count = 100;
					cout << "State X: " << c_x  << ", " << c_theta * (float(180)/float(M_PI)) << " and error: " << error << " at time t: " << i << endl;
				}
			}
		} 
};


int main(){
	inverted_pendulum_cart_model P1(0.2, 0.5, 0.3, 0.006, 0.1);
	P1.run(float(10), float(0.1), float(0), float(10), float(0.01));

	return 0;
}
