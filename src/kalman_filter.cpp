#include "kalman_filter.h"
#include <iostream>
constexpr auto PI = 3.14159;

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
	MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
	x_ = x_in;
	P_ = P_in;
	F_ = F_in;
	H_ = H_in;
	R_ = R_in;
	Q_ = Q_in;
}

void KalmanFilter::Predict() {
	/**
	TODO:
	  * predict the state
	*/
	x_ = F_ * x_;
	MatrixXd Ft = F_.transpose();
	P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
	/**
	TODO:
	  * update the state by using Kalman Filter equations
	*/
	// standard KF equations
	VectorXd z_pred = H_ * x_;
	VectorXd y = z - z_pred;
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
	/**
	TODO:
	  * update the state by using Extended Kalman Filter equations
	*/
	// convert prediction state to polar using measurement function: h(x') = {h_rho; h_phi; h_rho_dot}
	float h_rho = sqrt(x_[0] * x_[0] + x_[1] * x_[1]);			// sqrt(px^2 + py^2)
	if (h_rho < 0.0001) { h_rho = 0.0001; }  // avoid division by zero
	float h_phi = atan2(x_[1], x_[0]);							// arctan(py/px)
	float h_rho_dot = (x_[0] * x_[2] + x_[1] * x_[3]) / h_rho;	// (px*vx +py*vy)/sqrt(px^2 + py^2)

	// create new z prediction with results from h(x')
	VectorXd z_pred = VectorXd(3);
	z_pred << h_rho, h_phi, h_rho_dot;	// z_pred <-- h(x')
	VectorXd y = z - z_pred;
	
	// normalize phi in 'y' between [-PI, PI]
	y[1] = y[1] - (2*PI * std::floor((y[1] + PI) / (2 * PI)));

	// standard KF equations
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;

}
