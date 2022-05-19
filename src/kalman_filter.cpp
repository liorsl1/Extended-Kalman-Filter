#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

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
   * TODO: predict the state - Prediction step is the final step,
   * Where we output new x, and the final P uncertainty matrix.
   */
    x_ = F_ * x_; // noise = 0
    P_ = F_* P_* F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * TODO: update the state by using Kalman Filter equations
   */
    
    VectorXd z_pred = H_*x_;
    VectorXd y = z-z_pred;
    MeasureUpdate(y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */
  Eigen::VectorXd z_pred;
  z_pred = VectorXd(3);
  /*
    Creating Hj Matrix according to h(x), with polar coordinates of x.
    So, h(x) outputs the predicted z state in polar form.
  */
  z_pred(0) = (sqrt(x_(0)*x_(0)+x_(1)*x_(1)));
  z_pred(1) = (atan2(x_(1),x_(0)));
  
  if(abs(z(1) - z_pred(1)) > M_PI) {
      if(z_pred(1) > 0){
        z_pred(1) = z_pred(1) - 2*M_PI;
      }
      else {
        z_pred(1) = z_pred(1) + 2*M_PI;
      }
  }
	z_pred(2) = ((x_(0)*x_(2)+x_(1)*x_(3))/z_pred(0));

	VectorXd y = z - z_pred;
    MeasureUpdate(y);
}

void KalmanFilter::MeasureUpdate(const VectorXd &y) {
  /**
   * TODO: update the measurement according to the correct sensor's y measure.
   */
    MatrixXd K,S;
    S = H_ * P_ * H_.transpose() + R_;
    K = P_ * H_.transpose() * S.inverse();
    // New updated State of x (estimated velocity and position) and P (uncertainty)
    x_ = x_ + (K*y);
    int x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K*H_)*P_;
}


