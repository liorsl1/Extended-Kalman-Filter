#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

// Constructor
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.P_ = MatrixXd(4, 4);
  // measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  // measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;
  // Projection matrix from 4d vector to 2d position vector
  H_laser_ << 1,0,0,0,
              0,1,0,0;
  // F state transition matrix
  ekf_.F_<< 1,0,1,0,
            0,1,0,1,
            0,0,1,0,
            0,0,0,1;
  // P state uncertainty matrix - for all variables in state (P(1,1) = uncertainty x,P(2,2) = uncertainty y,P(3,3) = uncertainty vx,P(4,4) = uncertainty vy)
  ekf_.P_ << 0.2, 0, 0, 0,
    0, 0.2, 0, 0,
    0, 0, 0.2, 0,
    0, 0, 0, 0.2;
  // set acceleration noise
  noise_ax = 9;
  noise_ay = 9;
}

// Destructor
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  
  if (!is_initialized_) {
    /**
     * TODO: Initialize the state ekf_.x_ with the first measurement.
     * TODO: Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // TODO: Convert radar from polar to cartesian coordinates 
      //         and initialize state.
      double m_Rho;
      double m_Theta;
      double m_RhoRate;
      double m_px;
      double m_py;
      double m_vx = 0;
      double m_vy = 0;

      m_Rho = measurement_pack.raw_measurements_[0];
      m_Theta = measurement_pack.raw_measurements_[1];
      m_RhoRate = measurement_pack.raw_measurements_[2];
      // Set theta to be in [-p,pi] range
      while(m_Theta < -M_PI){
        m_Theta = m_Theta + (2*M_PI);
      }
      while(m_Theta > M_PI){
          m_Theta = m_Theta - (2*M_PI);
      }

        m_px = m_Rho *sin(m_Theta);
        m_py = m_Rho *cos(m_Theta);
        m_vx = m_Rho * sin(m_Theta);
        m_vy = m_Rho * cos(m_Theta);
        ekf_.x_ << m_px,m_py,0,0;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }
    previous_timestamp_ = measurement_pack.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

   double d_t = (measurement_pack.timestamp_ - previous_timestamp_) / 10e5; //d_t - from miliseconds to seconds.
  // Update the state transition matrix F according to the new elapsed time.
  ekf_.F_(0,2) = d_t;
  ekf_.F_(1,3) = d_t;

  // Update the process noise covariance matrix.
  double dt_4 = pow(d_t,4)/4;
  double dt_3 = pow(d_t,3)/2;
  double dt_2 = pow(d_t,2);
  ekf_.Q_ << dt_4*noise_ax, 0, dt_3*noise_ax, 0,
                0, dt_4*noise_ay, 0, dt_3*noise_ay,
                dt_3*noise_ax, 0, dt_2*noise_ax, 0,
                0, dt_3*noise_ay, 0, dt_2*noise_ay;
  // ===== Prediction =====
  ekf_.Predict();
  // ===== Update =====
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates, update Hj, R_radar
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.R_ = R_radar_;
    // Update measurement for Radar (non-linear)
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    // Update measurement for Laser (linear)
    ekf_.Update(measurement_pack.raw_measurements_);
  }
  // Update previous to current time
    previous_timestamp_ = measurement_pack.timestamp_;
  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
