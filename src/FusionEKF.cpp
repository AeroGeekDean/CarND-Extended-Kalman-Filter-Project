#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF. [done]
    * Set the process and measurement noises [done]
  */

  // KF measurement matrix (laser)
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  // init state covariance matrix
  ekf_.P_ = MatrixXd(4,4);
  ekf_.P_ << 1, 0, 0, 0,  // init guess. code will optimize
             0, 1, 0, 0,
             0, 0, 1000, 0,
             0, 0, 0, 1000;

  // init state transition matrix (will update w/ dt later)
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0,
             0, 1, 0, 1,
             0, 0, 1, 0,
             0, 0, 0, 1;

  // init process covariance matrix Q
  ekf_.Q_ = MatrixXd(4, 4);

  //set the acceleration noise components (for Q calc)
  noise_ax = 9;
  noise_ay = 9;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement. [done]
      * Create the covariance matrix. [done]
      * Remember: you'll need to convert radar from polar to cartesian coordinates. [done]
    */
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;
    ekf_.I_ = MatrixXd::Identity(4,4);

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to Cartesian coordinates and initialize state.
      */
      cout << "R:" << endl;
      float rho = measurement_pack.raw_measurements_[0];
      float theta = measurement_pack.raw_measurements_[1];
      float rhod = measurement_pack.raw_measurements_[2];

      ekf_.x_ << rho*cos(theta), // Px
                 rho*sin(theta), // Py
                 rhod*cos(theta), // Vx
                 rhod*sin(theta); // Vy
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      cout << "L:" << endl;
      ekf_.x_ << measurement_pack.raw_measurements_[0], // Px
                 measurement_pack.raw_measurements_[1], // Py
                 0, 0; // Vx, Vy
    }

    cout << "x_(init) = " << ekf_.x_ << endl;

    // save timestamp for next frame
    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time. [done]
      - Time is measured in seconds. [done]
     * Update the process noise covariance matrix. [done]
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix. [done]
   */
  float dt =  (measurement_pack.timestamp_-previous_timestamp_)/1000000.0; // [sec]

  // save timestamp for next frame
  previous_timestamp_ = measurement_pack.timestamp_;

  // update state transition matrix F
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  // update process covariance matrix Q
  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;
  ekf_.Q_ << dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
            0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
            dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
            0, dt_3/2*noise_ay, 0, dt_2*noise_ay;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step. [done]
     * Update the state and covariance matrices. [done]
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    cout << "R:" << endl;
    ekf_.R_ = R_radar_;
    Hj_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.H_ = Hj_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    cout << "L:" << endl;
    ekf_.R_ = R_laser_;
    ekf_.H_ = H_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }
 
  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
