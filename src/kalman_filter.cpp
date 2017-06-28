#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

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
  long x_size = x_.size();
  I_ = MatrixXd::Identity(x_size,x_size);
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state [done]
  */
  x_ = F_ * x_;
  P_ = F_*P_*F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations [done]
  */
  VectorXd y = z - H_*x_;

  // call the common part of update
  UpdateCommon(y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  
  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);

  // save off temp vars
  float c1 = px*px+py*py;
  float c2 = sqrt(c1);
  
  // Apply non-linear h(x) function
  VectorXd hx = VectorXd(3,1);
  hx << c2,
        atan2(py, px),
        (px*vx+py*vy)/c2;

  VectorXd y = z - hx;

  // keep phi angle y(1) between +/- pi
  while (y(1) < -M_PI) y(1) += M_PI;
  while (y(1) >  M_PI) y(1) -= M_PI;

  // call the common part of update
  UpdateCommon(y);
}

void KalmanFilter::UpdateCommon(const Eigen::VectorXd &y) {

  MatrixXd S = H_*P_*H_.transpose() + R_;
  MatrixXd K = P_*H_.transpose()*S.inverse();
  
  // new state
  x_ = x_ + K*y;
  P_ = (I_-K*H_)*P_;

}

