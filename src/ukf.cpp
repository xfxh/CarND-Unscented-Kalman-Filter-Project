#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 2;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  is_initialized_ = false;
  time_us_ = 0;
  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3-n_aug_;
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  weights_ = VectorXd(2*n_aug_+1);
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  if(!is_initialized_)
  {
      P_<<1,0,0,0,0,
          0,1,0,0,0,
          0,0,1,0,0,
          0,0,0,1,0,
          0,0,0,0,1;
      if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
      {
          x_<<meas_package.raw_measurements_[0]*cos(meas_package.raw_measurements_[1]),
              meas_package.raw_measurements_[0]*sin(meas_package.raw_measurements_[1]),
              0,0,0;  
      }
      else if (meas_package.sensor_type_ == MeasurementPackage::LASER)
      {
          x_<<meas_package.raw_measurements_[0],meas_package.raw_measurements_[1],0,0,0;
      }
      else
      {
          //do nothing
      }
      time_us_ = meas_package.timestamp_;
      is_initialized_ = true;
      return;
  }
  
  float dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;
  Prediction(dt);
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
  {
     UpdateRadar(meas_package);
  }
  else
  {
     UpdateLidar(meas_package);
  }

  // print the output
  cout << "x_ = " << x_ << endl;
  cout << "P_ = " << P_ << endl;
  
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  //Generating Sigma Points - Augmentation
  //create augmented mean vector
  VectorXd x_aug = VectorXd(7);
  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);
  //create example sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  //create augmented mean state
  x_aug.head(n_x_) = x_;
  x_aug(n_x_) = 0;
  x_aug(n_x_+1) = 0;
  //create augmented covariance matrix
  P_aug = MatrixXd::Zero(7,7);
  P_aug.topLeftCorner(n_x_,n_x_) = P_;
  MatrixXd Q = MatrixXd(2, 2);
  Q<< std_a_*std_a_,0,0,std_yawdd_*std_yawdd_;
  P_aug.bottomRightCorner(2,2) = Q;
  //create square root matrix
  MatrixXd A = P_aug.llt().matrixL();
  //create augmented sigma points
   Xsig_aug.col(0) = x_aug;
  for(int i=0; i<n_aug_; i++)
  {
      Xsig_aug.col(i+1) = x_aug + sqrt(lambda_+n_aug_)*A.col(i);
      Xsig_aug.col(i+n_aug_+1) = x_aug - sqrt(lambda_+n_aug_)*A.col(i);
  }

  //Sigma Points Prediction
  for(int i=0;i<2 * n_aug_ + 1;i++)
  {

      VectorXd a1=VectorXd(5);
      VectorXd b1=VectorXd(5);
      if(Xsig_aug.col(i)(n_x_-1) != 0)
      {
         a1<< Xsig_aug.col(i)(2)*(sin(Xsig_aug.col(i)(3)+Xsig_aug.col(i)(4)*delta_t)-sin(Xsig_aug.col(i)(3)))/Xsig_aug.col(i)(4),
              Xsig_aug.col(i)(2)*(cos(Xsig_aug.col(i)(3))-cos(Xsig_aug.col(i)(3)+Xsig_aug.col(i)(4)*delta_t))/Xsig_aug.col(i)(4),
              0,Xsig_aug.col(i)(4)*delta_t,0;
      }
      else
      {
          a1<< Xsig_aug.col(i)(2)*cos(Xsig_aug.col(i)(3))*delta_t,
               Xsig_aug.col(i)(2)*sin(Xsig_aug.col(i)(3))*delta_t,
               0,Xsig_aug.col(i)(4)*delta_t,0;
          
      }
      b1<< 0.5*delta_t*delta_t*cos(Xsig_aug.col(i)(3))*Xsig_aug.col(i)(5),
              0.5*delta_t*delta_t*sin(Xsig_aug.col(i)(3))*Xsig_aug.col(i)(5),
              Xsig_aug.col(i)(5)*delta_t,
              0.5*delta_t*delta_t*Xsig_aug.col(i)(6),
              delta_t*Xsig_aug.col(i)(6);
      Xsig_pred_.col(i) = Xsig_aug.col(i).head(n_x_)+a1+b1;
  }

  //Predicted Mean and Covariance
  //create vector for predicted state
  VectorXd x = VectorXd(n_x_);
  //create covariance matrix for prediction
  MatrixXd P = MatrixXd(n_x_, n_x_);
  //set weights
  weights_(0)=lambda_/(lambda_+n_aug_);
  for(int i=1;i<2*n_aug_+1;i++)
  {
      weights_(i)=1/(2*(lambda_+n_aug_));
  }
  //predict state mean
  x<<0,0,0,0,0;
  for(int i=0;i<2*n_aug_+1;i++)
  {
      x=x+weights_(i)*Xsig_pred_.col(i);
  }
  //predict state covariance matrix
  P<<0,0,0,0,0,
     0,0,0,0,0,
     0,0,0,0,0,
     0,0,0,0,0,
     0,0,0,0,0;
  for(int i=0;i<2*n_aug_+1;i++)
  {
      VectorXd x_diff = Xsig_pred_.col(i) - x;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

      P=P+weights_(i)*x_diff*x_diff.transpose();
  }

  x_ = x;
  P_ = P;

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 2;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z); 
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  
  //transform sigma points into measurement space
  for(int i=0;i<2 * n_aug_ + 1;i++)
  {
      VectorXd a=VectorXd(n_z);
      a<< Xsig_pred_.col(i)(0),
          Xsig_pred_.col(i)(1);
      Zsig.col(i)=a;
  }
  //calculate mean predicted measurement
  z_pred<<0,0;
  for(int i=0;i<2 * n_aug_ + 1;i++)
  {
      z_pred = z_pred+weights_(i)*Zsig.col(i);
  }
  //calculate innovation covariance matrix S
  S<<0,0,0,0;
  for(int i=0;i<2 * n_aug_ + 1;i++)
  {
       //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }
  MatrixXd R = MatrixXd(n_z,n_z);
  R<<std_laspx_*std_laspx_,0,
     0,std_laspy_*std_laspy_;
  S= S+R;
  
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  //calculate cross correlation matrix
  Tc.fill(0);
  for(int i=0;i<2 * n_aug_ + 1;i++)
  {
      //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    Tc=Tc+weights_(i)*x_diff*z_diff.transpose();
  }
  //calculate Kalman gain K;
  MatrixXd K = MatrixXd(n_x_, n_z);
  K=Tc*S.inverse();
  //update state mean and covariance matrix
  //residual
  VectorXd z_diff = meas_package.raw_measurements_ - z_pred;

  x_=x_+K*z_diff;
  P_=P_-K*S*K.transpose();

  //NIS
  e_l=z_diff.transpose()*S.inverse()*z_diff;
   // print the NIS
  cout << "Lidar NIS e = " << e_l << endl;

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z); 
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);


  //transform sigma points into measurement space
  for(int i=0;i<2 * n_aug_ + 1;i++)
  {
      VectorXd a=VectorXd(n_z);
      a<< sqrt(Xsig_pred_.col(i)(0)*Xsig_pred_.col(i)(0)+Xsig_pred_.col(i)(1)*Xsig_pred_.col(i)(1)),
          atan2(Xsig_pred_.col(i)(1),Xsig_pred_.col(i)(0)),
          (Xsig_pred_.col(i)(0)*cos(Xsig_pred_.col(i)(3))*Xsig_pred_.col(i)(2)+Xsig_pred_.col(i)(1)*sin(Xsig_pred_.col(i)(3))*Xsig_pred_.col(i)(2))/sqrt(Xsig_pred_.col(i)(0)*Xsig_pred_.col(i)(0)+Xsig_pred_.col(i)(1)*Xsig_pred_.col(i)(1));
      Zsig.col(i)=a;
  }
  //calculate mean predicted measurement
  z_pred<<0,0,0;
  for(int i=0;i<2 * n_aug_ + 1;i++)
  {
      z_pred = z_pred+weights_(i)*Zsig.col(i);
  }
  //calculate innovation covariance matrix S
  S<<0,0,0,
     0,0,0,
     0,0,0;
  for(int i=0;i<2 * n_aug_ + 1;i++)
  {
       //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }
  MatrixXd R = MatrixXd(n_z,n_z);
  R<<std_radr_*std_radr_,0,0,
     0,std_radphi_*std_radphi_,0,
     0,0,std_radrd_*std_radrd_;
  S= S+R;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  //calculate cross correlation matrix
  Tc.fill(0);
  for(int i=0;i<2 * n_aug_ + 1;i++)
  {
      //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
      Tc=Tc+weights_(i)*x_diff*z_diff.transpose();
  }
  //calculate Kalman gain K;
  MatrixXd K = MatrixXd(n_x_, n_z);
  K=Tc*S.inverse();
  //update state mean and covariance matrix
  //residual
  VectorXd z_diff = meas_package.raw_measurements_ - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
  x_=x_+K*z_diff;
  P_=P_-K*S*K.transpose();

  //NIS
  e_r=z_diff.transpose()*S.inverse()*z_diff;
  // print the NIS
  cout << "Radar NIS e = " << e_r << endl;

}
