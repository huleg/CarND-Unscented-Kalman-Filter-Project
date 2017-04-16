#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
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
  // TODO: to fine-tune the practical value for a bicyle
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  // TODO: to fine-tune the practical value for a bicyle
  std_yawdd_ = 30;

  // Below deviation values should be left unchanged/untuned, as they are
  // provided by manufacturers who had done their calibration.
  //
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

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  // x_, time_us_ are left uninitialized in contructor.
  // They'll be initialized in first measurement arrives.
  is_initialized_ = false;

  n_x_ = 5; // px, py, velocity, psi, psi_dot
  n_aug_ = 7; // two extra augmented: nu_a(nu_acc_longitude), nu_yawdd(nu_acc_yaw)
  lambda_ = 3;  // Value from video. More on this, refer to
                // Slide 22 of http://ais.informatik.uni-freiburg.de/teaching/ws12/mapping/pdf/slam05-ukf.pdf

  // Co-variance matrix start with identity
  P_ << 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, 1, 0,
          0, 0, 0, 0, 1;

  // sigma points
  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);

  //set weights
  double sum_lambda_n_aug = lambda_+n_aug_;
  weights_ = VectorXd(2*n_aug_+1);
  weights_(0) = lambda_/sum_lambda_n_aug;
  //std::cout << "weights_: " << weights_(0);
  for (int i=1; i<2*n_aug_+1; i++){
    weights_(i) = 1/(2*sum_lambda_n_aug);
    //std::cout << " " << weights_(i);
  }
  //std::cout << "\n";

}

UKF::~UKF() {}

/**
 * ProcessFirstMeasurement
 * @param meas_package The 1st measurement data of either radar or laser
 * Return true if no error occurs, otherwise false.
 */
bool UKF::ProcessFirstMeasurement(MeasurementPackage meas_package) {
  // Initialize x_, Xsig_pred_, P_
  if(meas_package.sensor_type_ == MeasurementPackage::LASER) {
    x_(0) = meas_package.raw_measurements_(0); // px
    x_(1) = meas_package.raw_measurements_(1); // py
  } else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    x_(0) = meas_package.raw_measurements_(0)*cos(meas_package.raw_measurements_(1));
    x_(1) = meas_package.raw_measurements_(0)*sin(meas_package.raw_measurements_(1));
  } else {
    std::cerr << "1st measurement: Unrecognized type of measurement data!" << endl;
    return false;
  }

  x_(2) = 0; // longitudinal velocity
  x_(3) = 0; // psi. Since x_(2) = 0, so as rotational speed
  x_(4) = 0; // psi_dot.

  // Init timestamp
  time_us_ = meas_package.timestamp_;

  return true;
}

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

  // 1st measurement
  if (!is_initialized_) {
    is_initialized_ = UKF::ProcessFirstMeasurement(meas_package);
  } else {
    // Timestamp and delta T
    double delta_t = meas_package.timestamp_ - time_us_;
    time_us_ = meas_package.timestamp_;
    //std::cout << "delta T = " << delta_t << endl;

    // Prediction
    UKF::Prediction(delta_t);

    // Update
    if(meas_package.sensor_type_ == MeasurementPackage::LASER) {
      UKF::UpdateLidar(meas_package);
    } else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      UKF::UpdateRadar(meas_package);
    } else {
      std::cerr << "Unrecognized type of measurement data!" << endl;
    }
  }

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
  // Augmented sigma points generation
  // State prediction with sigma points
  // State prediction mean and co-variance
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
  // Measurements prediction
  // Update state
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
  // Measurements prediction
  // Update state
}
