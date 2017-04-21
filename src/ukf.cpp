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
  std_a_ = 0.8; // roughly 0.4m/s of speed of error. In sample-1.txt, speed is roughly 1m/s

  // Process noise standard deviation yaw acceleration in rad/s^2
  // TODO: to fine-tune the practical value for a bicyle
  std_yawdd_ = 0.5; // roughly 25 degrees per second change

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
  lambda_ = 3-n_aug_;  // Value from video. More on this, refer to
                // Slide 22 of http://ais.informatik.uni-freiburg.de/teaching/ws12/mapping/pdf/slam05-ukf.pdf

  // Co-variance matrix start with identity
  P_ << 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, 1, 0,
          0, 0, 0, 0, 1;

  // TODO: sigma points initialized in Prediction step
  // Xsig_pred_
  // Sigma points
  Xsig_pred_ = MatrixXd(n_aug_, 2*n_aug_+1);
  Xsig_pred_.fill(0);
  std::cout << "Xsig_pred_ = " << Xsig_pred_ << endl;

  //set weights
  sum_lambda_n_aug_ = lambda_+n_aug_;
  weights_ = VectorXd(2*n_aug_+1);
  weights_(0) = lambda_/sum_lambda_n_aug_;
  //std::cout << "weights_: " << weights_(0);
  for (int i=1; i<2*n_aug_+1; i++){
    weights_(i) = 0.5/sum_lambda_n_aug_;
    //std::cout << " " << weights_(i);
  }
  //std::cout << "\n";

  NIS_laser_ = 0.0;
  NIS_radar_ = 0.0;

}

UKF::~UKF() {}

/**
 * Prediction step 1
 */
void UKF::GenerateSigmaPoints() {
  // Augmented state vector
  VectorXd *x_aug = new VectorXd(n_aug_);
  // Augmented covariance matrix
  MatrixXd *P_aug = new MatrixXd(n_aug_, n_aug_);

  //create augmented mean state
  x_aug->head(5) = x_;
  x_aug->coeffRef(5) = 0; // nu_a
  x_aug->coeffRef(6) = 0; // nu_yawdd

  //create augmented covariance matrix
  P_aug->fill(0.0);
  P_aug->topLeftCorner(5,5) = P_;
  P_aug->coeffRef(5,5) = std_a_*std_a_;
  P_aug->coeffRef(6,6) = std_yawdd_*std_yawdd_;

  //calculate square root of P - Chelosky decomposition
  MatrixXd A = P_aug->llt().matrixL();
  std::cout << "GenSigmaPts: A = " << A << std::endl;

  //set first column of sigma point matrix
  Xsig_pred_.fill(0.0); // Clear data from previous round, esp. reusing it for predicted X.
  Xsig_pred_.col(0)  = *x_aug;
  //set remaining sigma points - reusing X state prediction matrix
  double sqrt_sum_lambda_n_aug = sqrt(sum_lambda_n_aug_);
  for (int i = 0; i < n_aug_; i++)
  {
    Xsig_pred_.col(i+1)     = *x_aug + sqrt_sum_lambda_n_aug * A.col(i);
    Xsig_pred_.col(i+1+n_aug_) = *x_aug - sqrt_sum_lambda_n_aug * A.col(i);
  }
  // Debug
  //std::cout << "GenSigmaPts: Xsig_pred_ = " << std::endl << Xsig_pred_ << std::endl;
}

/**
 * Prediction step 2 - State Prediction
 *   Re-using Xsig_pred_ to save predicted states. It was saving sigma points.
 */
void UKF::StatePrediction(double delta_t) {
  //State prediction with sigma points
  for (int i = 0; i< 2*n_aug_+1; i++)
  {
    //extract values for better readability
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);
    double yawd = Xsig_pred_(4,i);
    double nu_a = Xsig_pred_(5,i);
    double nu_yawdd = Xsig_pred_(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
    // TODO: nu_a, nu_yawdd seems no need to write back to Xsig_pred_
  }
  // Debug
  //std::cout << "StatePred: Xsig_pred_ = " << Xsig_pred_ << std::endl;
}

/**
 * Prediction step 3 - Predicated states mean and covariance
 *   Directly update x_ and P_, both have only n_x_ elements in each col.
 *   However, Xsig_pred_ has n_aug_ elements each col. Pay attention to
 *   dimention change.
 */
void UKF::StatePredictionMeanCovariance() {
  //predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x_ += weights_(i) * Xsig_pred_.col(i).segment(0,5);
  }

  //predicted state covariance matrix
  P_.fill(0.0);
  std::cout << "P_ = " << P_ << std::endl;
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i).segment(0,5) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P_ += weights_(i) * x_diff * x_diff.transpose() ;
    //std::cout << "P_ = " << P_ << std::endl;
  }
  std::cout << "P_ = " << P_ << std::endl;
}

/**
 * ProcessFirstMeasurement
 * @param meas_package The 1st measurement data of either radar or laser
 * Return true if no error occurs, otherwise false.
 */
bool UKF::ProcessFirstMeasurement(MeasurementPackage meas_package) {
  // Initialize x_
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
  UKF::GenerateSigmaPoints();
  // State prediction with sigma points
  UKF::StatePrediction(delta_t);
  // State prediction mean and co-variance
  UKF::StatePredictionMeanCovariance();
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

  // Lidar measurements are from linear model, thus simpler measurements prediction.
  double z_pred_px = x_(0);
  double z_pred_py = x_(1);
  // Z pred covariance
  MatrixXd z_pred_covar = MatrixXd(2, 2);
  z_pred_covar = P_.topLeftCorner(2, 2);
  // Add measurement error covariances
  z_pred_covar(0,0) += std_laspx_ * std_laspx_;
  z_pred_covar(1,1) += std_laspy_ * std_laspy_;

  // Update state

  // Kalman gain
  MatrixXd T = MatrixXd(5,2);
  T.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i).segment(0,5) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    T += weights_(i) * x_diff * x_diff.segment(0,2).transpose() ;
  }
  //std::cout << "T = " << T << std::endl;

  MatrixXd K = MatrixXd(5,2);
  MatrixXd z_pred_covar_inv = z_pred_covar.inverse();
  std::cout << "z_pred_covar = " << z_pred_covar << std::endl;
  K = T * z_pred_covar_inv;
  //std::cout << "K = " << K << std::endl;

  // State update with measured data
  VectorXd z_err = VectorXd(2);
  z_err(0) = meas_package.raw_measurements_(0) - z_pred_px;
  z_err(1) = meas_package.raw_measurements_(1) - z_pred_py;
  x_ += K*z_err;
  std::cout << "x_ = " << x_ << std::endl;

  // State covariance matrix update
  P_ -= K*z_pred_covar*K.transpose();
  std::cout << "P_ = " << P_ << std::endl;

  // NIS
  NIS_laser_ = z_err.transpose()*z_pred_covar_inv*z_err;
  std::cout << "Radar NIS: NIS_laser_ = " << NIS_laser_ << std::endl;
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
  VectorXd z_pred_rho = VectorXd(15);
  VectorXd z_pred_phi = VectorXd(15);
  VectorXd z_pred_rho_dot = VectorXd(15);

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    double sqrt_px2_n_py2 = sqrt(Xsig_pred_.col(i)[0]*Xsig_pred_.col(i)[0] + Xsig_pred_.col(i)[1]*Xsig_pred_.col(i)[1]);
    if (sqrt_px2_n_py2 < 0.001) // If too small, avoid div-by-0.
        sqrt_px2_n_py2 = 0.001;
    z_pred_rho[i] = sqrt_px2_n_py2;
    z_pred_phi[i] = atan2(Xsig_pred_.col(i)[1], Xsig_pred_.col(i)[0]);
    double v1 = cos(Xsig_pred_.col(i)[3])*Xsig_pred_.col(i)[2];
    double v2 = sin(Xsig_pred_.col(i)[3])*Xsig_pred_.col(i)[2];
    z_pred_rho_dot[i] = (Xsig_pred_.col(i)[0]*v1 +
            Xsig_pred_.col(i)[1]*v2)/sqrt_px2_n_py2;
  }
  std::cout << "Radar update: z_pred_rho = " << z_pred_rho << std::endl;
  std::cout << "Radar update: z_pred_phi = " << z_pred_phi << std::endl;
  std::cout << "Radar update: z_pred_rho_dot = " << z_pred_rho_dot << std::endl;

  // z_pred mean calculation.
  double z_pred_rho_mean = 0.0;
  double z_pred_phi_mean = 0.0;
  double z_pred_rho_dot_mean = 0.0;

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    z_pred_rho_mean += weights_(i) * z_pred_rho[i];
    z_pred_phi_mean += weights_(i) * z_pred_phi[i];
    z_pred_rho_dot_mean += weights_(i) * z_pred_rho_dot[i];
  }
  std::cout << "Radar update: z_pred_rho_mean = " << z_pred_rho_mean << std::endl;
  std::cout << "Radar update: z_pred_phi_mean = " << z_pred_phi_mean << std::endl;
  std::cout << "Radar update: z_pred_rho_dot_mean = " << z_pred_rho_dot_mean << std::endl;

  // z_pred_covar calculation.
  MatrixXd z_pred_covar = MatrixXd(3, 3);
  z_pred_covar.fill(0.0);
  std::cout << "Radar update: z_pred_covar(init) = " << z_pred_covar << std::endl;

  for (int i = 1; i < 2 * n_aug_ + 1; i++) {
    // state difference
    VectorXd z_diff = VectorXd(3);
    z_diff(0) = z_pred_rho[i]-z_pred_rho_mean;
    z_diff(1) = z_pred_phi[i]-z_pred_phi_mean;
    z_diff(2) = z_pred_rho_dot[i]-z_pred_rho_dot_mean;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    z_pred_covar += weights_(i) * z_diff * z_diff.transpose() ;
    //std::cout << "Radar update: z_pred_covar = " << z_pred_covar << std::endl;
  }

  // Add measurement error covariances
  z_pred_covar(0,0) += std_radr_ * std_radr_;
  z_pred_covar(1,1) += std_radphi_ * std_radphi_;
  z_pred_covar(2,2) += std_radrd_* std_radrd_;
  std::cout << "Radar update: z_pred_covar (final) = " << z_pred_covar << std::endl;


  // Update state
  MatrixXd z_pred_covar_inv = z_pred_covar.inverse();
  // Kalman gain
  MatrixXd T = MatrixXd(5,3);
  T.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    // measurement state difference
    VectorXd z_diff = VectorXd(3);
    z_diff(0) = z_pred_rho[i]-z_pred_rho_mean;
    z_diff(1) = z_pred_phi[i]-z_pred_phi_mean;
    z_diff(2) = z_pred_rho_dot[i]-z_pred_rho_dot_mean;
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i).segment(0,5) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    T += weights_(i) * x_diff * z_diff.transpose() ;
  }
  std::cout << "Radar update: T = " << T << std::endl;

  MatrixXd K = MatrixXd(5,3);


  K = T * z_pred_covar_inv;
  std::cout << "Radar update: K = " << K << std::endl;

  // State update with measured data
  VectorXd z_err = VectorXd(3);
  z_err(0) = meas_package.raw_measurements_(0) - z_pred_rho_mean;
  z_err(1) = meas_package.raw_measurements_(1) - z_pred_phi_mean;
  z_err(2) = meas_package.raw_measurements_(2) - z_pred_rho_dot_mean;
  x_ += K*z_err;
  std::cout << "Radar update: x_ = " << x_ << std::endl;

  // State covariance matrix update
  P_ -= K*z_pred_covar*K.transpose();
  std::cout << "Radar update: P_ = " << P_ << std::endl;

  // NIS
  NIS_radar_ = z_err.transpose()*z_pred_covar_inv*z_err;
  std::cout << "Radar NIS: NIS_radar_ = " << NIS_radar_ << std::endl;

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
    // Update
    if(use_laser_ && (meas_package.sensor_type_ == MeasurementPackage::LASER)) {
      // Timestamp and delta T
      double delta_t = (meas_package.timestamp_ - time_us_)/1000000.0;
      time_us_ = meas_package.timestamp_;
      std::cout << "delta T = " << delta_t << endl;

      // Prediction
      UKF::Prediction(delta_t);

      UKF::UpdateLidar(meas_package);
    } else if (use_radar_ && (meas_package.sensor_type_ == MeasurementPackage::RADAR)) {
      // Timestamp and delta T
      double delta_t = (meas_package.timestamp_ - time_us_)/1000000.0;
      time_us_ = meas_package.timestamp_;
      std::cout << "delta T = " << delta_t << endl;

      // Prediction
      UKF::Prediction(delta_t);

      UKF::UpdateRadar(meas_package);
    } else {
      std::cerr << "Unrecognized type of measurement data!" << endl;
    }
  }

}

