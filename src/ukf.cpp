#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

#define EPS 1e-6

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
  std_a_ = 1.0;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;
  
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
  lambda_ = 3 - n_x_;
  n_aug_ = n_x_ + 2;
  n_sig_ = 2 * n_aug_ + 1;
  R_radar_ = MatrixXd(3, 3);
  R_radar_ << std_radr_*std_radr_, 0, 0,
	  0, std_radphi_*std_radphi_, 0,
	  0, 0, std_radrd_*std_radrd_;
  R_lidar_ = MatrixXd(2, 2);
  R_lidar_ << std_laspx_*std_laspx_, 0,
	  0, std_laspy_*std_laspy_;
  Xsig_pred_ = MatrixXd(n_x_, n_sig_); //create matrix with predicted sigma points as columns
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
	// Initialize
	if (!is_initialized_) {
        P_ << 1, 0, 0, 0, 0,
              0, 1, 0, 0, 0,
              0, 0, 1, 0, 0,
              0, 0, 0, 1, 0,
              0, 0, 0, 0, 1; 
		if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
			double rho = meas_package.raw_measurements_[0]; // range
			double phi = meas_package.raw_measurements_[1]; // bearing
			double rho_dot = meas_package.raw_measurements_[2]; // range rate
			double px = rho * cos(phi);
			double py = rho * sin(phi);
			double vx = rho_dot * cos(phi);
			double vy = rho_dot * sin(phi);
			double v = sqrt(vx*vx + vy*vy);
			x_ << px, py, v, 0, 0;
		}
		else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
			x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
		}
		weights_ = VectorXd(n_sig_);
		weights_.fill(0.5 / (lambda_ + n_aug_));
		weights_(0) = lambda_ / (lambda_ + n_aug_);
		time_us_ = meas_package.timestamp_;
		is_initialized_ = true;
		return;
	}

	double dt = (meas_package.timestamp_ - time_us_) / 1.0e6; // seconds
	time_us_ = meas_package.timestamp_;

	// Predict
	Prediction(dt);

	// Measurement updates
	if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
		UpdateRadar(meas_package);
	}
	else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
		UpdateLidar(meas_package);
	}
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
	//*** Step #1: Generate Sigma Points ***//
	VectorXd x_aug = VectorXd(n_aug_); // create augmented mean vector, 7
	MatrixXd P_aug = MatrixXd(n_aug_, n_aug_); // create augmented state covariance, 7x7
	MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sig_); // create sigma point matrix, 7x15
	
	//create augmented mean state
	x_aug.head(5) = x_;
	x_aug(5) = 0;
	x_aug(6) = 0;

	//create augmented covariance matrix
	P_aug.fill(0.0);
	P_aug.topLeftCorner(5, 5) = P_;
	P_aug(5, 5) = std_a_*std_a_;
	P_aug(6, 6) = std_yawdd_*std_yawdd_;

	//create square root matrix
	MatrixXd L = P_aug.llt().matrixL();

	//create augmented sigma points
	Xsig_aug.col(0) = x_aug;
	for (int i = 0; i< n_aug_; i++)
	{
		Xsig_aug.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
		Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
	}

	//*** Step #2: Predict Sigma Points ***//
	double dt2_half = 0.5 * delta_t * delta_t;
	for (int i = 0; i < n_sig_; ++i) {
		double v = Xsig_aug(2, i);
		double rho = Xsig_aug(3, i);
		double rho_dot = Xsig_aug(4, i);
		double nu_acc = Xsig_aug(5, i);
		double nu_rho_dd = Xsig_aug(6, i);
		if (fabs(rho_dot) > EPS) {
			Xsig_pred_(0, i) = v / rho_dot * (sin(rho + rho_dot*delta_t) - sin(rho));
			Xsig_pred_(1, i) = v / rho_dot * (-cos(rho + rho_dot*delta_t) + cos(rho));
		}
		else {
			Xsig_pred_(0, i) = v * cos(rho) * delta_t;
			Xsig_pred_(1, i) = v * sin(rho) * delta_t;
		}
		Xsig_pred_(0, i) += dt2_half*cos(rho)*nu_acc;
		Xsig_pred_(1, i) += dt2_half*sin(rho)*nu_acc;
		Xsig_pred_(2, i) = delta_t * nu_acc;
		Xsig_pred_(3, i) = rho_dot * delta_t + dt2_half * nu_rho_dd;
		Xsig_pred_(4, i) = delta_t * nu_rho_dd;
	}
	Xsig_pred_ += Xsig_aug.topLeftCorner(n_x_, 2 * n_aug_ + 1);

	//*** Step #3: Predict Mean and Covariance ***//
	x_.fill(0.0);
	P_.fill(0.0);
	
	for (int i = 0; i < n_sig_; ++i) {
		x_ += weights_(i) * Xsig_pred_.col(i);
	}

	for (int i = 0; i < n_sig_; ++i) {
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		NormalizeAngle(x_diff(3)); // angle normalization
		P_ += weights_(i) * x_diff * x_diff.transpose();
	}
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
	//*** Step #4: Predict Measurement ***//
	int n_z = 2;
	MatrixXd Zsig = MatrixXd(n_z, n_sig_); //matrix for sigma points in measurement space, 2x15
	VectorXd z_pred = VectorXd(n_z); //mean predicted measurement
	MatrixXd S = MatrixXd(n_z, n_z); //measurement covariance matrix S

	//transform sigma points into measurement space
	for (int i = 0; i < n_sig_; ++i) {
		double px = Xsig_pred_(0, i);
		double py = Xsig_pred_(1, i);
		Zsig(0, i) = px;
		Zsig(1, i) = py;
	}

	//calculate mean predicted measurement
	z_pred.fill(0.0);
	for (int i = 0; i < n_sig_; ++i) {
		z_pred += weights_(i) * Zsig.col(i);
	}

	//calculate innovation covariance matrix S
	S.fill(0.0);
	for (int i = 0; i < n_sig_; ++i) {
		VectorXd z_diff = Zsig.col(i) - z_pred;
		S += weights_(i) * z_diff * z_diff.transpose();
	}

	S += R_lidar_;

	//*** Step #5: Update States ***//
	MatrixXd Tc = MatrixXd(n_x_, n_z); //create matrix for cross correlation Tc

									   //calculate cross correlation matrix								   
	Tc.fill(0.0);
	for (int i = 0; i < n_sig_; ++i) {
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		NormalizeAngle(x_diff(3)); // angle normalization

		VectorXd z_diff = Zsig.col(i) - z_pred;
		Tc += weights_(i) * x_diff * z_diff.transpose();
	}

	//calculate Kalman gain K;
	MatrixXd K = MatrixXd(n_x_, n_z);
	K = Tc * S.inverse();

	//update state mean and covariance matrix
	VectorXd z = meas_package.raw_measurements_;
	VectorXd z_diff = z - z_pred;
	NormalizeAngle(z_diff(1)); // angle normalization
	x_ += K * z_diff;
	P_ -= K * S * K.transpose();

	// calculate NIS
	NIS_lidar_ = z_diff.transpose() * S.inverse() * z_diff;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  //*** Step #4: Predict Measurement ***//
	int n_z = 3;
	MatrixXd Zsig = MatrixXd(n_z, n_sig_); //matrix for sigma points in measurement space, 3x15
	VectorXd z_pred = VectorXd(n_z); //mean predicted measurement
	MatrixXd S = MatrixXd(n_z, n_z); //measurement covariance matrix S

	//transform sigma points into measurement space
	for (int i = 0; i < n_sig_; ++i) {
		double px = Xsig_pred_(0, i);
		double py = Xsig_pred_(1, i);
		double v = Xsig_pred_(2, i);
		double rho = Xsig_pred_(3, i);
		double px2py2_sq = sqrt(px*px + py*py);
		Zsig(0, i) = px2py2_sq;
		Zsig(1, i) = atan2(py, px);
		if (fabs(px2py2_sq) > EPS) {
			Zsig(2, i) = (px * cos(rho) * v + py * sin(rho) * v) / px2py2_sq;
		}
		else {
			Zsig(2, i) = 0.0;
		}
	}

	//calculate mean predicted measurement
	z_pred.fill(0.0);
	for (int i = 0; i < n_sig_; ++i) {
		z_pred += weights_(i) * Zsig.col(i);
	}

	//calculate innovation covariance matrix S
	S.fill(0.0);
	for (int i = 0; i < n_sig_; ++i) {
		VectorXd z_diff = Zsig.col(i) - z_pred;
		NormalizeAngle(z_diff(1)); // angle normalization
		S += weights_(i) * z_diff * z_diff.transpose();
	}

	S += R_radar_;

  //*** Step #5: Update States ***//
	MatrixXd Tc = MatrixXd(n_x_, n_z); //create matrix for cross correlation Tc
	
	//calculate cross correlation matrix								   
	Tc.fill(0.0);
	for (int i = 0; i < n_sig_; ++i) {
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		NormalizeAngle(x_diff(3)); // angle normalization

		VectorXd z_diff = Zsig.col(i) - z_pred;
		NormalizeAngle(z_diff(1)); // angle normalization
		Tc += weights_(i) * x_diff * z_diff.transpose();
	}

	//calculate Kalman gain K;
	MatrixXd K = MatrixXd(n_x_, n_z);
	K = Tc * S.inverse();

	//update state mean and covariance matrix
	VectorXd z = meas_package.raw_measurements_;
	VectorXd z_diff = z - z_pred;
	NormalizeAngle(z_diff(1)); // angle normalization
	x_ += K * z_diff;
	P_ -= K * S * K.transpose();

	// calculate NIS
	NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;
}

void UKF::NormalizeAngle(double &angle) {
	while (angle > M_PI) { angle -= M_PI; }
	while (angle < -M_PI) { angle += M_PI; }
}
