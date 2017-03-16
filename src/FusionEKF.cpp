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

  //previous_timestamp_ = 0;
  

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
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    // initial transition matrix F_
    ekf_.F_ = MatrixXd(4,4);
    ekf_.F_ << 1,0,1,0,
               0,1,0,1,
               0,0,1,0,
               0,0,0,1;

    // initial state covariance matrix P (error in the estimate/prediction)
    ekf_.P_ = MatrixXd (4,4);
    ekf_.P_ << 1,0,0,0,
               0,1,0,0,
               0,0,1000,0,
               0,0,0,1000;

      
    // Laser measurement matrix
    H_laser_ << 1,0,0,0,
                0,1,0,0;

    // Radar measurement matrix
    Hj_ << 1,1,0,0,
           1,1,0,0,
           1,1,1,1;

    // Use sigma_ax = 9 and sigma_ay = 9 for your Q matrix 
    sigma_ax = 9;
    sigma_ay = 9;

    // initial timestamp 
    previous_timestamp_ = measurement_pack.timestamp_;

  

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float px = measurement_pack.raw_measurements_[0]*cos(measurement_pack.raw_measurements_[1]);
      //px = rho * cos(phi)
      float py = measurement_pack.raw_measurements_[0]*sin(measurement_pack.raw_measurements_[1]);
      //py = rho * sin(phi)
      float vx = measurement_pack.raw_measurements_[2]*cos(measurement_pack.raw_measurements_[1]);
      //vx = rho_dot * cos(phi)
      float vy = measurement_pack.raw_measurements_[2]*sin(measurement_pack.raw_measurements_[1]);
      //vy = rho_dot * sin(phi)


      if (measurement_pack.raw_measurements_[0] == 0 or measurement_pack.raw_measurements_[1] == 0){
        ekf_.x_<< 0, 0, 0, 0;
        return;
      }
      
      ekf_.x_ <<px, py, vx, vy;
       
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      px = measurement_pack.raw_measurements_[0];
      py = measurement_pack.raw_measurements_[1];

      if (px == 0 or py == 0){
        return;
      }
      ekf_.x_<<px, py, 0, 0;
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use sigma_ax = 9 and sigma_ay = 9 for your Q matrix.
   */
  //dt - convert delta time to seconds
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0; 
  previous_timestamp_ = measurement_pack.timestamp_;
  cout << dt << endl;
  
  //update the transision Matrix F
  ekf_.F_ << 1, 0, dt, 0,
             0, 1, 0, dt,
             0, 0, 1, 0,
             0, 0, 0, 1;

  float dt_2 = pow(dt, 2);
  float dt_3 = pow(dt, 3);
  float dt_4 = pow(dt, 4);
  
  //set the process noise covariance matrix Q
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ <<  (dt_4/4)*sigma_ax, 0, (dt_3/2)*sigma_ax, 0,
              0, (dt_4/4)*sigma_ay, 0, (dt_3/2)*sigma_ay,
              (dt_3/2)*sigma_ax, 0, dt_2*sigma_ax, 0,
              0, (dt_3/2)*sigma_ay, 0, dt_2*sigma_ay;


  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    Hj_ << tools.CalculateJacobian(ekf_.x_);
    ekf_.H_ = Hj_;
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);  
  } 
  else {
    // Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }
  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
