#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * TODO: Calculate the RMSE here.
   */
  VectorXd rmse(4);
  rmse << 0,0,0,0;
    
  //accumulate squared residuals
  for(unsigned int i=0; i < estimations.size(); ++i){
    VectorXd residuals = estimations[i] - ground_truth[i];
    residuals = residuals.array() * residuals.array();
    rmse += residuals;
  }

  rmse = rmse / estimations.size();
  rmse = sqrt(rmse.array());
  
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * TODO:
   * Calculate a Jacobian here.
   */
  MatrixXd Hj(3,4);
    
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  float sum = px * px + py * py;

  //check division by zero
  if(sum < .00001) {
    px += .001;
    py += .001;
    sum = px * px + py * py;
  }

  float sqrt_sum = sqrt(sum);

  Hj << px/sqrt_sum, py/sqrt_sum, 0, 0,
        -py/sum, px/sum, 0, 0,
        (py*(vx*py - vy*px))/(sum*sqrt_sum), (px*(vy*px - vx*py))/(sum*sqrt_sum), px/sqrt_sum, py/sqrt_sum;

  return Hj;
}
