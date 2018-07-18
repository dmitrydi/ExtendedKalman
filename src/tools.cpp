#include <iostream>
#include "tools.h"


using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
    * Calculate the RMSE here.
  */
  long data_size = estimations.size();
  if (ground_truth.size() != data_size || data_size == 0) {
  	throw invalid_argument("number of estimations is uneaqual to number of groundtruth or zero length data");
  }

	VectorXd rmse(4);
	rmse << 0,0,0,0;
	for (long i = 0, i < data_size, i++) {
		VectorXd d = estimations[i] - groundtruth[i];
		d = d.array()*d.array();
		rmse += d;
  }
  rmse /= data_size;
  rmse = rmse.array().sqrt();
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
    * Calculate a Jacobian here.
  */
  const float TINY = numeric_limits<float>::min(); //small number to avoid division by zero
  float px, py, vx, vy;
  px = x_state(0);
  py = x_state(1);
  vx = x_state(2);
  vy = x_state(3);
  MatrixXd Jac;
  Jac = MatrixXd(3,4);
  float den = pow(px,2) + pow(py,2) + TINY;
  float sqrt_den = sqrt(den);

  Jac << px/sqrt_den, py/sqrt_den, 0, 0,
  			-py/den, px/den, 0, 0,
  			py*(vx*py - vy*px)/(den*sqrt_den), px*(vy*px - vx*py)/(den*sqrt_den), px/sqrt_den, py/sqrt_den;
  return Jac;
}

VectroXd Tools::Calculate_h(const VectroXd& x_state) {
	 /**
    * Calculate radar measurement function h(x) here.
  */
	const double pi = 3.1415926535897932384626433832795
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);
	float rho = sqrt(px*px+py*py);
	if (rho == 0.) {
		VectorXd h(0.,0.,0.);
	} else {
			if(px == 0.) {
				float phi = 0.5*pi;
			} else {
				phi = atan(py/px);
			}
			float rhodot = (px*vx+py*vy)/rho;
			VectorXd h(rho, phi, rhodot);
	}
	return h;
}