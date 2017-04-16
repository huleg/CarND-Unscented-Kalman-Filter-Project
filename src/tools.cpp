#include <iostream>
#include "tools.h"
#include <assert.h>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
    * Calculate the RMSE here.
    *
    * RMSE = SQRT(Sigma((estimation - ground_truth)^2))/num_samples
    *
  */
  // Data sizes
  int est_size = estimations.size();
  int gr_size = ground_truth.size();
  assert (est_size == gr_size);
  std::cout << "Total number of estimations: " << est_size << "\n";

  // RMSE
  VectorXd rmse(4);
  rmse << 0,0,0,0;
  for(int i=0; i< est_size; i++) {
    VectorXd tmp;
    // Element-wise operations
    tmp = estimations[i] - ground_truth[i];
    tmp = (tmp.array()*tmp.array());
    // Accumulation into rmse
    rmse += tmp;
  }
  rmse = rmse.array()/est_size;
  return rmse.array().sqrt();
}

