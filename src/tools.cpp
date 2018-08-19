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
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  if (estimations.size()==0)
    {
        cout<<"the estimation vector size should not be zero"<<endl;
        return rmse;
    }
  if (estimations.size() != ground_truth.size())
    {
        cout<<"the estimation vector size should equal ground truth vector size"<<endl;
        return rmse;
    }

  //accumulate squared residuals
  for(unsigned int i=0; i < estimations.size(); ++i){
    //coefficient-wise multiplication
     VectorXd residual = estimations[i] - ground_truth[i];
              residual = residual.array()*residual.array();
		 rmse += residual;
	}

  //calculate the mean
  rmse = rmse/estimations.size();
  //calculate the squared root
  rmse = rmse.array().sqrt();

  return rmse;
  
}
