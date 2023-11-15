// This file is part of RIO - Radar Inertial Odometry and Radar ego velocity estimation.
// Copyright (C) 2021  Christopher Doer <christopher.doer@kit.edu>

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

#include <boost/math/distributions/chi_squared.hpp>

#include <rio_utils/math_helper.h>

#include <ekf_rio/ekf_rio_filter.h>

using namespace rio;

// EKF滤波器初始化，状态、姿态、偏置、协方差矩阵
bool EkfRioFilter::init(const std::vector<ImuDataStamped>& imu_init_vec, const Real& baro_h0)
{
  x_error_ = Vector::Zero(error_idx_.base_state_length);

  nav_sol_.setPosition_n_b(init_struct_.p_0);
  nav_sol_.v_n_b = init_struct_.v_0;

  // init att from acc
  Vector3 a_accu(0, 0, 0), w_accu(0, 0, 0);
  for (const auto& imu_data : imu_init_vec)
  {
    a_accu += imu_data.a_b_ib;
    w_accu += imu_data.w_b_ib;
  }
  const Vector3 acc_mean = a_accu / imu_init_vec.size();
  const Vector3 w_mean   = w_accu / imu_init_vec.size();

  EulerAngles roll_pitch = math_helper::initFromAcc(acc_mean, init_struct_.gravity);

  ROS_INFO_STREAM(kStreamingPrefix << "Initialized attitude: " << roll_pitch.to_degrees().x() << "deg, "
                                   << roll_pitch.to_degrees().y() << "deg");

  nav_sol_.setEuler_n_b(EulerAngles(roll_pitch.roll(), roll_pitch.pitch(), init_struct_.yaw_0));

  bias_.acc  = init_struct_.b_a_0;
  bias_.gyro = init_struct_.omega_calibration ? w_mean + init_struct_.b_w_0 : init_struct_.b_w_0;
  bias_.alt  = baro_h0;

  T_b_r_.translation() = init_struct_.l_b_r_0;
  T_b_r_.linear()      = Matrix3(init_struct_.q_b_r_0);

  ROS_INFO_STREAM(kStreamingPrefix << "Initialized w_bias: " << EulerAngles(bias_.gyro).to_degrees().transpose());

  covariance_ = init_struct_.P_kk_0;

  strapdown_ = Strapdown(init_struct_.gravity);

  time_stamp_ = imu_init_vec.back().time_stamp;

  return true;
}

// 状态传播，更新状态估计和协方差矩阵
bool EkfRioFilter::propagate(const ImuDataStamped& imu)
{
  time_stamp_ = imu.time_stamp;

  const Vector3 a_b_ib_corrected = imu.a_b_ib - bias_.acc;
  const Vector3 w_b_ib_corrected = imu.w_b_ib - bias_.gyro;
  nav_sol_                       = strapdown_.propagate(nav_sol_, a_b_ib_corrected, w_b_ib_corrected, imu.dt);

  const Matrix Phi = getPhi(a_b_ib_corrected, imu.dt);
  const Matrix G   = system_noise_.getG(nav_sol_.getC_n_b(), error_idx_, covariance_.rows());

  // TODO: consider more efficient implementation!!
  // 更新协方差矩阵，通过状态传播矩阵Phi与系统噪声矩阵G
  covariance_ = Phi * covariance_ * Phi.transpose() + G * system_noise_.getQ(imu.dt) * G.transpose();

  return true;
}

// 计算状态传播矩阵，描述了如何根据先前的状态和输入数据来预测下一个状态
Matrix EkfRioFilter::getPhi(Vector3 a_b_ib, const Real& T) const
{
  Matrix F = Matrix::Zero(error_idx_.prob_noise_state_length, error_idx_.prob_noise_state_length);
  F.block(error_idx_.position, error_idx_.velocity, 3, 3)  = Matrix::Identity(3, 3);
  F.block(error_idx_.velocity, error_idx_.attitude, 3, 3)  = -1 * math_helper::skewVec(nav_sol_.getC_n_b() * a_b_ib);
  F.block(error_idx_.velocity, error_idx_.bias_acc, 3, 3)  = -nav_sol_.getC_n_b();
  F.block(error_idx_.attitude, error_idx_.bias_gyro, 3, 3) = -nav_sol_.getC_n_b();

  Matrix Phi = Matrix::Identity(covariance_.rows(), covariance_.cols());
  Phi.block(0, 0, F.rows(), F.cols()) += F * T;

  return Phi;
}

//EKF状态更新
bool EkfRioFilter::kfUpdate(const Vector& r, const Matrix& H, const Vector& R_diag)
{
  const Matrix R = R_diag.asDiagonal();
  const Matrix S = H * covariance_ * H.transpose() + R;

  Matrix S_inv = Matrix::Identity(S.rows(), S.cols());
  S.llt().solveInPlace(S_inv);

  const Matrix K       = covariance_ * H.transpose() * S_inv; //卡尔曼增益
  const Vector x_error = K * r; // 状态估计校准误差

  covariance_ = covariance_ - K * H * covariance_; //更新协方差矩阵
  correctNominalState(x_error); 

  return true;  // TODO add check for invalid result
}

//校正估计的状态以匹配测量结果
bool EkfRioFilter::correctNominalState(const Vector x_error)
{
  //导航解状态校正
  nav_sol_.setPosition_n_b(nav_sol_.getPosition_n_b() - x_error.segment(error_idx_.position, 3));
  nav_sol_.v_n_b -= x_error.segment(error_idx_.velocity, 3);
  nav_sol_.setQuaternion(getCorrectedQuaternion(x_error.segment(error_idx_.attitude, 3), nav_sol_.getQuaternion_n_b()));

  //陀螺仪加速度计偏置校正
  bias_.acc -= x_error.segment(error_idx_.bias_acc, 3);
  bias_.gyro -= x_error.segment(error_idx_.bias_gyro, 3);
  bias_.alt -= x_error(error_idx_.bias_alt);

  //传感器到机体变换校正
  T_b_r_.translation() = T_b_r_.translation() - x_error.segment(error_idx_.l_b_r, 3);
  T_b_r_.linear()      = getCorrectedQuaternion(x_error.segment(error_idx_.eul_b_r, 3), Quaternion(T_b_r_.linear()))
                        .normalized()
                        .toRotationMatrix();

  if (covariance_.rows() > error_idx_.base_state_length)
  {
    radar_clone_state_.nav_sol.setPosition_n_b(radar_clone_state_.nav_sol.getPosition_n_b() -
                                               x_error.segment(error_idx_.sc_position, 3));
    radar_clone_state_.nav_sol.v_n_b -= x_error.segment(error_idx_.sc_velocity, 3);

    radar_clone_state_.nav_sol.setQuaternion(getCorrectedQuaternion(x_error.segment(error_idx_.sc_attitude, 3),
                                                                    radar_clone_state_.nav_sol.getQuaternion_n_b()));

    radar_clone_state_.T_b_r.translation() =
        radar_clone_state_.T_b_r.translation() - x_error.segment(error_idx_.sc_l_b_r, 3);
    radar_clone_state_.T_b_r.linear() =
        getCorrectedQuaternion(x_error.segment(error_idx_.sc_eul_b_r, 3), Quaternion(radar_clone_state_.T_b_r.linear()))
            .normalized()
            .toRotationMatrix();
  }

  return true;
}

// 添加雷达状态克隆
bool EkfRioFilter::addRadarStateClone(const ros::Time& trigger_stamp)
{
  if (covariance_.rows() > error_idx_.base_state_length)
  {
    removeRadarStateClone();
    ROS_WARN_STREAM(kStreamingPrefix << "Last radar clone was not removed, removed now!");
  }

  Matrix J_                               = Matrix::Zero(18, error_idx_.base_state_length);
  J_.block(0, 0, 9, 9)                    = Matrix::Identity(9, 9);
  J_.block(9, error_idx_.bias_gyro, 3, 3) = Matrix::Identity(3, 3);
  J_.block(12, error_idx_.l_b_r, 6, 6)    = Matrix::Identity(6, 6);

  Matrix J = Matrix::Zero(error_idx_.base_state_length + error_idx_.radar_clone_length, error_idx_.base_state_length);
  J.block(0, 0, error_idx_.base_state_length, error_idx_.base_state_length) =
      Matrix::Identity(error_idx_.base_state_length, error_idx_.base_state_length);
  J.block(error_idx_.base_state_length, 0, J_.rows(), J_.cols()) = J_;

  covariance_ = J * covariance_ * J.transpose();

  radar_clone_state_.time_stamp          = time_stamp_;
  radar_clone_state_.trigger_time_stamp  = trigger_stamp;
  radar_clone_state_.nav_sol             = nav_sol_;
  radar_clone_state_.offset_gyro         = bias_.gyro;
  radar_clone_state_.T_b_r.linear()      = T_b_r_.linear();
  radar_clone_state_.T_b_r.translation() = T_b_r_.translation();

  return true;
}

// 删除雷达状态的克隆
bool EkfRioFilter::removeRadarStateClone()
{
  if (covariance_.rows() > error_idx_.base_state_length)
  {
    if (covariance_.rows() > error_idx_.base_state_length + error_idx_.radar_clone_length)
    {
      ROS_ERROR_STREAM(kStreamingPrefix << "Dims of covariance larger than one clone, THIS SHOULD NOT HAPPEN!");
    }
    radar_clone_state_.time_stamp         = ros::TIME_MIN;
    radar_clone_state_.trigger_time_stamp = ros::TIME_MIN;

    // workaround: first row gets randomly zeros if no tmp used
    const Matrix covariance_tmp = covariance_;
    covariance_ = covariance_tmp.topLeftCorner(error_idx_.base_state_length, error_idx_.base_state_length);

    return true;
  }

  ROS_WARN_STREAM(kStreamingPrefix << "No radar clone state to be removed!");

  return false;
}

// 更新高度计的数据
bool EkfRioFilter::updateAltimeter(const Real neg_rel_h, const Real& sigma)
{
  Matrix H = Matrix::Zero(1, getCovarianceMatrix().cols());
  Vector r(1);
  Vector R_diag(1);

  H(0, getErrorIdx().position + 2) = 1;
  H(0, getErrorIdx().bias_alt)     = 1;
  const Real h_filter              = getNavigationSolution().getPosition_n_b().z();
  const Real h_meas                = neg_rel_h - getBias().alt;
  r(0)                             = h_filter - h_meas;
  R_diag(0)                        = sigma * sigma;
  kfUpdate(r, H, R_diag);

  return true;
}

// 更新雷达自身速度数据
bool EkfRioFilter::updateRadarEgoVelocity(const Vector3 v_r,
                                          const Vector3 sigma_v_r,
                                          const Vector3 w,
                                          const Real outlier_rejection_thresh)
{
  Matrix H = Matrix::Zero(3, getCovarianceMatrix().cols()); //构建雷达自身速度的观测模型
  Vector r(3); //存储雷达自身速度的观测残差
  Vector R_diag(3); //存储雷达自身速度的观测噪声方差

  // 获取克隆状态中的矩阵和矢量
  const Matrix C_n_b   = getRadarCloneState().nav_sol.getC_n_b(); 
  const Matrix C_b_r   = getRadarCloneState().T_b_r.linear();
  const Vector3 l_b_br = getRadarCloneState().T_b_r.translation();

  //计算雷达自身速度的各个分量
  const Vector3 v_w   = math_helper::skewVec(w - getRadarCloneState().offset_gyro) * l_b_br; //计算雷达的速度在雷达坐标系下的投影
  const Vector3 v_n_b = getRadarCloneState().nav_sol.v_n_b; //
  const Vector3 v_b   = C_n_b.transpose() * v_n_b; //雷达坐标系中的速度

  //计算观测模型的各个分量
  const Matrix3 H_v     = C_b_r.transpose() * C_n_b.transpose(); //雷达自身速度与导航速度之间的关系
  const Matrix3 H_q     = C_b_r.transpose() * C_n_b.transpose() * math_helper::skewVec(v_n_b);
  const Matrix3 H_bg    = -C_b_r.transpose() * math_helper::skewVec(l_b_br);
  const Matrix3 H_l_b_r = C_b_r.transpose() * math_helper::skewVec(w);
  const Matrix3 H_q_b_r = C_b_r.transpose() * math_helper::skewVec(v_w + v_b);

  H.block(0, getErrorIdx().sc_velocity, 3, 3)  = H_v;
  H.block(0, getErrorIdx().sc_attitude, 3, 3)  = H_q;
  H.block(0, getErrorIdx().sc_bias_gyro, 3, 3) = H_bg;
  H.block(0, getErrorIdx().sc_l_b_r, 3, 3)     = H_l_b_r;
  H.block(0, getErrorIdx().sc_eul_b_r, 3, 3)   = H_q_b_r;

  //雷达自身速度的估计值
  const Vector3 v_r_filter = C_b_r.transpose() * (v_w + v_b);

  r      = v_r_filter - v_r; //速度残差
  R_diag = (sigma_v_r).array().square(); //观测噪声协方差

  // outlier rejection 异常值排除
  // 马哈拉诺比斯距离
  if (outlier_rejection_thresh > 0.001)
  {
    // 计算马哈拉诺比斯距离，用于观测是否偏离了预测，以进行异常值检测
    
    const Real gamma =
        r.transpose() * (H * getCovarianceMatrix() * H.transpose() + Matrix(R_diag.asDiagonal())).inverse() * r;
    boost::math::chi_squared chiSquaredDist(3.0);
    // std::cout << "r" <<r.transpose() << std::endl;
    // std::cout << "2:"<<(H * getCovarianceMatrix() * H.transpose() + Matrix(R_diag.asDiagonal())).inverse() << std::endl;
    // std::cout << "H:" << H * getCovarianceMatrix() * H.transpose() << std::endl;
    // std::cout <<  "3:" <<(H * getCovarianceMatrix() * H.transpose() + Matrix(R_diag.asDiagonal())).inverse() * r << "\n" << std::endl;

    // 使用卡方分布计算阈值
    const double gamma_thresh = boost::math::quantile(chiSquaredDist, 1 - outlier_rejection_thresh);
    
    if (gamma < gamma_thresh)
    {
      kfUpdate(r, H, R_diag);
    }
    else
    {
      ROS_INFO_STREAM(kStreamingPrefix << "Outlier radar");
      return false;
    }
  }
  else
    kfUpdate(r, H, R_diag);

  return true;
}

Quaternion EkfRioFilter::getCorrectedQuaternion(const Vector3& err_euler, const Quaternion& q) const
{
  return math_helper::quaternionMultiplicationHamilton(
      Quaternion(1, -0.5 * err_euler.x(), -0.5 * err_euler.y(), -0.5 * err_euler.z()), q);
}
