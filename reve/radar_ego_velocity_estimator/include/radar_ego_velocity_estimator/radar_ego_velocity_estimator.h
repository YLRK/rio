// This file is part of REVE - Radar Ego Velocity Estimator
// Copyright (C) 2021  Christopher Doer <christopher.doer@kit.edu>

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

#pragma once

#include <sensor_msgs/PointCloud2.h>

#include <pcl/PCLPointCloud2.h>
#include <pcl/pcl_macros.h>
#include <pcl/point_types.h>
#include <pcl_conversions/pcl_conversions.h>

#include <pcl_ros/transforms.h>

#include <radar_ego_velocity_estimator/data_types.h>
#include <radar_ego_velocity_estimator/radar_point_cloud.h>
#include <radar_ego_velocity_estimator/ros_helper.h>

#include <radar_ego_velocity_estimator/RadarEgoVelocityEstimatorConfig.h>

namespace reve
{
struct RadarEgoVelocityEstimatorIndices
{
  uint azimuth   = 0;
  uint elevation = 1;
  uint x_r       = 2;
  uint y_r       = 3;
  uint z_r       = 4;
  uint peak_db   = 5;
  uint r_x       = 6;
  uint r_y       = 7;
  uint r_z       = 8;
  uint v_d       = 9;
  uint noise_db  = 10;
};

class RadarEgoVelocityEstimator
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  /**
   * @brief RadarEgoVelocityEstimator constructor
   */
  RadarEgoVelocityEstimator() {}

  /**
   * @brief Reconfigure callback
   * @param config  has to contain RadarEgoVelocityEstimatorConfig
   * @return
   */
  template <class ConfigContainingRadarEgoVelocityEstimatorConfig>
  bool configure(ConfigContainingRadarEgoVelocityEstimatorConfig& config);

  /**
   * @brief Estimates the radar ego velocity based on a single radar scan
   * @param[in] radar_scan_msg       radar scan
   * @param[out] v_r                 estimated radar ego velocity
   * @param[out] sigma_v_r           estimated sigmas of ego velocity
   * @param[out] inlier_radar_scan   inlier point cloud
   * @returns true if estimation successful
   */
  bool estimate(const sensor_msgs::PointCloud2& radar_scan_msg, Vector3& v_r, Matrix3& P_v_r);
  bool estimate(const sensor_msgs::PointCloud2& radar_scan_msg, Vector3& v_r, Vector3& sigma_v_r);
  bool estimate(const sensor_msgs::PointCloud2& radar_scan_msg,
                Vector3& v_r,
                Matrix3& P_v_r,
                sensor_msgs::PointCloud2& inlier_radar_msg);
  bool estimate(const sensor_msgs::PointCloud2& radar_scan_msg,
                Vector3& v_r,
                Vector3& sigma_v_r,
                sensor_msgs::PointCloud2& inlier_radar_msg);
  bool estimate(const sensor_msgs::PointCloud2& radar_scan_msg,
                Vector3& v_r,
                Matrix3& P_v_r,
                pcl::PointCloud<RadarPointCloudType>& inlier_radar,
                const Matrix3& C_stab_r = Matrix3::Identity());

private:
  /**
   * @brief Implementation of the ransac based estimation
   * @param[in] radar_data          matrix of parsed radar scan --> see
   * RadarEgoVelocityEstimatorIndices
   * @param[out] v_r                estimated radar ego velocity
   * @param[out] sigma_v_r          estimated sigmas of ego velocity
   * @param[out] inlier_idx_best    idices of inlier
   * @returns true if estimation successful
   */
  bool solve3DLsqRansac(const Matrix& radar_data, Vector3& v_r, Matrix3& P_v_r, std::vector<uint>& inlier_idx_best);

  /**
   * @brief Estimates the radar ego velocity using all mesurements provided in
   * radar_data
   * @param[in] radar_data          matrix of parsed radar scan --> see
   * RadarEgoVelocityEstimatorIndices
   * @param[out] v_r                estimated radar ego velocity
   * @param[out] sigma_v_r          estimated sigmas of ego velocity
   * @param estimate_sigma          if true sigma will be estimated as well
   * @returns true if estimation successful
   */
  bool solve3DLsq(const Matrix& radar_data, Vector3& v_r, Matrix3& P_v_r, bool estimate_sigma = true);

  bool solve3DOdr(const Matrix& radar_data, Vector3& v_r, Matrix3& P_v_r);
  /**
   * @brief Helper function which estiamtes the number of RANSAC iterations
   */
  void setRansacIter()
  {
    ransac_iter_ = uint((std::log(1.0 - config_.success_prob)) /
                        std::log(1.0 - std::pow(1.0 - config_.outlier_prob, config_.N_ransac_points)));
  }

  const std::string kPrefix = "[RadarEgoVelocityEstimator]: ";
  const RadarEgoVelocityEstimatorIndices idx_;

  radar_ego_velocity_estimation::RadarEgoVelocityEstimatorConfig config_;
  uint ransac_iter_ = 0;
};

template <class ConfigContainingRadarEgoVelocityEstimatorConfig>
bool RadarEgoVelocityEstimator::configure(ConfigContainingRadarEgoVelocityEstimatorConfig& config)
{
  config_.min_dist                           = config.min_dist;
  config_.max_dist                           = config.max_dist;
  config_.min_db                             = config.min_db;
  config_.elevation_thresh_deg               = config.elevation_thresh_deg;
  config_.azimuth_thresh_deg                 = config.azimuth_thresh_deg;
  config_.filter_min_z                       = config.filter_min_z;
  config_.filter_max_z                       = config.filter_max_z;
  config_.doppler_velocity_correction_factor = config.doppler_velocity_correction_factor;

  config_.thresh_zero_velocity       = config.thresh_zero_velocity;
  config_.allowed_outlier_percentage = config.allowed_outlier_percentage;
  config_.sigma_zero_velocity_x      = config.sigma_zero_velocity_x;
  config_.sigma_zero_velocity_y      = config.sigma_zero_velocity_y;
  config_.sigma_zero_velocity_z      = config.sigma_zero_velocity_z;

  config_.sigma_offset_radar_x = config.sigma_offset_radar_x;
  config_.sigma_offset_radar_y = config.sigma_offset_radar_y;
  config_.sigma_offset_radar_z = config.sigma_offset_radar_z;

  config_.max_sigma_x                    = config.max_sigma_x;
  config_.max_sigma_y                    = config.max_sigma_y;
  config_.max_sigma_z                    = config.max_sigma_z;
  config_.max_r_cond                     = config.max_r_cond;
  config_.use_cholesky_instead_of_bdcsvd = config.use_cholesky_instead_of_bdcsvd;

  config_.use_ransac      = config.use_ransac;
  config_.outlier_prob    = config.outlier_prob;
  config_.success_prob    = config.success_prob;
  config_.N_ransac_points = config.N_ransac_points;
  config_.inlier_thresh   = config.inlier_thresh;

  config_.use_odr   = config.use_odr;
  config_.sigma_v_d = config.sigma_v_d;
  config_.min_speed_odr = config.min_speed_odr;
  config_.model_noise_offset_deg = config.model_noise_offset_deg;
  config_.model_noise_scale_deg = config.model_noise_scale_deg;

  setRansacIter();

  ROS_INFO_STREAM(kPrefix << "Number of Ransac iterations: " << ransac_iter_);

  return true;
}
}  // namespace reve
