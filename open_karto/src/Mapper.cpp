/*
 * Copyright 2010 SRI International
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <sstream>
#include <fstream>
#include <stdexcept>
#include <set>
#include <list>
#include <iterator>

#include <math.h>
#include <assert.h>

#include "open_karto/Mapper.h"

namespace karto
{

  // enable this for verbose debug information
  // #define KARTO_DEBUG

  #define MAX_VARIANCE            500.0
  #define DISTANCE_PENALTY_GAIN   0.2
  #define ANGLE_PENALTY_GAIN      0.2

  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////

  void MapperSensorManager::RegisterSensor(const Name& rSensorName)
  {
    if (GetScanManager(rSensorName) == NULL)
    {
      m_ScanManagers[rSensorName] = new ScanManager(m_RunningBufferMaximumSize, m_RunningBufferMaximumDistance);
    }
  }


  /**
   * Gets scan from given device with given ID
   * @param rSensorName
   * @param scanNum
   * @return localized range scan
   */
  LocalizedRangeScan* MapperSensorManager::GetScan(const Name& rSensorName, kt_int32s scanIndex)
  {
    ScanManager* pScanManager = GetScanManager(rSensorName);
    if (pScanManager != NULL)
    {
      return pScanManager->GetScans().at(scanIndex);
    }

    assert(false);
    return NULL;
  }

  /**
   * Adds scan to scan vector of device that recorded scan
   * @param pScan
   */
  void MapperSensorManager::AddScan(LocalizedRangeScan* pScan)
  {
    GetScanManager(pScan)->AddScan(pScan, m_NextScanId);
    m_Scans.push_back(pScan);
    m_NextScanId++;
  }

  /**
   * Gets all scans of all devices
   * @return all scans of all devices
   */
  LocalizedRangeScanVector MapperSensorManager::GetAllScans()
  {
    LocalizedRangeScanVector scans;

    forEach(ScanManagerMap, &m_ScanManagers)
    {
      LocalizedRangeScanVector& rScans = iter->second->GetScans();

      scans.insert(scans.end(), rScans.begin(), rScans.end());
    }

    return scans;
  }

  /**
   * Deletes all scan managers of all devices
   */
  void MapperSensorManager::Clear()
  {
//    SensorManager::Clear();

    forEach(ScanManagerMap, &m_ScanManagers)
    {
      delete iter->second;
    }

    m_ScanManagers.clear();
  }

  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////

  ScanMatcher::~ScanMatcher()
  {
    delete m_pCorrelationGrid;
    delete m_pSearchSpaceProbs;
    delete m_pGridLookup;
  }

  ScanMatcher* ScanMatcher::Create(Mapper* pMapper, kt_double searchSize, kt_double resolution,
                                   kt_double smearDeviation, kt_double rangeThreshold)
  {
    // invalid parameters
    if (resolution <= 0)
    {
      return NULL;
    }
    if (searchSize <= 0)
    {
      return NULL;
    }
    if (smearDeviation < 0)
    {
      return NULL;
    }
    if (rangeThreshold <= 0)
    {
      return NULL;
    }

    assert(math::DoubleEqual(math::Round(searchSize / resolution), (searchSize / resolution)));

    // calculate search space in grid coordinates
    kt_int32u searchSpaceSideSize = static_cast<kt_int32u>(math::Round(searchSize / resolution) + 1);

    // compute requisite size of correlation grid (pad grid so that scan points can't fall off the grid
    // if a scan is on the border of the search space)
    kt_int32u pointReadingMargin = static_cast<kt_int32u>(ceil(rangeThreshold / resolution));

    kt_int32s gridSize = searchSpaceSideSize + 2 * pointReadingMargin;

    // create correlation grid
    assert(gridSize % 2 == 1);
    CorrelationGrid* pCorrelationGrid = CorrelationGrid::CreateGrid(gridSize, gridSize, resolution, smearDeviation);

    // create search space probabilities
    Grid<kt_double>* pSearchSpaceProbs = Grid<kt_double>::CreateGrid(searchSpaceSideSize,
                                                                     searchSpaceSideSize, resolution);

    ScanMatcher* pScanMatcher = new ScanMatcher(pMapper);
    pScanMatcher->m_pCorrelationGrid = pCorrelationGrid;
    pScanMatcher->m_pSearchSpaceProbs = pSearchSpaceProbs;
    pScanMatcher->m_pGridLookup = new GridIndexLookup<kt_int8u>(pCorrelationGrid);

    return pScanMatcher;
  }

  /**
   * 基于相关方法的scan-to-map匹配
   * 1、用前一帧优化前后的位姿变换，来校正当前帧优化前的位姿，得到一个初始位姿
   * 2、提取局部map中当前帧也可以看到的点，对应栅格设置为占用
   *    剔除相对于当前帧属于物体背面的点，也就是近邻帧与当前帧在物体两面
   * 3、scan-to-map匹配
   *    1) 创建旋转角度-激光点旋转后相对于当前帧位置的偏移量
   *       创建索引表，计算当前帧激光点，在各个旋转角度变换下，新的激光点位置相对于当前帧位置（机器人的位置）的坐标偏移量
   *       不管当前帧在什么位置，只要旋转角度一定，那么坐标偏移量也是确定的。
   *       目的是后面在搜索空间中对当前帧位姿施加不同的x平移、y平移之后再施加一个旋转的操作，不用重新计算每个搜索位置处旋转某个角度之后对应的点坐标，只需要加上对应旋转角度下算好的坐标偏移量即可。
   *    2) 遍历pose搜索空间计算响应值     
   *       a. 当前帧激光点集合在某个pose变换下得到新的位置，统计在局部map栅格中这些新位置处占用的数量，占用越多响应值越高，表示匹配越好
   *       b. 惩罚一下搜索偏移量，越远响应值折扣越多（当前帧搜索距离、角度偏移）
   *    3) 找到响应值最大的pose，如果有多个最佳响应pose，那么对pose求个均值，得到最终pose
   *    4) 计算位姿的协方差
   * 4、如果响应值是0，没匹配上，扩大角度搜索范围，增加20°、40°、60°，再执行一次scan-to-map匹配
   * 5、在scan-to-map优化之后的位姿基础上，再缩小搜索空间优化一次，搜索范围减半，角度区间减半，执行scan-to-map匹配
   * @param pScan         当前帧
   * @param rBaseScans    局部map（当前帧时空维度上相邻的帧集合）
   * @param rMean         输出优化后位姿
   * @param rCovariance   输出协方差矩阵，dx，dy，dtheta
   * @param doPenalize    是否对响应值做搜索距离上的惩罚
   * @param doRefineMatch 是否细化搜索空间，执行第二次优化
   * @return              返回响应值
   */
  kt_double ScanMatcher::MatchScan(LocalizedRangeScan* pScan, const LocalizedRangeScanVector& rBaseScans, Pose2& rMean,
                                   Matrix3& rCovariance, kt_bool doPenalize, kt_bool doRefineMatch)
  {
    // 用前一帧优化前后的位姿变换，来校正当前帧优化前的位姿，得到一个初始位姿
    Pose2 scanPose = pScan->GetSensorPose();

    // 当前帧激光点数为零
    if (pScan->GetNumberOfRangeReadings() == 0)
    {
      rMean = scanPose;

      // 协方差越大，位姿的不确定度越大
      rCovariance(0, 0) = MAX_VARIANCE;  // XX
      rCovariance(1, 1) = MAX_VARIANCE;  // YY
      rCovariance(2, 2) = 4 * math::Square(m_pMapper->m_pCoarseAngleResolution->GetValue());  // TH*TH

      return 0.0;
    }

    // 2. get size of grid
    Rectangle2<kt_int32s> roi = m_pCorrelationGrid->GetROI();

    // 3. compute offset (in meters - lower left corner)
    Vector2<kt_double> offset;
    offset.SetX(scanPose.GetX() - (0.5 * (roi.GetWidth() - 1) * m_pCorrelationGrid->GetResolution()));
    offset.SetY(scanPose.GetY() - (0.5 * (roi.GetHeight() - 1) * m_pCorrelationGrid->GetResolution()));

    // 4. set offset
    m_pCorrelationGrid->GetCoordinateConverter()->SetOffset(offset);

    ///////////////////////////////////////

    // 提取局部map中当前帧也可以看到的点，对应栅格设置为占用
    AddScans(rBaseScans, scanPose.GetPosition());

    // 搜索矩形框
    Vector2<kt_double> searchDimensions(m_pSearchSpaceProbs->GetWidth(), m_pSearchSpaceProbs->GetHeight());
    Vector2<kt_double> coarseSearchOffset(0.5 * (searchDimensions.GetX() - 1) * m_pCorrelationGrid->GetResolution(),
                                          0.5 * (searchDimensions.GetY() - 1) * m_pCorrelationGrid->GetResolution());

    // 搜索步长，2个栅格
    Vector2<kt_double> coarseSearchResolution(2 * m_pCorrelationGrid->GetResolution(),
                                              2 * m_pCorrelationGrid->GetResolution());

    // 基于相关方法的scan-to-map匹配
    // 1、创建旋转角度-激光点旋转后相对于当前帧位置的偏移量
    //    创建索引表，计算当前帧激光点，在各个旋转角度变换下，新的激光点位置相对于当前帧位置（机器人的位置）的坐标偏移量
    //    不管当前帧在什么位置，只要旋转角度一定，那么坐标偏移量也是确定的。
    //    目的是后面在搜索空间中对当前帧位姿施加不同的x平移、y平移之后再施加一个旋转的操作，不用重新计算每个搜索位置处旋转某个角度之后对应的点坐标，只需要加上对应旋转角度下算好的坐标偏移量即可。
    // 2、遍历pose搜索空间计算响应值
    //    1) 当前帧激光点集合在某个pose变换下得到新的位置，统计在局部map栅格中这些新位置处占用的数量，占用越多响应值越高，表示匹配越好
    //    2) 惩罚一下搜索偏移量，越远响应值折扣越多（当前帧搜索距离、角度偏移）
    // 3、找到响应值最大的pose，如果有多个最佳响应pose，那么对pose求个均值，得到最终pose
    // 4、计算位姿的协方差
    kt_double bestResponse = CorrelateScan(pScan, scanPose, coarseSearchOffset, coarseSearchResolution,
                                           m_pMapper->m_pCoarseSearchAngleOffset->GetValue(),
                                           m_pMapper->m_pCoarseAngleResolution->GetValue(),
                                           doPenalize, rMean, rCovariance, false);

    if (m_pMapper->m_pUseResponseExpansion->GetValue() == true)
    {
      // 如果响应值是0，没匹配上，扩大角度搜索范围
      if (math::DoubleEqual(bestResponse, 0.0))
      {
#ifdef KARTO_DEBUG
        std::cout << "Mapper Info: Expanding response search space!" << std::endl;
#endif
        // try and increase search angle offset with 20 degrees and do another match
        // 扩大搜索角度范围，增加20°、40°、60°
        kt_double newSearchAngleOffset = m_pMapper->m_pCoarseSearchAngleOffset->GetValue();
        for (kt_int32u i = 0; i < 3; i++)
        {
          newSearchAngleOffset += math::DegreesToRadians(20);
          // 再执行一次scan-to-map匹配
          bestResponse = CorrelateScan(pScan, scanPose, coarseSearchOffset, coarseSearchResolution,
                                       newSearchAngleOffset, m_pMapper->m_pCoarseAngleResolution->GetValue(),
                                       doPenalize, rMean, rCovariance, false);

          // 响应值还是0，退出
          if (math::DoubleEqual(bestResponse, 0.0) == false)
          {
            break;
          }
        }

#ifdef KARTO_DEBUG
        if (math::DoubleEqual(bestResponse, 0.0))
        {
          std::cout << "Mapper Warning: Unable to calculate response!" << std::endl;
        }
#endif
      }
    }

    // 在scan-to-map优化之后的位姿基础上，再缩小搜索空间优化一次
    if (doRefineMatch)
    {
      // 搜索范围减半，角度区间减半
      Vector2<kt_double> fineSearchOffset(coarseSearchResolution * 0.5);
      Vector2<kt_double> fineSearchResolution(m_pCorrelationGrid->GetResolution(), m_pCorrelationGrid->GetResolution());
      // 执行一次scan-to-map匹配
      bestResponse = CorrelateScan(pScan, rMean, fineSearchOffset, fineSearchResolution,
                                   0.5 * m_pMapper->m_pCoarseAngleResolution->GetValue(),
                                   m_pMapper->m_pFineSearchAngleOffset->GetValue(),
                                   doPenalize, rMean, rCovariance, true);
    }

#ifdef KARTO_DEBUG
    std::cout << "  BEST POSE = " << rMean << " BEST RESPONSE = " << bestResponse << ",  VARIANCE = "
              << rCovariance(0, 0) << ", " << rCovariance(1, 1) << std::endl;
#endif
    assert(math::InRange(rMean.GetHeading(), -KT_PI, KT_PI));

    return bestResponse;
  }

  /**
   * Finds the best pose for the scan centering the search in the correlation grid
   * at the given pose and search in the space by the vector and angular offsets
   * in increments of the given resolutions
   * 基于相关方法的scan-to-map匹配
   * 1、创建旋转角度-激光点旋转后相对于当前帧位置的偏移量
   *    创建索引表，计算当前帧激光点，在各个旋转角度变换下，新的激光点位置相对于当前帧位置（机器人的位置）的坐标偏移量
   *    不管当前帧在什么位置，只要旋转角度一定，那么坐标偏移量也是确定的。
   *    目的是后面在搜索空间中对当前帧位姿施加不同的x平移、y平移之后再施加一个旋转的操作，不用重新计算每个搜索位置处旋转某个角度之后对应的点坐标，只需要加上对应旋转角度下算好的坐标偏移量即可。
   * 2、遍历pose搜索空间计算响应值
   *    1) 当前帧激光点集合在某个pose变换下得到新的位置，统计在局部map栅格中这些新位置处占用的数量，占用越多响应值越高，表示匹配越好
   *    2) 惩罚一下搜索偏移量，越远响应值折扣越多（当前帧搜索距离、角度偏移）
   * 3、找到响应值最大的pose，如果有多个最佳响应pose，那么对pose求个均值，得到最终pose
   * 4、计算位姿的协方差
   * @param pScan                   当前帧
   * @param rSearchCenter           当前帧位姿（搜索空间的中心）
   * @param rSearchSpaceOffset      搜索范围，矩形半径
   * @param rSearchSpaceResolution  搜索步长
   * @param searchAngleOffset       搜索角度范围
   * @param searchAngleResolution   角度步长
   * @param doPenalize              惩罚搜索距离较远的位姿，降低响应值
   * @param rMean                   输出匹配结果位姿
   * @param rCovariance             输出匹配结果方差
   * @param doingFineMatch          是否执行二次优化
   * @return                        返回响应值
   */
  kt_double ScanMatcher::CorrelateScan(LocalizedRangeScan* pScan, const Pose2& rSearchCenter,
                                       const Vector2<kt_double>& rSearchSpaceOffset,
                                       const Vector2<kt_double>& rSearchSpaceResolution,
                                       kt_double searchAngleOffset, kt_double searchAngleResolution,
                                       kt_bool doPenalize, Pose2& rMean, Matrix3& rCovariance, kt_bool doingFineMatch)
  {
    assert(searchAngleResolution != 0.0);

    // 创建索引表，计算当前帧激光点，在各个旋转角度变换下，新的激光点位置相对于当前帧位置（机器人的位置）的坐标偏移量
    // 不管当前帧在什么位置，只要旋转角度一定，那么坐标偏移量也是确定的。
    // 目的是后面在搜索空间中对当前帧位姿施加不同的x平移、y平移之后再施加一个旋转的操作，不用重新计算每个搜索位置处旋转某个角度之后对应的点坐标，只需要加上对应旋转角度下算好的坐标偏移量即可。
    m_pGridLookup->ComputeOffsets(pScan, rSearchCenter.GetHeading(), searchAngleOffset, searchAngleResolution);

    // only initialize probability grid if computing positional covariance (during coarse match)
    if (!doingFineMatch)
    {
      m_pSearchSpaceProbs->Clear();

      // position search grid - finds lower left corner of search grid
      Vector2<kt_double> offset(rSearchCenter.GetPosition() - rSearchSpaceOffset);
      m_pSearchSpaceProbs->GetCoordinateConverter()->SetOffset(offset);
    }

    // calculate position arrays

    std::vector<kt_double> xPoses;
    // x轴搜索步数 
    kt_int32u nX = static_cast<kt_int32u>(math::Round(rSearchSpaceOffset.GetX() *
                                          2.0 / rSearchSpaceResolution.GetX()) + 1);
    // x轴搜索起始偏移
    kt_double startX = -rSearchSpaceOffset.GetX();
    // x轴搜索偏移集合
    for (kt_int32u xIndex = 0; xIndex < nX; xIndex++)
    {
      xPoses.push_back(startX + xIndex * rSearchSpaceResolution.GetX());
    }
    assert(math::DoubleEqual(xPoses.back(), -startX));

    std::vector<kt_double> yPoses;
    // y轴搜索步数
    kt_int32u nY = static_cast<kt_int32u>(math::Round(rSearchSpaceOffset.GetY() *
                                          2.0 / rSearchSpaceResolution.GetY()) + 1);
    // y轴搜索起始偏移
    kt_double startY = -rSearchSpaceOffset.GetY();
    // y轴搜索偏移集合
    for (kt_int32u yIndex = 0; yIndex < nY; yIndex++)
    {
      yPoses.push_back(startY + yIndex * rSearchSpaceResolution.GetY());
    }
    assert(math::DoubleEqual(yPoses.back(), -startY));

    // calculate pose response array size
    // 搜索角度数量
    kt_int32u nAngles = static_cast<kt_int32u>(math::Round(searchAngleOffset * 2.0 / searchAngleResolution) + 1);
    // 搜索空间数
    kt_int32u poseResponseSize = static_cast<kt_int32u>(xPoses.size() * yPoses.size() * nAngles);

    // allocate array
    // 搜索空间
    std::pair<kt_double, Pose2>* pPoseResponse = new std::pair<kt_double, Pose2>[poseResponseSize];
    // 搜索起点
    Vector2<kt_int32s> startGridPoint = m_pCorrelationGrid->WorldToGrid(Vector2<kt_double>(rSearchCenter.GetX()
                                                                        + startX, rSearchCenter.GetY() + startY));

    // 按照y方向、x方向、角度，遍历搜索空间
    kt_int32u poseResponseCounter = 0;
    forEachAs(std::vector<kt_double>, &yPoses, yIter)
    {
      // y偏移量
      kt_double y = *yIter;
      // 当前帧偏移y的位置
      kt_double newPositionY = rSearchCenter.GetY() + y;
      // 搜索距离平方
      kt_double squareY = math::Square(y);

      forEachAs(std::vector<kt_double>, &xPoses, xIter)
      {
        // x偏移量
        kt_double x = *xIter;
        // 当前帧偏移x的位置
        kt_double newPositionX = rSearchCenter.GetX() + x;
        // 搜索距离平方
        kt_double squareX = math::Square(x);

        // 当前帧偏移(x,y)之后所在的网格点坐标，网格索引
        Vector2<kt_int32s> gridPoint = m_pCorrelationGrid->WorldToGrid(Vector2<kt_double>(newPositionX, newPositionY));
        kt_int32s gridIndex = m_pCorrelationGrid->GridIndex(gridPoint);
        assert(gridIndex >= 0);

        kt_double angle = 0.0;
        // 起始角度，heading是0°朝向，searchAngleOffset是一半搜索角度偏移量
        kt_double startAngle = rSearchCenter.GetHeading() - searchAngleOffset;
        for (kt_int32u angleIndex = 0; angleIndex < nAngles; angleIndex++)
        {
          // 角度
          angle = startAngle + angleIndex * searchAngleResolution;

          // 计算搜索pose的响应值
          // 当前帧激光点集合在某个pose变换下得到新的位置，统计在局部map栅格中这些新位置处占用的数量，占用越多响应值越高，表示匹配越好
          kt_double response = GetResponse(angleIndex, gridIndex);
          // 惩罚一下搜索偏移量，越远响应值折扣越多（当前帧搜索距离、角度偏移）
          if (doPenalize && (math::DoubleEqual(response, 0.0) == false))
          {
            // simple model (approximate Gaussian) to take odometry into account
            // 距离平方
            kt_double squaredDistance = squareX + squareY;
            // 距离惩罚项
            kt_double distancePenalty = 1.0 - (DISTANCE_PENALTY_GAIN *
                                               squaredDistance / m_pMapper->m_pDistanceVariancePenalty->GetValue());
            distancePenalty = math::Maximum(distancePenalty, m_pMapper->m_pMinimumDistancePenalty->GetValue());

            // 角度惩罚项
            kt_double squaredAngleDistance = math::Square(angle - rSearchCenter.GetHeading());
            kt_double anglePenalty = 1.0 - (ANGLE_PENALTY_GAIN *
                                            squaredAngleDistance / m_pMapper->m_pAngleVariancePenalty->GetValue());
            anglePenalty = math::Maximum(anglePenalty, m_pMapper->m_pMinimumAnglePenalty->GetValue());

            // 响应打个折扣，距离、角度偏移越大，响应越小
            // 意思是理论上纠正的位姿与当前的位姿不会变化很大，如果特别大就不太靠谱
            response *= (distancePenalty * anglePenalty);
          }

          // store response and pose
          // 存响应值、新pose
          pPoseResponse[poseResponseCounter] = std::pair<kt_double, Pose2>(response, Pose2(newPositionX, newPositionY,
                                                                           math::NormalizeAngle(angle)));
          poseResponseCounter++;
        }

        assert(math::DoubleEqual(angle, rSearchCenter.GetHeading() + searchAngleOffset));
      }
    }

    assert(poseResponseSize == poseResponseCounter);

    // find value of best response (in [0; 1])
    // 找到响应值最大的pose
    kt_double bestResponse = -1;
    for (kt_int32u i = 0; i < poseResponseSize; i++)
    {
      bestResponse = math::Maximum(bestResponse, pPoseResponse[i].first);

      // will compute positional covariance, save best relative probability for each cell
      if (!doingFineMatch)
      {
        const Pose2& rPose = pPoseResponse[i].second;
        Vector2<kt_int32s> grid = m_pSearchSpaceProbs->WorldToGrid(rPose.GetPosition());

        // Changed (kt_double*) to the reinterpret_cast - Luc
        kt_double* ptr = reinterpret_cast<kt_double*>(m_pSearchSpaceProbs->GetDataPointer(grid));
        if (ptr == NULL)
        {
          throw std::runtime_error("Mapper FATAL ERROR - Index out of range in probability search!");
        }

        *ptr = math::Maximum(pPoseResponse[i].first, *ptr);
      }
    }

    // average all poses with same highest response
    // 有多个最佳响应pose，那么对pose求个均值，得到最终pose
    Vector2<kt_double> averagePosition;
    kt_double thetaX = 0.0;
    kt_double thetaY = 0.0;
    kt_int32s averagePoseCount = 0;
    for (kt_int32u i = 0; i < poseResponseSize; i++)
    {
      if (math::DoubleEqual(pPoseResponse[i].first, bestResponse))
      {
        averagePosition += pPoseResponse[i].second.GetPosition();

        kt_double heading = pPoseResponse[i].second.GetHeading();
        thetaX += cos(heading);
        thetaY += sin(heading);

        averagePoseCount++;
      }
    }

    Pose2 averagePose;
    if (averagePoseCount > 0)
    {
      averagePosition /= averagePoseCount;

      thetaX /= averagePoseCount;
      thetaY /= averagePoseCount;

      averagePose = Pose2(averagePosition, atan2(thetaY, thetaX));
    }
    else
    {
      throw std::runtime_error("Mapper FATAL ERROR - Unable to find best position");
    }

    // delete pose response array
    delete [] pPoseResponse;

#ifdef KARTO_DEBUG
    std::cout << "bestPose: " << averagePose << std::endl;
    std::cout << "bestResponse: " << bestResponse << std::endl;
#endif

    if (!doingFineMatch)
    {
      // 计算新的位姿的协方差，位置
      ComputePositionalCovariance(averagePose, bestResponse, rSearchCenter, rSearchSpaceOffset,
                                  rSearchSpaceResolution, searchAngleResolution, rCovariance);
    }
    else
    {
      // 计算新的位姿的协方差，角度
      ComputeAngularCovariance(averagePose, bestResponse, rSearchCenter,
                              searchAngleOffset, searchAngleResolution, rCovariance);
    }

    // 最终pose
    rMean = averagePose;

#ifdef KARTO_DEBUG
    std::cout << "bestPose: " << averagePose << std::endl;
#endif

    if (bestResponse > 1.0)
    {
      bestResponse = 1.0;
    }

    assert(math::InRange(bestResponse, 0.0, 1.0));
    assert(math::InRange(rMean.GetHeading(), -KT_PI, KT_PI));

    // 最佳响应
    return bestResponse;
  }

  /**
   * 最大响应值只有一个（或者一部分），还有另外一部分不是最大响应值，但与最大值很接近，把这些响应值对应的pose加在一起，
   * 计算协方差。协方差的变量是dx,dy，什么意思呢，dx = newPose - oldPose，也就是优化后pose与优化前pose的位置差，
   * 这个值作为协方差的变量。根据方差公式D = (E - x)^2，E就是均值，这里用最佳pose-oldPose作为均值，x就是前面说到的
   * 这些pose-oldPose。当然这个变量还可以乘上一个reponse，累加计算均值，最后再把reponse除掉，等于是计算过程中考虑了
   * 响应这一权重。
   * 其实是不是可以理解为优化后的所有pose，是否聚集在最优pose周围，也就是方差大小
   * @param rBestPose               最佳候选pose
   * @param bestResponse            最佳响应值
   * @param rSearchCenter           搜索中心
   * @param rSearchSpaceOffset      搜索范围偏移
   * @param rSearchSpaceResolution  搜索步长
   * @param searchAngleResolution   搜索角度步长
   * @param rCovariance             协方差，x,y,theta
   */
  void ScanMatcher::ComputePositionalCovariance(const Pose2& rBestPose, kt_double bestResponse,
                                                const Pose2& rSearchCenter,
                                                const Vector2<kt_double>& rSearchSpaceOffset,
                                                const Vector2<kt_double>& rSearchSpaceResolution,
                                                kt_double searchAngleResolution, Matrix3& rCovariance)
  {
    // reset covariance to identity matrix
    rCovariance.SetToIdentity();

    // if best response is vary small return max variance
    // 响应值太小，协方差设置max
    if (bestResponse < KT_TOLERANCE)
    {
      rCovariance(0, 0) = MAX_VARIANCE;  // XX
      rCovariance(1, 1) = MAX_VARIANCE;  // YY
      rCovariance(2, 2) = 4 * math::Square(searchAngleResolution);  // TH*TH

      return;
    }

    kt_double accumulatedVarianceXX = 0;
    kt_double accumulatedVarianceXY = 0;
    kt_double accumulatedVarianceYY = 0;
    kt_double norm = 0;

    // 位置偏差
    kt_double dx = rBestPose.GetX() - rSearchCenter.GetX();
    kt_double dy = rBestPose.GetY() - rSearchCenter.GetY();

    // 搜索范围偏移
    kt_double offsetX = rSearchSpaceOffset.GetX();
    kt_double offsetY = rSearchSpaceOffset.GetY();

    kt_int32u nX = static_cast<kt_int32u>(math::Round(offsetX * 2.0 / rSearchSpaceResolution.GetX()) + 1);
    kt_double startX = -offsetX;
    assert(math::DoubleEqual(startX + (nX - 1) * rSearchSpaceResolution.GetX(), -startX));

    kt_int32u nY = static_cast<kt_int32u>(math::Round(offsetY * 2.0 / rSearchSpaceResolution.GetY()) + 1);
    kt_double startY = -offsetY;
    assert(math::DoubleEqual(startY + (nY - 1) * rSearchSpaceResolution.GetY(), -startY));
     // 按照搜索步长划分搜索空间，遍历
    for (kt_int32u yIndex = 0; yIndex < nY; yIndex++)
    {
      kt_double y = startY + yIndex * rSearchSpaceResolution.GetY();

      for (kt_int32u xIndex = 0; xIndex < nX; xIndex++)
      {
        kt_double x = startX + xIndex * rSearchSpaceResolution.GetX();

        Vector2<kt_int32s> gridPoint = m_pSearchSpaceProbs->WorldToGrid(Vector2<kt_double>(rSearchCenter.GetX() + x,
                                                                                           rSearchCenter.GetY() + y));
        // 该位置处响应
        kt_double response = *(m_pSearchSpaceProbs->GetDataPointer(gridPoint));

        // response is not a low response
        // 响应与最佳响应值接近
        if (response >= (bestResponse - 0.1))
        {
          // dx、dy相当于均值，优化后pose与优化前pose距离之差，x等同
          // 计算方差
          norm += response;
          accumulatedVarianceXX += (math::Square(x - dx) * response);
          accumulatedVarianceXY += ((x - dx) * (y - dy) * response);
          accumulatedVarianceYY += (math::Square(y - dy) * response);
        }
      }
    }

    if (norm > KT_TOLERANCE)
    {
      // 计算平均方差
      kt_double varianceXX = accumulatedVarianceXX / norm;
      kt_double varianceXY = accumulatedVarianceXY / norm;
      kt_double varianceYY = accumulatedVarianceYY / norm;
      kt_double varianceTHTH = 4 * math::Square(searchAngleResolution);

      // lower-bound variances so that they are not too small;
      // ensures that links are not too tight
      kt_double minVarianceXX = 0.1 * math::Square(rSearchSpaceResolution.GetX());
      kt_double minVarianceYY = 0.1 * math::Square(rSearchSpaceResolution.GetY());
      varianceXX = math::Maximum(varianceXX, minVarianceXX);
      varianceYY = math::Maximum(varianceYY, minVarianceYY);

      // increase variance for poorer responses
      kt_double multiplier = 1.0 / bestResponse;
      rCovariance(0, 0) = varianceXX * multiplier;
      rCovariance(0, 1) = varianceXY * multiplier;
      rCovariance(1, 0) = varianceXY * multiplier;
      rCovariance(1, 1) = varianceYY * multiplier;
      rCovariance(2, 2) = varianceTHTH;  // this value will be set in ComputeAngularCovariance
    }

    // if values are 0, set to MAX_VARIANCE
    // values might be 0 if points are too sparse and thus don't hit other points
    if (math::DoubleEqual(rCovariance(0, 0), 0.0))
    {
      rCovariance(0, 0) = MAX_VARIANCE;
    }

    if (math::DoubleEqual(rCovariance(1, 1), 0.0))
    {
      rCovariance(1, 1) = MAX_VARIANCE;
    }
  }

  /**
   * 计算角度的协方差
   * @param rBestPose
   * @param bestResponse
   * @param rSearchCenter
   * @param rSearchAngleOffset
   * @param searchAngleResolution
   * @param rCovariance
   */
  void ScanMatcher::ComputeAngularCovariance(const Pose2& rBestPose,
                                             kt_double bestResponse,
                                             const Pose2& rSearchCenter,
                                             kt_double searchAngleOffset,
                                             kt_double searchAngleResolution,
                                             Matrix3& rCovariance)
  {
    // NOTE: do not reset covariance matrix

    // normalize angle difference
    kt_double bestAngle = math::NormalizeAngleDifference(rBestPose.GetHeading(), rSearchCenter.GetHeading());

    Vector2<kt_int32s> gridPoint = m_pCorrelationGrid->WorldToGrid(rBestPose.GetPosition());
    kt_int32s gridIndex = m_pCorrelationGrid->GridIndex(gridPoint);

    kt_int32u nAngles = static_cast<kt_int32u>(math::Round(searchAngleOffset * 2 / searchAngleResolution) + 1);

    kt_double angle = 0.0;
    kt_double startAngle = rSearchCenter.GetHeading() - searchAngleOffset;

    kt_double norm = 0.0;
    kt_double accumulatedVarianceThTh = 0.0;
    for (kt_int32u angleIndex = 0; angleIndex < nAngles; angleIndex++)
    {
      angle = startAngle + angleIndex * searchAngleResolution;
      kt_double response = GetResponse(angleIndex, gridIndex);

      // response is not a low response
      if (response >= (bestResponse - 0.1))
      {
        norm += response;
        accumulatedVarianceThTh += (math::Square(angle - bestAngle) * response);
      }
    }
    assert(math::DoubleEqual(angle, rSearchCenter.GetHeading() + searchAngleOffset));

    if (norm > KT_TOLERANCE)
    {
      if (accumulatedVarianceThTh < KT_TOLERANCE)
      {
        accumulatedVarianceThTh = math::Square(searchAngleResolution);
      }

      accumulatedVarianceThTh /= norm;
    }
    else
    {
      accumulatedVarianceThTh = 1000 * math::Square(searchAngleResolution);
    }

    rCovariance(2, 2) = accumulatedVarianceThTh;
  }

  /**
   * 提取局部map中当前帧也可以看到的点，对应栅格设置为占用
   * @param rScans    局部map（当前帧时空维度上相邻的帧集合）
   * @param viewPoint 当前帧位置
   */
  void ScanMatcher::AddScans(const LocalizedRangeScanVector& rScans, Vector2<kt_double> viewPoint)
  {
    m_pCorrelationGrid->Clear();

    // 在局部近邻帧激光点中集合中提取近邻帧和当前帧都可以看到的点（剔除相对于当前帧属于物体背面的点，也就是近邻帧与当前帧在物体两面），对应的栅格设置占用
    const_forEach(LocalizedRangeScanVector, &rScans)
    {
      AddScan(*iter, viewPoint);
    }
  }

  /**
   * 在局部近邻帧激光点中集合中提取近邻帧和当前帧都可以看到的点（剔除相对于当前帧属于物体背面的点，也就是近邻帧与当前帧在物体两面），对应的栅格设置占用
   * @param pScan      局部近邻帧
   * @param rViewPoint 当前帧位置
   */
  void ScanMatcher::AddScan(LocalizedRangeScan* pScan, const Vector2<kt_double>& rViewPoint, kt_bool doSmear)
  {
    // 在局部近邻帧激光点中集合中提取近邻帧和当前帧都可以看到的点
    // 例如近邻帧看到一堵墙的一面，下面代码判断出来当前帧在这堵墙的另一面，那么近邻帧看到的这些墙面点，当前帧都是看不到的，墙是有厚度的
    PointVectorDouble validPoints = FindValidPoints(pScan, rViewPoint);

    // 遍历局部近邻帧有效点，对应栅格（局部map对应的栅格）设置为占用，值为GridStates_Occupied=100，未占用默认为0
    const_forEach(PointVectorDouble, &validPoints)
    {
      Vector2<kt_int32s> gridPoint = m_pCorrelationGrid->WorldToGrid(*iter);
      if (!math::IsUpTo(gridPoint.GetX(), m_pCorrelationGrid->GetROI().GetWidth()) ||
          !math::IsUpTo(gridPoint.GetY(), m_pCorrelationGrid->GetROI().GetHeight()))
      {
        // point not in grid
        continue;
      }

      int gridIndex = m_pCorrelationGrid->GridIndex(gridPoint);

      // set grid cell as occupied
      if (m_pCorrelationGrid->GetDataPointer()[gridIndex] == GridStates_Occupied)
      {
        // value already set
        continue;
      }

      m_pCorrelationGrid->GetDataPointer()[gridIndex] = GridStates_Occupied;

      // smear grid
      if (doSmear == true)
      {
        m_pCorrelationGrid->SmearPoint(gridPoint);
      }
    }
  }

  /**
   * 在局部近邻帧激光点中集合中提取近邻帧和当前帧都可以看到的点
   * 例如近邻帧看到一堵墙的一面，下面代码判断出来当前帧在这堵墙的另一面，那么近邻帧看到的这些墙面点，当前帧都是看不到的，墙是有厚度的
   * @param pScan       局部近邻帧
   * @param rViewPoint  当前帧位置
   * @return            返回局部近邻帧和当前帧可以同时看到的那些点
   */
  PointVectorDouble ScanMatcher::FindValidPoints(LocalizedRangeScan* pScan, const Vector2<kt_double>& rViewPoint) const
  {
    // 局部近邻帧的激光点
    const PointVectorDouble& rPointReadings = pScan->GetPointReadings();

    // points must be at least 10 cm away when making comparisons of inside/outside of viewpoint
    const kt_double minSquareDistance = math::Square(0.1);  // in m^2

    // this iterator lags from the main iterator adding points only when the points are on
    // the same side as the viewpoint
    PointVectorDouble::const_iterator trailingPointIter = rPointReadings.begin();
    PointVectorDouble validPoints;

    Vector2<kt_double> firstPoint;
    kt_bool firstTime = true;
    // 遍历局部近邻帧的激光点
    const_forEach(PointVectorDouble, &rPointReadings)
    {
      Vector2<kt_double> currentPoint = *iter;

      if (firstTime && !std::isnan(currentPoint.GetX()) && !std::isnan(currentPoint.GetY()))
      {
        firstPoint = currentPoint;
        firstTime = false;
      }

      // 与前一个检查激光点的距离，要求至少大于10cm
      // 距离大于10cm的时候，计算两个点相对于当前帧位置是顺时针还是逆时针的顺序，如果是逆时针，那么两个点之间的所有点都会加进来，否则都不加
      Vector2<kt_double> delta = firstPoint - currentPoint;
      if (delta.SquaredLength() > minSquareDistance)
      {
        // This compute the Determinant (viewPoint FirstPoint, viewPoint currentPoint)
        // Which computes the direction of rotation, if the rotation is counterclock
        // wise then we are looking at data we should keep. If it's negative rotation
        // we should not included in in the matching
        // have enough distance, check viewpoint
        // 平面向量的叉乘，a=(x1,y1),b=(x2,y2),axb=x1y2-x2y1,如果大于0，b在a的逆时针方向，如果小于0，b在a的顺时针方向
        // 例子：比如一堵墙，近邻帧在墙的一边，两个点按时序是逆时针的顺序，如果以当前帧位置为参考，两个点是顺时针的顺序，那么当前帧在墙的另一边
        // 显然这些点当前帧是看不到的，不需要加进来
        double a = rViewPoint.GetY() - firstPoint.GetY();
        double b = firstPoint.GetX() - rViewPoint.GetX();
        double c = firstPoint.GetY() * rViewPoint.GetX() - firstPoint.GetX() * rViewPoint.GetY();
        double ss = currentPoint.GetX() * a + currentPoint.GetY() * b + c;

        // reset beginning point
        firstPoint = currentPoint;

        if (ss < 0.0)  // wrong side, skip and keep going
        {
          trailingPointIter = iter;
        }
        else
        {
          for (; trailingPointIter != iter; ++trailingPointIter)
          {
            validPoints.push_back(*trailingPointIter);
          }
        }
      }
    }

    return validPoints;
  }

  /**
   * Get response at given position for given rotation (only look up valid points)
   * 计算搜索pose的响应值
   * 当前帧激光点集合在某个pose变换下得到新的位置，统计在局部map栅格中这些新位置处占用的数量，占用越多响应值越高，表示匹配越好
   * @param angleIndex        旋转角度索引
   * @param gridPositionIndex 新的激光帧网格位置索引
   * @return                  响应值
   */
  kt_double ScanMatcher::GetResponse(kt_int32u angleIndex, kt_int32s gridPositionIndex) const
  {
    kt_double response = 0.0;

    // add up value for each point
    // 局部map的栅格数据，存的是100占用，默认值0非占用，加上gridPositionIndex之后，pByte指向当前新的激光帧网格位置处
    kt_int8u* pByte = m_pCorrelationGrid->GetDataPointer() + gridPositionIndex;

    // 索引表内容为旋转角-旋转后激光点相对帧位置的坐标偏移量（转网格索引）之间的对应关系，取出当前旋转角下对应的偏移数据
    const LookupArray* pOffsets = m_pGridLookup->GetLookupArray(angleIndex);
    assert(pOffsets != NULL);

    // get number of points in offset list
    // 激光点数量
    kt_int32u nPoints = pOffsets->GetSize();
    if (nPoints == 0)
    {
      return response;
    }

    // calculate response
    kt_int32s* pAngleIndexPointer = pOffsets->GetArrayPointer();
    // 遍历激光点
    for (kt_int32u i = 0; i < nPoints; i++)
    {
      // ignore points that fall off the grid
      kt_int32s pointGridIndex = gridPositionIndex + pAngleIndexPointer[i];
      if (!math::IsUpTo(pointGridIndex, m_pCorrelationGrid->GetDataSize()) || pAngleIndexPointer[i] == INVALID_SCAN)
      {
        continue;
      }

      // uses index offsets to efficiently find location of point in the grid
      // pByte指向局部map栅格中新（经过平移）的激光帧位置，再加上pAngleIndexPointer[i]是当前旋转之后新的激光点相对于激光帧位置的偏移量，
      // 就等于新的激光点的位置，取它的栅格值。
      // 也就是当前帧pose施加一个平移、旋转，激光点落在的栅格对应值相加，即为当前帧pose的响应值
      // 100占用，0非占用。注：虽然定义了0未知、100占用、255free，但这三个是用于更新地图的，这里计算响应只用到100、0
      response += pByte[pAngleIndexPointer[i]];
    }

    // 归一化[0,1]，响应值越大越好，表示scan-to-map匹配越好
    // 仅通过占用的数量就能判定匹配的好不好，这个在激光里面是合理的，因为激光打到物体表面上，占用点就是薄薄的一层，甚至1个像素，
    // free点就特别多了，能匹配到占用点，说明很接近了。
    // 如果点云是平面，那么匹配可以很好的纠正旋转，平移就不好说了；如果是墙角这种形状，旋转和平移都能纠正。
    // 类似于VO里面的直接法，两帧图像直接用像素灰度值配准，当然它有很强的灰度不变假设。
    // 上面仅是个人理解。
    response /= (nPoints * GridStates_Occupied);
    assert(fabs(response) <= 1.0);

    return response;
  }


  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////

  MapperGraph::MapperGraph(Mapper* pMapper, kt_double rangeThreshold)
    : m_pMapper(pMapper)
  {
    m_pLoopScanMatcher = ScanMatcher::Create(pMapper, m_pMapper->m_pLoopSearchSpaceDimension->GetValue(),
                                             m_pMapper->m_pLoopSearchSpaceResolution->GetValue(),
                                             m_pMapper->m_pLoopSearchSpaceSmearDeviation->GetValue(), rangeThreshold);
    assert(m_pLoopScanMatcher);

    m_pTraversal = new BreadthFirstTraversal<LocalizedRangeScan>(this);
  }

  MapperGraph::~MapperGraph()
  {
    delete m_pLoopScanMatcher;
    m_pLoopScanMatcher = NULL;

    delete m_pTraversal;
    m_pTraversal = NULL;
  }

  void MapperGraph::AddVertex(LocalizedRangeScan* pScan)
  {
    assert(pScan);

    if (pScan != NULL)
    {
      Vertex<LocalizedRangeScan>* pVertex = new Vertex<LocalizedRangeScan>(pScan);
      Graph<LocalizedRangeScan>::AddVertex(pScan->GetSensorName(), pVertex);
      if (m_pMapper->m_pScanOptimizer != NULL)
      {
        m_pMapper->m_pScanOptimizer->AddNode(pVertex);
      }
    }
  }

  /**
   * 构建BA图的边
   * 1、与前一帧建立一条边
   * 2、在局部map帧集合中选取与当前帧最接近的一帧，建立边
   * 3、在相邻帧的chain中取与当前帧最接近的一帧，建立边
   *    1) 在顶点图中广搜找当前帧的相邻帧，然后对每个相邻帧截取前后时间范围内、且距离小于阈值的一段帧集合（称为chain）
   *    2) 用当前帧与chain进行scan-to-map匹配
   *    3) 如果匹配较好，在chain里面选最近一帧，建立边
   *    4) 对所有数量足够的chain执行前面2、3步
   *    注：相当于当前帧跟比较近，但又没有直接关系的帧（时间上相邻帧，局部匹配帧，闭环帧）建立边
   * 如果闭环成功，也会添加一条边。（注：不在这个方法里面）
   * 4、最后设置当前帧pose为第3步中计算得到的优化pose集合（包含当前帧pose）的均值
   * @param pScan       当前帧
   * @param rCovariance 当前帧pose协方差
  */
  void MapperGraph::AddEdges(LocalizedRangeScan* pScan, const Matrix3& rCovariance)
  {
    MapperSensorManager* pSensorManager = m_pMapper->m_pMapperSensorManager;

    const Name& rSensorName = pScan->GetSensorName();

    // link to previous scan
    // 1.连接前一帧
    kt_int32s previousScanNum = pScan->GetStateId() - 1;
    if (pSensorManager->GetLastScan(rSensorName) != NULL)
    {
      assert(previousScanNum >= 0);
      // 两帧之间创建一条边，添加到图中（BA），边的信息包括：前一帧位姿，当前帧位姿，当前帧位姿协方差
      LinkScans(pSensorManager->GetScan(rSensorName, previousScanNum), pScan, pScan->GetSensorPose(), rCovariance);
    }

    Pose2Vector means;
    std::vector<Matrix3> covariances;

    // first scan (link to first scan of other robots)
    if (pSensorManager->GetLastScan(rSensorName) == NULL)
    {
      assert(pSensorManager->GetScans(rSensorName).size() == 1);

      std::vector<Name> deviceNames = pSensorManager->GetSensorNames();
      forEach(std::vector<Name>, &deviceNames)
      {
        const Name& rCandidateSensorName = *iter;

        // skip if candidate device is the same or other device has no scans
        if ((rCandidateSensorName == rSensorName) || (pSensorManager->GetScans(rCandidateSensorName).empty()))
        {
          continue;
        }

        Pose2 bestPose;
        Matrix3 covariance;
        kt_double response = m_pMapper->m_pSequentialScanMatcher->MatchScan(pScan,
                                                                  pSensorManager->GetScans(rCandidateSensorName),
                                                                  bestPose, covariance);
        LinkScans(pSensorManager->GetScan(rCandidateSensorName, 0), pScan, bestPose, covariance);

        // only add to means and covariances if response was high "enough"
        if (response > m_pMapper->m_pLinkMatchMinimumResponseFine->GetValue())
        {
          means.push_back(bestPose);
          covariances.push_back(covariance);
        }
      }
    }
    else
    {
      // link to running scans
      // 当前帧位姿
      Pose2 scanPose = pScan->GetSensorPose();
      means.push_back(scanPose);
      covariances.push_back(rCovariance);
      // 2.在局部map帧集合中选取与当前帧最接近的一帧，建立边
      LinkChainToScan(pSensorManager->GetRunningScans(rSensorName), pScan, scanPose, rCovariance);
    }

    // link to other near chains (chains that include new scan are invalid)
    // 3、在相邻帧的chain中取与当前帧最接近的一帧，建立边
    //    1) 在顶点图中广搜找当前帧的相邻帧，然后对每个相邻帧截取前后时间范围内、且距离小于阈值的一段帧集合（称为chain）
    //    2) 用当前帧与chain进行scan-to-map匹配
    //    3) 如果匹配较好，在chain里面选最近一帧，建立边
    //    4) 对所有数量足够的chain执行前面2、3步
    //    注：相当于当前帧跟比较近，但又没有直接关系的帧（时间上相邻帧，局部匹配帧，闭环帧）建立边
    LinkNearChains(pScan, means, covariances);

    if (!means.empty())
    {
      // 给定多个pose和对应的协方差矩阵，计算平均pose
      pScan->SetSensorPose(ComputeWeightedMean(means, covariances));
    }
  }

  /**
   * 闭环检测
   * 1、遍历所有历史帧，提取候选闭环帧集合
   *    1) 在已构造的顶点图中，广搜找与当前顶点帧临近的所有顶点帧，记为临近帧集合
   *    2) 按时序遍历所有历史帧，距离够近的帧加入候选闭环集合，一旦发现一帧在临近帧集合中，清空候选闭环集合，接着往下找
   *    3) 直到距离超过阈值，如果候选闭环集合数量够多，认为合格
   * 2、当前帧与候选闭环帧集合，执行scan-to-map匹配
   * 3、如果响应够大，方差够小，执行二次scan-to-map匹配
   * 4、如果确认闭环，在候选闭环帧集合中选取与当前帧最接近的一帧，建立边，同时更新当前帧位姿
   * 5、循环提取下一个闭环候选帧集合，寻找闭环
  */
  kt_bool MapperGraph::TryCloseLoop(LocalizedRangeScan* pScan, const Name& rSensorName)
  {
    kt_bool loopClosed = false;

    kt_int32u scanIndex = 0;
    // 遍历所有历史帧，提取候选闭环帧集合
    // 1、在已构造的顶点图中，广搜找与当前顶点帧临近的所有顶点帧，记为临近帧集合
    // 2、按时序遍历所有历史帧，距离够近的帧加入候选闭环集合，一旦发现一帧在临近帧集合中，清空候选闭环集合，接着往下找
    // 3、直到距离超过阈值，如果候选闭环集合数量够多，认为合格
    LocalizedRangeScanVector candidateChain = FindPossibleLoopClosure(pScan, rSensorName, scanIndex);

    while (!candidateChain.empty())
    {
      Pose2 bestPose;
      Matrix3 covariance;
      // 当前帧与候选闭环帧集合，执行scan-to-map匹配
      kt_double coarseResponse = m_pLoopScanMatcher->MatchScan(pScan, candidateChain,
                                                               bestPose, covariance, false, false);

      std::stringstream stream;
      stream << "COARSE RESPONSE: " << coarseResponse
             << " (> " << m_pMapper->m_pLoopMatchMinimumResponseCoarse->GetValue() << ")"
             << std::endl;
      stream << "            var: " << covariance(0, 0) << ",  " << covariance(1, 1)
             << " (< " << m_pMapper->m_pLoopMatchMaximumVarianceCoarse->GetValue() << ")";

      m_pMapper->FireLoopClosureCheck(stream.str());

      // 如果响应够大，方差够小，执行二次匹配
      if ((coarseResponse > m_pMapper->m_pLoopMatchMinimumResponseCoarse->GetValue()) &&
          (covariance(0, 0) < m_pMapper->m_pLoopMatchMaximumVarianceCoarse->GetValue()) &&
          (covariance(1, 1) < m_pMapper->m_pLoopMatchMaximumVarianceCoarse->GetValue()))
      {
        // 复制当前帧
        LocalizedRangeScan tmpScan(pScan->GetSensorName(), pScan->GetRangeReadingsVector());
        tmpScan.SetUniqueId(pScan->GetUniqueId());
        tmpScan.SetTime(pScan->GetTime());
        tmpScan.SetStateId(pScan->GetStateId());
        tmpScan.SetCorrectedPose(pScan->GetCorrectedPose());
        tmpScan.SetSensorPose(bestPose);  // This also updates OdometricPose.
        // 用新的位姿，再执行一次scan-to-map匹配
        kt_double fineResponse = m_pMapper->m_pSequentialScanMatcher->MatchScan(&tmpScan, candidateChain,
                                                                                bestPose, covariance, false);

        std::stringstream stream1;
        stream1 << "FINE RESPONSE: " << fineResponse << " (>"
                << m_pMapper->m_pLoopMatchMinimumResponseFine->GetValue() << ")" << std::endl;
        m_pMapper->FireLoopClosureCheck(stream1.str());

        // 响应值不够大，认为闭环不成立
        if (fineResponse < m_pMapper->m_pLoopMatchMinimumResponseFine->GetValue())
        {
          m_pMapper->FireLoopClosureCheck("REJECTED!");
        }
        else
        {
          // 确认闭环
          m_pMapper->FireBeginLoopClosure("Closing loop...");

          // 更新当前帧pose
          pScan->SetSensorPose(bestPose);
          // 在候选闭环帧集合中选取与当前帧最接近的一帧，建立边
          LinkChainToScan(candidateChain, pScan, bestPose, covariance);
          CorrectPoses();

          m_pMapper->FireEndLoopClosure("Loop closed!");

          loopClosed = true;
        }
      }

      // 用pose更新之后的当前帧，在剩下（scanIndex在累加）的历史帧中继续寻找闭环帧集合
      candidateChain = FindPossibleLoopClosure(pScan, rSensorName, scanIndex);
    }

    return loopClosed;
  }

  /**
   * 在局部map中，找到与当前帧pose最接近的一帧
   * @param rScans 局部map
   * @param rPose  当前帧pose
  */
  LocalizedRangeScan* MapperGraph::GetClosestScanToPose(const LocalizedRangeScanVector& rScans,
                                                        const Pose2& rPose) const
  {
    LocalizedRangeScan* pClosestScan = NULL;
    kt_double bestSquaredDistance = DBL_MAX;

    // 遍历局部map
    const_forEach(LocalizedRangeScanVector, &rScans)
    {
      // 局部帧的pose
      Pose2 scanPose = (*iter)->GetReferencePose(m_pMapper->m_pUseScanBarycenter->GetValue());
      
      // 与当前帧pose之间的距离
      // 前面的pose有两种表示，一种是点云中心位置，一种是激光帧pose。这里计算距离，也就可以是两种表示了
      kt_double squaredDistance = rPose.GetPosition().SquaredDistance(scanPose.GetPosition());
      if (squaredDistance < bestSquaredDistance)
      {
        bestSquaredDistance = squaredDistance;
        pClosestScan = *iter;
      }
    }

    return pClosestScan;
  }

  Edge<LocalizedRangeScan>* MapperGraph::AddEdge(LocalizedRangeScan* pSourceScan,
                                                 LocalizedRangeScan* pTargetScan, kt_bool& rIsNewEdge)
  {
    // check that vertex exists
    assert(pSourceScan->GetStateId() < (kt_int32s)m_Vertices[pSourceScan->GetSensorName()].size());
    assert(pTargetScan->GetStateId() < (kt_int32s)m_Vertices[pTargetScan->GetSensorName()].size());

    Vertex<LocalizedRangeScan>* v1 = m_Vertices[pSourceScan->GetSensorName()][pSourceScan->GetStateId()];
    Vertex<LocalizedRangeScan>* v2 = m_Vertices[pTargetScan->GetSensorName()][pTargetScan->GetStateId()];

    // see if edge already exists
    const_forEach(std::vector<Edge<LocalizedRangeScan>*>, &(v1->GetEdges()))
    {
      Edge<LocalizedRangeScan>* pEdge = *iter;

      if (pEdge->GetTarget() == v2)
      {
        rIsNewEdge = false;
        return pEdge;
      }
    }

    Edge<LocalizedRangeScan>* pEdge = new Edge<LocalizedRangeScan>(v1, v2);
    Graph<LocalizedRangeScan>::AddEdge(pEdge);
    rIsNewEdge = true;
    return pEdge;
  }

  /**
   * 两帧之间创建一条边，添加到图中（BA），边的信息包括：前一帧位姿，当前帧位姿，当前帧位姿协方差
   * @param pFromScan   前一帧
   * @param pToScan     当前帧
   * @param rMean       当前帧pose
   * @param rCovariance 当前帧pose协方差
  */
  void MapperGraph::LinkScans(LocalizedRangeScan* pFromScan, LocalizedRangeScan* pToScan,
                              const Pose2& rMean, const Matrix3& rCovariance)
  {
    kt_bool isNewEdge = true;
    Edge<LocalizedRangeScan>* pEdge = AddEdge(pFromScan, pToScan, isNewEdge);

    // only attach link information if the edge is new
    if (isNewEdge == true)
    {
      pEdge->SetLabel(new LinkInfo(pFromScan->GetSensorPose(), rMean, rCovariance));
      if (m_pMapper->m_pScanOptimizer != NULL)
      {
        m_pMapper->m_pScanOptimizer->AddConstraint(pEdge);
      }
    }
  }

  /**
   * 1、在顶点图中广搜找当前帧的相邻帧，然后对每个相邻帧截取前后时间范围内、且距离小于阈值的一段帧集合（称为chain）
   * 2、用当前帧与chain进行scan-to-map匹配
   * 3、如果匹配较好，在chain里面选最近一帧，建立边
   * 4、对所有数量足够的chain执行前面2、3步
   * 注：相当于当前帧跟比较近，但又没有直接关系的帧（时间上相邻帧，局部匹配帧，闭环帧）建立边
   * @param pScan         当前帧
   * @param rMeans        当前帧与各个临近帧的chain匹配，优化后的pose集合
   * @param rCovariances  当前帧与各个临近帧的chain匹配，优化后的pose对应协方差集合
  */
  void MapperGraph::LinkNearChains(LocalizedRangeScan* pScan, Pose2Vector& rMeans, std::vector<Matrix3>& rCovariances)
  {
    // 在顶点图中当前帧位置广搜得到相邻帧集合，然后对于每个相邻帧，取前后时间范围内（距离也要满足）的一部分帧作为该帧的chain，相当于时空维度上临近帧集合
    // 返回{chain1,chain2,chain3,...}，每个chain是当前帧的相邻帧对应的chain
    const std::vector<LocalizedRangeScanVector> nearChains = FindNearChains(pScan);
    const_forEach(std::vector<LocalizedRangeScanVector>, &nearChains)
    {
      // chain中激光帧的数量需要大于一个阈值
      if (iter->size() < m_pMapper->m_pLoopMatchMinimumChainSize->GetValue())
      {
        continue;
      }

      Pose2 mean;
      Matrix3 covariance;
      // match scan against "near" chain
      // 当前帧与chain进行scan-to-map匹配
      kt_double response = m_pMapper->m_pSequentialScanMatcher->MatchScan(pScan, *iter, mean, covariance, false);
      // 响应值够大，说明匹配程度比较高
      if (response > m_pMapper->m_pLinkMatchMinimumResponseFine->GetValue() - KT_TOLERANCE)
      {
        // 与chain匹配后的pose加入集合
        rMeans.push_back(mean);
        rCovariances.push_back(covariance);
        // 在chain集合中选取与当前帧最接近的一帧，建立边
        LinkChainToScan(*iter, pScan, mean, covariance);
      }
    }
  }

  /**
   * 在局部map（当前帧的局部map、候选闭环帧map、相邻帧的chain）帧集合中选取与当前帧最接近的一帧，建立边
   * @param rChain        局部map
   * @param pScan         当前帧scan
   * @param rMean         当前帧pose
   * @param rCovariance   当前帧pose协方差
  */
  void MapperGraph::LinkChainToScan(const LocalizedRangeScanVector& rChain, LocalizedRangeScan* pScan,
                                    const Pose2& rMean, const Matrix3& rCovariance)
  {
    // 当前帧pose有两种表示方式，一种是用点云中心位置表示，朝向是0°，一种是用激光帧自身pose表示
    Pose2 pose = pScan->GetReferencePose(m_pMapper->m_pUseScanBarycenter->GetValue());
    // 在局部map中，与当前帧pose最接近的一帧
    LocalizedRangeScan* pClosestScan = GetClosestScanToPose(rChain, pose);
    assert(pClosestScan != NULL);
    // 这一帧的pose
    Pose2 closestScanPose = pClosestScan->GetReferencePose(m_pMapper->m_pUseScanBarycenter->GetValue());
    // 计算与当前帧的pose距离
    kt_double squaredDistance = pose.GetPosition().SquaredDistance(closestScanPose.GetPosition());
    // 如果距离不算远,就添加一条边
    if (squaredDistance < math::Square(m_pMapper->m_pLinkScanMaximumDistance->GetValue()) + KT_TOLERANCE)
    {
      LinkScans(pClosestScan, pScan, rMean, rCovariance);
    }
  }

  /**
   * 在顶点图中当前帧位置广搜得到相邻帧集合，然后对于每个相邻帧，取前后时间范围内（距离也要满足）的一部分帧作为该帧的chain，相当于时空维度上临近帧集合
   * 返回{chain1,chain2,chain3,...}，每个chain是当前帧的相邻帧对应的chain
   * @param pScan 当前帧
  */
  std::vector<LocalizedRangeScanVector> MapperGraph::FindNearChains(LocalizedRangeScan* pScan)
  {
    std::vector<LocalizedRangeScanVector> nearChains;
    // 当前帧pose
    Pose2 scanPose = pScan->GetReferencePose(m_pMapper->m_pUseScanBarycenter->GetValue());

    // to keep track of which scans have been added to a chain
    LocalizedRangeScanVector processed;

    // 在已构造的顶点图中，广搜找与当前顶点帧临近的所有顶点帧
    const LocalizedRangeScanVector nearLinkedScans = FindNearLinkedScans(pScan,
                                                     m_pMapper->m_pLinkScanMaximumDistance->GetValue());
    // 遍历临近帧集合
    const_forEach(LocalizedRangeScanVector, &nearLinkedScans)
    {
      LocalizedRangeScan* pNearScan = *iter;

      if (pNearScan == pScan)
      {
        continue;
      }

      // scan has already been processed, skip
      if (find(processed.begin(), processed.end(), pNearScan) != processed.end())
      {
        continue;
      }

      processed.push_back(pNearScan);

      // build up chain
      kt_bool isValidChain = true;
      std::list<LocalizedRangeScan*> chain;

      // add scans before current scan being processed
      // 遍历临近帧时间之前的那些帧，选取一定距离范围内的帧加入chain集合
      for (kt_int32s candidateScanNum = pNearScan->GetStateId() - 1; candidateScanNum >= 0; candidateScanNum--)
      {
        LocalizedRangeScan* pCandidateScan = m_pMapper->m_pMapperSensorManager->GetScan(pNearScan->GetSensorName(),
                                                                                        candidateScanNum);

        // chain is invalid--contains scan being added
        if (pCandidateScan == pScan)
        {
          isValidChain = false;
        }
        // 计算与当前帧pose之间的距离
        Pose2 candidatePose = pCandidateScan->GetReferencePose(m_pMapper->m_pUseScanBarycenter->GetValue());
        kt_double squaredDistance = scanPose.GetPosition().SquaredDistance(candidatePose.GetPosition());

        // 如果小于阈值，加入chain集合
        if (squaredDistance < math::Square(m_pMapper->m_pLinkScanMaximumDistance->GetValue()) + KT_TOLERANCE)
        {
          chain.push_front(pCandidateScan);
          processed.push_back(pCandidateScan);
        }
        else
        {
          break;
        }
      }

      // 加入临近帧
      chain.push_back(pNearScan);

      // add scans after current scan being processed
      // 遍历临近帧之后的那些帧，选取一定距离范围内的帧加入chain集合
      kt_int32u end = static_cast<kt_int32u>(m_pMapper->m_pMapperSensorManager->GetScans(pNearScan->GetSensorName()).size());
      for (kt_int32u candidateScanNum = pNearScan->GetStateId() + 1; candidateScanNum < end; candidateScanNum++)
      {
        LocalizedRangeScan* pCandidateScan = m_pMapper->m_pMapperSensorManager->GetScan(pNearScan->GetSensorName(),
                                                                                        candidateScanNum);

        if (pCandidateScan == pScan)
        {
          isValidChain = false;
        }

        Pose2 candidatePose = pCandidateScan->GetReferencePose(m_pMapper->m_pUseScanBarycenter->GetValue());;
        kt_double squaredDistance = scanPose.GetPosition().SquaredDistance(candidatePose.GetPosition());

        if (squaredDistance < math::Square(m_pMapper->m_pLinkScanMaximumDistance->GetValue()) + KT_TOLERANCE)
        {
          chain.push_back(pCandidateScan);
          processed.push_back(pCandidateScan);
        }
        else
        {
          break;
        }
      }

      if (isValidChain)
      {
        // change list to vector
        LocalizedRangeScanVector tempChain;
        std::copy(chain.begin(), chain.end(), std::inserter(tempChain, tempChain.begin()));
        // add chain to collection
        nearChains.push_back(tempChain);
      }
    }

    return nearChains;
  }

  /**
   * 在已构造的顶点图中，广搜找与当前顶点帧临近的所有顶点帧
   * @param pScan       当前帧
   * @param maxDistance 最大距离 
  */
  LocalizedRangeScanVector MapperGraph::FindNearLinkedScans(LocalizedRangeScan* pScan, kt_double maxDistance)
  {
    NearScanVisitor* pVisitor = new NearScanVisitor(pScan, maxDistance, m_pMapper->m_pUseScanBarycenter->GetValue());
    // 广度优先遍历，找到与起始顶点pose距离小于阈值的所有顶点
    LocalizedRangeScanVector nearLinkedScans = m_pTraversal->Traverse(GetVertex(pScan), pVisitor);
    delete pVisitor;

    return nearLinkedScans;
  }

  /**
   * 给定多个pose和对应的协方差矩阵，计算平均pose
  */
  Pose2 MapperGraph::ComputeWeightedMean(const Pose2Vector& rMeans, const std::vector<Matrix3>& rCovariances) const
  {
    assert(rMeans.size() == rCovariances.size());

    // compute sum of inverses and create inverse list
    std::vector<Matrix3> inverses;
    inverses.reserve(rCovariances.size());

    Matrix3 sumOfInverses;
    const_forEach(std::vector<Matrix3>, &rCovariances)
    {
      Matrix3 inverse = iter->Inverse();
      inverses.push_back(inverse);

      sumOfInverses += inverse;
    }
    Matrix3 inverseOfSumOfInverses = sumOfInverses.Inverse();

    // compute weighted mean
    Pose2 accumulatedPose;
    kt_double thetaX = 0.0;
    kt_double thetaY = 0.0;

    Pose2Vector::const_iterator meansIter = rMeans.begin();
    const_forEach(std::vector<Matrix3>, &inverses)
    {
      Pose2 pose = *meansIter;
      kt_double angle = pose.GetHeading();
      thetaX += cos(angle);
      thetaY += sin(angle);

      Matrix3 weight = inverseOfSumOfInverses * (*iter);
      accumulatedPose += weight * pose;

      ++meansIter;
    }

    thetaX /= rMeans.size();
    thetaY /= rMeans.size();
    accumulatedPose.SetHeading(atan2(thetaY, thetaX));

    return accumulatedPose;
  }

  /**
   * 遍历所有历史帧，提取候选闭环帧集合
   * 1、在已构造的顶点图中，广搜找与当前顶点帧临近的所有顶点帧，记为临近帧集合
   * 2、按时序遍历所有历史帧，距离够近的帧加入候选闭环集合，一旦发现一帧在临近帧集合中，清空候选闭环集合，接着往下找
   * 3、直到距离超过阈值，如果候选闭环集合数量够多，认为合格
   * @param pScan       当前帧
   * @param rSensorName sensor名字
   * @param rStartNum   起始查找帧idx
  */
  LocalizedRangeScanVector MapperGraph::FindPossibleLoopClosure(LocalizedRangeScan* pScan,
                                                                const Name& rSensorName,
                                                                kt_int32u& rStartNum)
  {
    LocalizedRangeScanVector chain;  // return value
    // 当前帧pose
    Pose2 pose = pScan->GetReferencePose(m_pMapper->m_pUseScanBarycenter->GetValue());

    // possible loop closure chain should not include close scans that have a
    // path of links to the scan of interest
    // 在已构造的顶点图中，广搜找与当前顶点帧临近的所有顶点帧
    const LocalizedRangeScanVector nearLinkedScans =
          FindNearLinkedScans(pScan, m_pMapper->m_pLoopSearchMaximumDistance->GetValue());

    // 按时序遍历所有历史帧
    kt_int32u nScans = static_cast<kt_int32u>(m_pMapper->m_pMapperSensorManager->GetScans(rSensorName).size());
    for (; rStartNum < nScans; rStartNum++)
    {
      // 候选匹配帧
      LocalizedRangeScan* pCandidateScan = m_pMapper->m_pMapperSensorManager->GetScan(rSensorName, rStartNum);
      // 候选匹配帧pose
      Pose2 candidateScanPose = pCandidateScan->GetReferencePose(m_pMapper->m_pUseScanBarycenter->GetValue());

      // 候选匹配帧与当前帧pose差，距离小于阈值，添加可能的闭环帧
      kt_double squaredDistance = candidateScanPose.GetPosition().SquaredDistance(pose.GetPosition());
      if (squaredDistance < math::Square(m_pMapper->m_pLoopSearchMaximumDistance->GetValue()) + KT_TOLERANCE)
      {
        // a linked scan cannot be in the chain
        // 候选闭环帧不能出现在当前帧的临近帧集合中，闭环帧一定是与当前帧时间上相隔很远，距离比较接近的帧，临近帧在时间、空间上都临近
        // 一旦当前帧不是候选闭环帧，之前加入的候选帧都删掉。因为是按时序遍历的历史帧，如果当前帧不是候选闭环帧，之前加入的也不是
        if (find(nearLinkedScans.begin(), nearLinkedScans.end(), pCandidateScan) != nearLinkedScans.end())
        {
          chain.clear();
        }
        else
        {
          chain.push_back(pCandidateScan);
        }
      }
      // 当候选匹配帧与当前帧距离比较远了
      else
      {
        // return chain if it is long "enough"
        // 如果候选闭环帧集合数量够多了，认为找到了闭环
        if (chain.size() >= m_pMapper->m_pLoopMatchMinimumChainSize->GetValue())
        {
          return chain;
        }
        else
        {
          chain.clear();
        }
      }
    }

    return chain;
  }
  /**
   * 执行BA图优化，更新图中所有节点位姿
  */
  void MapperGraph::CorrectPoses()
  {
    // optimize scans!
    ScanSolver* pSolver = m_pMapper->m_pScanOptimizer;
    if (pSolver != NULL)
    {
      pSolver->Compute();

      const_forEach(ScanSolver::IdPoseVector, &pSolver->GetCorrections())
      {
        m_pMapper->m_pMapperSensorManager->GetScan(iter->first)->SetSensorPose(iter->second);
      }

      pSolver->Clear();
    }
  }

  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////

  /**
   * Default constructor
   */
  Mapper::Mapper()
    : Module("Mapper")
    , m_Initialized(false)
    , m_pSequentialScanMatcher(NULL)
    , m_pMapperSensorManager(NULL)
    , m_pGraph(NULL)
    , m_pScanOptimizer(NULL)
  {
    InitializeParameters();
  }

  /**
   * Default constructor
   */
  Mapper::Mapper(const std::string& rName)
    : Module(rName)
    , m_Initialized(false)
    , m_pSequentialScanMatcher(NULL)
    , m_pMapperSensorManager(NULL)
    , m_pGraph(NULL)
    , m_pScanOptimizer(NULL)
  {
    InitializeParameters();
  }

  /**
   * Destructor
   */
  Mapper::~Mapper()
  {
    Reset();

    delete m_pMapperSensorManager;
  }

  void Mapper::InitializeParameters()
  {
    m_pUseScanMatching = new Parameter<kt_bool>(
        "UseScanMatching",
        "When set to true, the mapper will use a scan matching algorithm. "
        "In most real-world situations this should be set to true so that the "
        "mapper algorithm can correct for noise and errors in odometry and "
        "scan data. In some simulator environments where the simulated scan "
        "and odometry data are very accurate, the scan matching algorithm can "
        "produce worse results. In those cases set this to false to improve "
        "results.",
        true,
        GetParameterManager());

    m_pUseScanBarycenter = new Parameter<kt_bool>(
        "UseScanBarycenter",
        "Use the barycenter of scan endpoints to define distances between "
        "scans.",
        true, GetParameterManager());

    m_pMinimumTimeInterval = new Parameter<kt_double>(
        "MinimumTimeInterval",
        "Sets the minimum time between scans. If a new scan's time stamp is "
        "longer than MinimumTimeInterval from the previously processed scan, "
        "the mapper will use the data from the new scan. Otherwise, it will "
        "discard the new scan if it also does not meet the minimum travel "
        "distance and heading requirements. For performance reasons, it is "
        "generally it is a good idea to only process scans if a reasonable "
        "amount of time has passed. This parameter is particularly useful "
        "when there is a need to process scans while the robot is stationary.",
        3600, GetParameterManager());

    m_pMinimumTravelDistance = new Parameter<kt_double>(
        "MinimumTravelDistance",
        "Sets the minimum travel between scans.  If a new scan's position is "
        "more than minimumTravelDistance from the previous scan, the mapper "
        "will use the data from the new scan. Otherwise, it will discard the "
        "new scan if it also does not meet the minimum change in heading "
        "requirement. For performance reasons, generally it is a good idea to "
        "only process scans if the robot has moved a reasonable amount.",
        0.2, GetParameterManager());

    m_pMinimumTravelHeading = new Parameter<kt_double>(
        "MinimumTravelHeading",
        "Sets the minimum heading change between scans. If a new scan's "
        "heading is more than MinimumTravelHeading from the previous scan, the "
        "mapper will use the data from the new scan.  Otherwise, it will "
        "discard the new scan if it also does not meet the minimum travel "
        "distance requirement. For performance reasons, generally it is a good "
        "idea to only process scans if the robot has moved a reasonable "
        "amount.",
        math::DegreesToRadians(10), GetParameterManager());

    m_pScanBufferSize = new Parameter<kt_int32u>(
        "ScanBufferSize",
        "Scan buffer size is the length of the scan chain stored for scan "
        "matching. \"ScanBufferSize\" should be set to approximately "
        "\"ScanBufferMaximumScanDistance\" / \"MinimumTravelDistance\". The "
        "idea is to get an area approximately 20 meters long for scan "
        "matching. For example, if we add scans every MinimumTravelDistance == "
        "0.3 meters, then \"scanBufferSize\" should be 20 / 0.3 = 67.)",
        70, GetParameterManager());

    m_pScanBufferMaximumScanDistance = new Parameter<kt_double>(
        "ScanBufferMaximumScanDistance",
        "Scan buffer maximum scan distance is the maximum distance between the "
        "first and last scans in the scan chain stored for matching.",
        20.0, GetParameterManager());

    m_pLinkMatchMinimumResponseFine = new Parameter<kt_double>(
        "LinkMatchMinimumResponseFine",
        "Scans are linked only if the correlation response value is greater "
        "than this value.",
        0.8, GetParameterManager());

    m_pLinkScanMaximumDistance = new Parameter<kt_double>(
        "LinkScanMaximumDistance",
        "Maximum distance between linked scans.  Scans that are farther apart "
        "will not be linked regardless of the correlation response value.",
        10.0, GetParameterManager());

    m_pLoopSearchMaximumDistance = new Parameter<kt_double>(
        "LoopSearchMaximumDistance",
        "Scans less than this distance from the current position will be "
        "considered for a match in loop closure.",
        4.0, GetParameterManager());

    m_pDoLoopClosing = new Parameter<kt_bool>(
        "DoLoopClosing",
        "Enable/disable loop closure.",
        true, GetParameterManager());

    m_pLoopMatchMinimumChainSize = new Parameter<kt_int32u>(
        "LoopMatchMinimumChainSize",
        "When the loop closure detection finds a candidate it must be part of "
        "a large set of linked scans. If the chain of scans is less than this "
        "value we do not attempt to close the loop.",
        10, GetParameterManager());

    m_pLoopMatchMaximumVarianceCoarse = new Parameter<kt_double>(
        "LoopMatchMaximumVarianceCoarse",
        "The co-variance values for a possible loop closure have to be less "
        "than this value to consider a viable solution. This applies to the "
        "coarse search.",
        math::Square(0.4), GetParameterManager());

    m_pLoopMatchMinimumResponseCoarse = new Parameter<kt_double>(
        "LoopMatchMinimumResponseCoarse",
        "If response is larger then this, then initiate loop closure search at "
        "the coarse resolution.",
        0.8, GetParameterManager());

    m_pLoopMatchMinimumResponseFine = new Parameter<kt_double>(
        "LoopMatchMinimumResponseFine",
        "If response is larger then this, then initiate loop closure search at "
        "the fine resolution.",
        0.8, GetParameterManager());

    //////////////////////////////////////////////////////////////////////////////
    //    CorrelationParameters correlationParameters;

    m_pCorrelationSearchSpaceDimension = new Parameter<kt_double>(
        "CorrelationSearchSpaceDimension",
        "The size of the search grid used by the matcher. The search grid will "
        "have the size CorrelationSearchSpaceDimension * "
        "CorrelationSearchSpaceDimension",
        0.3, GetParameterManager());

    m_pCorrelationSearchSpaceResolution = new Parameter<kt_double>(
        "CorrelationSearchSpaceResolution",
        "The resolution (size of a grid cell) of the correlation grid.",
        0.01, GetParameterManager());

    m_pCorrelationSearchSpaceSmearDeviation = new Parameter<kt_double>(
        "CorrelationSearchSpaceSmearDeviation",
        "The point readings are smeared by this value in X and Y to create a "
        "smoother response.",
        0.03, GetParameterManager());


    //////////////////////////////////////////////////////////////////////////////
    //    CorrelationParameters loopCorrelationParameters;

    m_pLoopSearchSpaceDimension = new Parameter<kt_double>(
        "LoopSearchSpaceDimension",
        "The size of the search grid used by the matcher.",
        8.0, GetParameterManager());

    m_pLoopSearchSpaceResolution = new Parameter<kt_double>(
        "LoopSearchSpaceResolution",
        "The resolution (size of a grid cell) of the correlation grid.",
        0.05, GetParameterManager());

    m_pLoopSearchSpaceSmearDeviation = new Parameter<kt_double>(
        "LoopSearchSpaceSmearDeviation",
        "The point readings are smeared by this value in X and Y to create a "
        "smoother response.",
        0.03, GetParameterManager());

    //////////////////////////////////////////////////////////////////////////////
    // ScanMatcherParameters;

    m_pDistanceVariancePenalty = new Parameter<kt_double>(
        "DistanceVariancePenalty",
        "Variance of penalty for deviating from odometry when scan-matching. "
        "The penalty is a multiplier (less than 1.0) is a function of the "
        "delta of the scan position being tested and the odometric pose.",
        math::Square(0.3), GetParameterManager());

    m_pAngleVariancePenalty = new Parameter<kt_double>(
        "AngleVariancePenalty",
        "See DistanceVariancePenalty.",
        math::Square(math::DegreesToRadians(20)), GetParameterManager());

    m_pFineSearchAngleOffset = new Parameter<kt_double>(
        "FineSearchAngleOffset",
        "The range of angles to search during a fine search.",
        math::DegreesToRadians(0.2), GetParameterManager());

    m_pCoarseSearchAngleOffset = new Parameter<kt_double>(
        "CoarseSearchAngleOffset",
        "The range of angles to search during a coarse search.",
        math::DegreesToRadians(20), GetParameterManager());

    m_pCoarseAngleResolution = new Parameter<kt_double>(
        "CoarseAngleResolution",
        "Resolution of angles to search during a coarse search.",
        math::DegreesToRadians(2), GetParameterManager());

    m_pMinimumAnglePenalty = new Parameter<kt_double>(
        "MinimumAnglePenalty",
        "Minimum value of the angle penalty multiplier so scores do not become "
        "too small.",
        0.9, GetParameterManager());

    m_pMinimumDistancePenalty = new Parameter<kt_double>(
        "MinimumDistancePenalty",
        "Minimum value of the distance penalty multiplier so scores do not "
        "become too small.",
        0.5, GetParameterManager());

    m_pUseResponseExpansion = new Parameter<kt_bool>(
        "UseResponseExpansion",
        "Whether to increase the search space if no good matches are initially "
        "found.",
        false, GetParameterManager());
  }
  /* Adding in getters and setters here for easy parameter access */

  // General Parameters

  bool Mapper::getParamUseScanMatching()
  {
    return static_cast<bool>(m_pUseScanMatching->GetValue());
  }

  bool Mapper::getParamUseScanBarycenter()
  {
    return static_cast<bool>(m_pUseScanBarycenter->GetValue());
  }

  double Mapper::getParamMinimumTimeInterval()
  {
    return static_cast<double>(m_pMinimumTimeInterval->GetValue());
  }

  double Mapper::getParamMinimumTravelDistance()
  {
    return static_cast<double>(m_pMinimumTravelDistance->GetValue());
  }

  double Mapper::getParamMinimumTravelHeading()
  {
    return math::RadiansToDegrees(static_cast<double>(m_pMinimumTravelHeading->GetValue()));
  }

  int Mapper::getParamScanBufferSize()
  {
    return static_cast<int>(m_pScanBufferSize->GetValue());
  }

  double Mapper::getParamScanBufferMaximumScanDistance()
  {
    return static_cast<double>(m_pScanBufferMaximumScanDistance->GetValue());
  }

  double Mapper::getParamLinkMatchMinimumResponseFine()
  {
    return static_cast<double>(m_pLinkMatchMinimumResponseFine->GetValue());
  }

  double Mapper::getParamLinkScanMaximumDistance()
  {
    return static_cast<double>(m_pLinkScanMaximumDistance->GetValue());
  }

  double Mapper::getParamLoopSearchMaximumDistance()
  {
    return static_cast<double>(m_pLoopSearchMaximumDistance->GetValue());
  }

  bool Mapper::getParamDoLoopClosing()
  {
    return static_cast<bool>(m_pDoLoopClosing->GetValue());
  }

  int Mapper::getParamLoopMatchMinimumChainSize()
  {
    return static_cast<int>(m_pLoopMatchMinimumChainSize->GetValue());
  }

  double Mapper::getParamLoopMatchMaximumVarianceCoarse()
  {
    return static_cast<double>(std::sqrt(m_pLoopMatchMaximumVarianceCoarse->GetValue()));
  }

  double Mapper::getParamLoopMatchMinimumResponseCoarse()
  {
    return static_cast<double>(m_pLoopMatchMinimumResponseCoarse->GetValue());
  }

  double Mapper::getParamLoopMatchMinimumResponseFine()
  {
    return static_cast<double>(m_pLoopMatchMinimumResponseFine->GetValue());
  }

  // Correlation Parameters - Correlation Parameters

  double Mapper::getParamCorrelationSearchSpaceDimension()
  {
    return static_cast<double>(m_pCorrelationSearchSpaceDimension->GetValue());
  }

  double Mapper::getParamCorrelationSearchSpaceResolution()
  {
    return static_cast<double>(m_pCorrelationSearchSpaceResolution->GetValue());
  }

  double Mapper::getParamCorrelationSearchSpaceSmearDeviation()
  {
    return static_cast<double>(m_pCorrelationSearchSpaceSmearDeviation->GetValue());
  }

  // Correlation Parameters - Loop Correlation Parameters

  double Mapper::getParamLoopSearchSpaceDimension()
  {
    return static_cast<double>(m_pLoopSearchSpaceDimension->GetValue());
  }

  double Mapper::getParamLoopSearchSpaceResolution()
  {
    return static_cast<double>(m_pLoopSearchSpaceResolution->GetValue());
  }

  double Mapper::getParamLoopSearchSpaceSmearDeviation()
  {
    return static_cast<double>(m_pLoopSearchSpaceSmearDeviation->GetValue());
  }

  // ScanMatcher Parameters

  double Mapper::getParamDistanceVariancePenalty()
  {
    return std::sqrt(static_cast<double>(m_pDistanceVariancePenalty->GetValue()));
  }

  double Mapper::getParamAngleVariancePenalty()
  {
    return std::sqrt(static_cast<double>(m_pAngleVariancePenalty->GetValue()));
  }

  double Mapper::getParamFineSearchAngleOffset()
  {
    return static_cast<double>(m_pFineSearchAngleOffset->GetValue());
  }

  double Mapper::getParamCoarseSearchAngleOffset()
  {
    return static_cast<double>(m_pCoarseSearchAngleOffset->GetValue());
  }

  double Mapper::getParamCoarseAngleResolution()
  {
    return static_cast<double>(m_pCoarseAngleResolution->GetValue());
  }

  double Mapper::getParamMinimumAnglePenalty()
  {
    return static_cast<double>(m_pMinimumAnglePenalty->GetValue());
  }

  double Mapper::getParamMinimumDistancePenalty()
  {
    return static_cast<double>(m_pMinimumDistancePenalty->GetValue());
  }

  bool Mapper::getParamUseResponseExpansion()
  {
    return static_cast<bool>(m_pUseResponseExpansion->GetValue());
  }

  /* Setters for parameters */
  // General Parameters
  void Mapper::setParamUseScanMatching(bool b)
  {
    m_pUseScanMatching->SetValue((kt_bool)b);
  }

  void Mapper::setParamUseScanBarycenter(bool b)
  {
    m_pUseScanBarycenter->SetValue((kt_bool)b);
  }

  void Mapper::setParamMinimumTimeInterval(double d)
  {
    m_pMinimumTimeInterval->SetValue((kt_double)d);
  }

  void Mapper::setParamMinimumTravelDistance(double d)
  {
    m_pMinimumTravelDistance->SetValue((kt_double)d);
  }

  void Mapper::setParamMinimumTravelHeading(double d)
  {
    m_pMinimumTravelHeading->SetValue((kt_double)d);
  }

  void Mapper::setParamScanBufferSize(int i)
  {
    m_pScanBufferSize->SetValue((kt_int32u)i);
  }

  void Mapper::setParamScanBufferMaximumScanDistance(double d)
  {
    m_pScanBufferMaximumScanDistance->SetValue((kt_double)d);
  }

  void Mapper::setParamLinkMatchMinimumResponseFine(double d)
  {
    m_pLinkMatchMinimumResponseFine->SetValue((kt_double)d);
  }

  void Mapper::setParamLinkScanMaximumDistance(double d)
  {
    m_pLinkScanMaximumDistance->SetValue((kt_double)d);
  }

  void Mapper::setParamLoopSearchMaximumDistance(double d)
  {
    m_pLoopSearchMaximumDistance->SetValue((kt_double)d);
  }

  void Mapper::setParamDoLoopClosing(bool b)
  {
    m_pDoLoopClosing->SetValue((kt_bool)b);
  }

  void Mapper::setParamLoopMatchMinimumChainSize(int i)
  {
    m_pLoopMatchMinimumChainSize->SetValue((kt_int32u)i);
  }

  void Mapper::setParamLoopMatchMaximumVarianceCoarse(double d)
  {
    m_pLoopMatchMaximumVarianceCoarse->SetValue((kt_double)math::Square(d));
  }

  void Mapper::setParamLoopMatchMinimumResponseCoarse(double d)
  {
    m_pLoopMatchMinimumResponseCoarse->SetValue((kt_double)d);
  }

  void Mapper::setParamLoopMatchMinimumResponseFine(double d)
  {
    m_pLoopMatchMinimumResponseFine->SetValue((kt_double)d);
  }

  // Correlation Parameters - Correlation Parameters
  void Mapper::setParamCorrelationSearchSpaceDimension(double d)
  {
    m_pCorrelationSearchSpaceDimension->SetValue((kt_double)d);
  }

  void Mapper::setParamCorrelationSearchSpaceResolution(double d)
  {
    m_pCorrelationSearchSpaceResolution->SetValue((kt_double)d);
  }

  void Mapper::setParamCorrelationSearchSpaceSmearDeviation(double d)
  {
    m_pCorrelationSearchSpaceSmearDeviation->SetValue((kt_double)d);
  }


  // Correlation Parameters - Loop Closure Parameters
  void Mapper::setParamLoopSearchSpaceDimension(double d)
  {
    m_pLoopSearchSpaceDimension->SetValue((kt_double)d);
  }

  void Mapper::setParamLoopSearchSpaceResolution(double d)
  {
    m_pLoopSearchSpaceResolution->SetValue((kt_double)d);
  }

  void Mapper::setParamLoopSearchSpaceSmearDeviation(double d)
  {
    m_pLoopSearchSpaceSmearDeviation->SetValue((kt_double)d);
  }


  // Scan Matcher Parameters
  void Mapper::setParamDistanceVariancePenalty(double d)
  {
    m_pDistanceVariancePenalty->SetValue((kt_double)math::Square(d));
  }

  void Mapper::setParamAngleVariancePenalty(double d)
  {
    m_pAngleVariancePenalty->SetValue((kt_double)math::Square(d));
  }

  void Mapper::setParamFineSearchAngleOffset(double d)
  {
    m_pFineSearchAngleOffset->SetValue((kt_double)d);
  }

  void Mapper::setParamCoarseSearchAngleOffset(double d)
  {
    m_pCoarseSearchAngleOffset->SetValue((kt_double)d);
  }

  void Mapper::setParamCoarseAngleResolution(double d)
  {
    m_pCoarseAngleResolution->SetValue((kt_double)d);
  }

  void Mapper::setParamMinimumAnglePenalty(double d)
  {
    m_pMinimumAnglePenalty->SetValue((kt_double)d);
  }

  void Mapper::setParamMinimumDistancePenalty(double d)
  {
    m_pMinimumDistancePenalty->SetValue((kt_double)d);
  }

  void Mapper::setParamUseResponseExpansion(bool b)
  {
    m_pUseResponseExpansion->SetValue((kt_bool)b);
  }



  void Mapper::Initialize(kt_double rangeThreshold)
  {
    if (m_Initialized == false)
    {
      // create sequential scan and loop matcher
      m_pSequentialScanMatcher = ScanMatcher::Create(this,
                                                    m_pCorrelationSearchSpaceDimension->GetValue(),
                                                    m_pCorrelationSearchSpaceResolution->GetValue(),
                                                    m_pCorrelationSearchSpaceSmearDeviation->GetValue(),
                                                    rangeThreshold);
      assert(m_pSequentialScanMatcher);

      m_pMapperSensorManager = new MapperSensorManager(m_pScanBufferSize->GetValue(),
                                                       m_pScanBufferMaximumScanDistance->GetValue());

      m_pGraph = new MapperGraph(this, rangeThreshold);

      m_Initialized = true;
    }
  }

  void Mapper::Reset()
  {
    delete m_pSequentialScanMatcher;
    m_pSequentialScanMatcher = NULL;

    delete m_pGraph;
    m_pGraph = NULL;

    delete m_pMapperSensorManager;
    m_pMapperSensorManager = NULL;

    m_Initialized = false;
  }

  kt_bool Mapper::Process(Object*  /*pObject*/)
  {
    return true;
  }

  /**
   * 主流程
   * 1、用前一帧优化前后的位姿变换，来校正当前帧优化前的位姿，得到一个初始位姿
   * 2、判断两帧间隔是否足够大，不够大则直接丢弃这一帧
   *    满足一下任意一个条件，返回true，都不满足则返回false
   *    a.第一帧
   *    b.两帧之间时间间隔足够长
   *    c.两帧yaw角差量足够大
   *    d.两帧位置差量足够大
   * 3、基于相关方法的scan-to-map匹配
   *    1) 用前一帧优化前后的位姿变换，来校正当前帧优化前的位姿，得到一个初始位姿
   *    2) 提取局部map中当前帧也可以看到的点，对应栅格设置为占用
   *       剔除相对于当前帧属于物体背面的点，也就是近邻帧与当前帧在物体两面
   *    3) scan-to-map匹配
   *       a. 创建旋转角度-激光点旋转后相对于当前帧位置的偏移量
   *          创建索引表，计算当前帧激光点，在各个旋转角度变换下，新的激光点位置相对于当前帧位置（机器人的位置）的坐标偏移量
   *          不管当前帧在什么位置，只要旋转角度一定，那么坐标偏移量也是确定的。
   *          目的是后面在搜索空间中对当前帧位姿施加不同的x平移、y平移之后再施加一个旋转的操作，不用重新计算每个搜索位置处旋转某个角度之后对应的点坐标，只需要加上对应旋转角度下算好的坐标偏移量即可。
   *       b. 遍历pose搜索空间计算响应值     
   *          i. 当前帧激光点集合在某个pose变换下得到新的位置，统计在局部map栅格中这些新位置处占用的数量，占用越多响应值越高，表示匹配越好
   *          ii. 惩罚一下搜索偏移量，越远响应值折扣越多（当前帧搜索距离、角度偏移）
   *       c. 找到响应值最大的pose，如果有多个最佳响应pose，那么对pose求个均值，得到最终pose
   *       d. 计算位姿的协方差
   *    4) 如果响应值是0，没匹配上，扩大角度搜索范围，增加20°、40°、60°，再执行一次scan-to-map匹配
   *    5) 在scan-to-map优化之后的位姿基础上，再缩小搜索空间优化一次，搜索范围减半，角度区间减半，执行scan-to-map匹配
   * 4、构建BA图的顶点、边
   *    1) 与前一帧建立一条边
   *    2) 在局部map帧集合中选取与当前帧最接近的一帧，建立边
   *    3) 在相邻帧的chain中取与当前帧最接近的一帧，建立边
   *       a. 在顶点图中广搜找当前帧的相邻帧，然后对每个相邻帧截取前后时间范围内、且距离小于阈值的一段帧集合（称为chain）
   *       b. 用当前帧与chain进行scan-to-map匹配
   *       c. 如果匹配较好，在chain里面选最近一帧，建立边
   *       d. 对所有数量足够的chain执行前面2、3步
   *       注：相当于当前帧跟比较近，但又没有直接关系的帧（时间上相邻帧，局部匹配帧，闭环帧）建立边
   *     如果闭环成功，也会添加一条边。（注：不在这个方法里面）
   *     4) 最后设置当前帧pose为第3步中计算得到的优化pose集合（包含当前帧pose）的均值
   * 5、添加一帧数据并更新局部map，维护当前帧时空维度上临近的帧集合
   * 6、闭环检测
   *     1) 遍历所有历史帧，提取候选闭环帧集合
   *         a. 在已构造的顶点图中，广搜找与当前顶点帧临近的所有顶点帧，记为临近帧集合
   *         b. 按时序遍历所有历史帧，距离够近的帧加入候选闭环集合，一旦发现一帧在临近帧集合中，清空候选闭环集合，接着往下找
   *         c. 直到距离超过阈值，如果候选闭环集合数量够多，认为合格
   *     2) 当前帧与候选闭环帧集合，执行scan-to-map匹配
   *     3) 如果响应够大，方差够小，执行二次scan-to-map匹配
   *     4) 如果确认闭环，在候选闭环帧集合中选取与当前帧最接近的一帧，建立边，同时更新当前帧位姿
   *     5) 循环提取下一个闭环候选帧集合，寻找闭环
   * 7、找到闭环，则执行BA图优化，更新所有顶点位姿
  */
  kt_bool Mapper::Process(LocalizedRangeScan* pScan)
  {
    if (pScan != NULL)
    {
      karto::LaserRangeFinder* pLaserRangeFinder = pScan->GetLaserRangeFinder();

      // validate scan
      if (pLaserRangeFinder == NULL || pScan == NULL || pLaserRangeFinder->Validate(pScan) == false)
      {
        return false;
      }

      if (m_Initialized == false)
      {
        // initialize mapper with range threshold from device
        Initialize(pLaserRangeFinder->GetRangeThreshold());
      }

      // 前一帧激光数据
      LocalizedRangeScan* pLastScan = m_pMapperSensorManager->GetLastScan(pScan->GetSensorName());

      // 用前一帧优化前后的位姿变换，来校正当前帧优化前的位姿，得到一个初始位姿
      if (pLastScan != NULL)
      {
        Transform lastTransform(pLastScan->GetOdometricPose(), pLastScan->GetCorrectedPose());
        pScan->SetCorrectedPose(lastTransform.TransformPose(pScan->GetOdometricPose()));
      }

      // 判断两帧间隔是否足够大，不够大则直接丢弃这一帧
      // 满足一下任意一个条件，返回true，都不满足则返回false
      // 1、第一帧
      // 2、两帧之间时间间隔足够长
      // 3、两帧yaw角差量足够大
      // 4、两帧位置差量足够大
      if (!HasMovedEnough(pScan, pLastScan))
      {
        return false;
      }

      Matrix3 covariance;
      covariance.SetToIdentity();

      // correct scan (if not first scan)
      if (m_pUseScanMatching->GetValue() && pLastScan != NULL)
      {
        Pose2 bestPose;
        // 基于相关方法的scan-to-map匹配
        // 1、用前一帧优化前后的位姿变换，来校正当前帧优化前的位姿，得到一个初始位姿
        // 2、提取局部map中当前帧也可以看到的点，对应栅格设置为占用
        //    剔除相对于当前帧属于物体背面的点，也就是近邻帧与当前帧在物体两面
        // 3、scan-to-map匹配
        //    1) 创建旋转角度-激光点旋转后相对于当前帧位置的偏移量
        //       创建索引表，计算当前帧激光点，在各个旋转角度变换下，新的激光点位置相对于当前帧位置（机器人的位置）的坐标偏移量
        //       不管当前帧在什么位置，只要旋转角度一定，那么坐标偏移量也是确定的。
        //       目的是后面在搜索空间中对当前帧位姿施加不同的x平移、y平移之后再施加一个旋转的操作，不用重新计算每个搜索位置处旋转某个角度之后对应的点坐标，只需要加上对应旋转角度下算好的坐标偏移量即可。
        //    2) 遍历pose搜索空间计算响应值     
        //       a. 当前帧激光点集合在某个pose变换下得到新的位置，统计在局部map栅格中这些新位置处占用的数量，占用越多响应值越高，表示匹配越好
        //       b. 惩罚一下搜索偏移量，越远响应值折扣越多（当前帧搜索距离、角度偏移）
        //    3) 找到响应值最大的pose，如果有多个最佳响应pose，那么对pose求个均值，得到最终pose
        //    4) 计算位姿的协方差
        // 4、如果响应值是0，没匹配上，扩大角度搜索范围，增加20°、40°、60°，再执行一次scan-to-map匹配
        // 5、在scan-to-map优化之后的位姿基础上，再缩小搜索空间优化一次，搜索范围减半，角度区间减半，执行scan-to-map匹配
        m_pSequentialScanMatcher->MatchScan(pScan,
                                            // 局部map，维护当前帧时空维度上临近的帧集合
                                            m_pMapperSensorManager->GetRunningScans(pScan->GetSensorName()),
                                                                                    bestPose,
                                                                                    covariance);
        // 更新当前帧位姿
        pScan->SetSensorPose(bestPose);
      }

      // 添加当前帧
      m_pMapperSensorManager->AddScan(pScan);

      if (m_pUseScanMatching->GetValue())
      {
        // add to graph
        // 加入BA图中，添加顶点、边
        m_pGraph->AddVertex(pScan);
        // 构建BA图的边
        // 1、与前一帧建立一条边
        // 2、在局部map帧集合中选取与当前帧最接近的一帧，建立边
        // 3、在相邻帧的chain中取与当前帧最接近的一帧，建立边
        //   1) 在顶点图中广搜找当前帧的相邻帧，然后对每个相邻帧截取前后时间范围内、且距离小于阈值的一段帧集合（称为chain）
        //   2) 用当前帧与chain进行scan-to-map匹配
        //   3) 如果匹配较好，在chain里面选最近一帧，建立边
        //   4) 对所有数量足够的chain执行前面2、3步
        //   注：相当于当前帧跟比较近，但又没有直接关系的帧（时间上相邻帧，局部匹配帧，闭环帧）建立边
        // 如果闭环成功，也会添加一条边。（注：不在这个方法里面）
        // 4、最后设置当前帧pose为第3步中计算得到的优化pose集合（包含当前帧pose）的均值
        m_pGraph->AddEdges(pScan, covariance);

        // 添加一帧数据并更新局部map，维护当前帧时空维度上临近的帧集合
        m_pMapperSensorManager->AddRunningScan(pScan);

        // 闭环检测，如果找到闭环，执行BA优化图顶点位姿
        if (m_pDoLoopClosing->GetValue())
        {
          std::vector<Name> deviceNames = m_pMapperSensorManager->GetSensorNames();
          const_forEach(std::vector<Name>, &deviceNames)
          {
            // 闭环检测
            // 1、遍历所有历史帧，提取候选闭环帧集合
            //     1) 在已构造的顶点图中，广搜找与当前顶点帧临近的所有顶点帧，记为临近帧集合
            //     2) 按时序遍历所有历史帧，距离够近的帧加入候选闭环集合，一旦发现一帧在临近帧集合中，清空候选闭环集合，接着往下找
            //     3) 直到距离超过阈值，如果候选闭环集合数量够多，认为合格
            // 2、当前帧与候选闭环帧集合，执行scan-to-map匹配
            // 3、如果响应够大，方差够小，执行二次scan-to-map匹配
            // 4、如果确认闭环，在候选闭环帧集合中选取与当前帧最接近的一帧，建立边，同时更新当前帧位姿
            // 5、循环提取下一个闭环候选帧集合，寻找闭环
            m_pGraph->TryCloseLoop(pScan, *iter);
          }
        }
      }

      m_pMapperSensorManager->SetLastScan(pScan);

      return true;
    }

    return false;
  }

  /**
   * 判断两帧间隔是否足够大，可以理解为关键帧
   * 满足一下任意一个条件，返回true，都不满足则返回false
   * 1、第一帧
   * 2、两帧之间时间间隔足够长
   * 3、两帧yaw角差量足够大
   * 4、两帧位置差量足够大
   */
  kt_bool Mapper::HasMovedEnough(LocalizedRangeScan* pScan, LocalizedRangeScan* pLastScan) const
  {
    // 第一帧
    if (pLastScan == NULL)
    {
      return true;
    }

    // 两帧之间时间间隔足够长
    kt_double timeInterval = pScan->GetTime() - pLastScan->GetTime();
    if (timeInterval >= m_pMinimumTimeInterval->GetValue())
    {
      return true;
    }

    // 两帧对应的优化前位姿
    Pose2 lastScannerPose = pLastScan->GetSensorAt(pLastScan->GetOdometricPose());
    Pose2 scannerPose = pScan->GetSensorAt(pScan->GetOdometricPose());

    // 两帧yaw角差量，是否足够大
    kt_double deltaHeading = math::NormalizeAngle(scannerPose.GetHeading() - lastScannerPose.GetHeading());
    if (fabs(deltaHeading) >= m_pMinimumTravelHeading->GetValue())
    {
      return true;
    }

    // 两帧位置差量，是否足够大
    kt_double squaredTravelDistance = lastScannerPose.GetPosition().SquaredDistance(scannerPose.GetPosition());
    if (squaredTravelDistance >= math::Square(m_pMinimumTravelDistance->GetValue()) - KT_TOLERANCE)
    {
      return true;
    }

    return false;
  }

  /**
   * Gets all the processed scans
   * @return all scans
   */
  const LocalizedRangeScanVector Mapper::GetAllProcessedScans() const
  {
    LocalizedRangeScanVector allScans;

    if (m_pMapperSensorManager != NULL)
    {
      allScans = m_pMapperSensorManager->GetAllScans();
    }

    return allScans;
  }

  /**
   * Adds a listener
   * @param pListener
   */
  void Mapper::AddListener(MapperListener* pListener)
  {
    m_Listeners.push_back(pListener);
  }

  /**
   * Removes a listener
   * @param pListener
   */
  void Mapper::RemoveListener(MapperListener* pListener)
  {
    std::vector<MapperListener*>::iterator iter = std::find(m_Listeners.begin(), m_Listeners.end(), pListener);
    if (iter != m_Listeners.end())
    {
      m_Listeners.erase(iter);
    }
  }

  void Mapper::FireInfo(const std::string& rInfo) const
  {
    const_forEach(std::vector<MapperListener*>, &m_Listeners)
    {
      (*iter)->Info(rInfo);
    }
  }

  void Mapper::FireDebug(const std::string& rInfo) const
  {
    const_forEach(std::vector<MapperListener*>, &m_Listeners)
    {
      MapperDebugListener* pListener = dynamic_cast<MapperDebugListener*>(*iter);

      if (pListener != NULL)
      {
        pListener->Debug(rInfo);
      }
    }
  }

  void Mapper::FireLoopClosureCheck(const std::string& rInfo) const
  {
    const_forEach(std::vector<MapperListener*>, &m_Listeners)
    {
      MapperLoopClosureListener* pListener = dynamic_cast<MapperLoopClosureListener*>(*iter);

      if (pListener != NULL)
      {
        pListener->LoopClosureCheck(rInfo);
      }
    }
  }

  void Mapper::FireBeginLoopClosure(const std::string& rInfo) const
  {
    const_forEach(std::vector<MapperListener*>, &m_Listeners)
    {
      MapperLoopClosureListener* pListener = dynamic_cast<MapperLoopClosureListener*>(*iter);

      if (pListener != NULL)
      {
        pListener->BeginLoopClosure(rInfo);
      }
    }
  }

  void Mapper::FireEndLoopClosure(const std::string& rInfo) const
  {
    const_forEach(std::vector<MapperListener*>, &m_Listeners)
    {
      MapperLoopClosureListener* pListener = dynamic_cast<MapperLoopClosureListener*>(*iter);

      if (pListener != NULL)
      {
        pListener->EndLoopClosure(rInfo);
      }
    }
  }

  void Mapper::SetScanSolver(ScanSolver* pScanOptimizer)
  {
    m_pScanOptimizer = pScanOptimizer;
  }

  MapperGraph* Mapper::GetGraph() const
  {
    return m_pGraph;
  }

  ScanMatcher* Mapper::GetSequentialScanMatcher() const
  {
    return m_pSequentialScanMatcher;
  }

  ScanMatcher* Mapper::GetLoopScanMatcher() const
  {
    return m_pGraph->GetLoopScanMatcher();
  }
}  // namespace karto
