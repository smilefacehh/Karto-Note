# Karto-Note
Karto Slam源码详细注释。

Karto是一种2D激光SLAM解决方案，它是基于稀疏图优化的方法，带闭环检测。

Karto的代码个人觉得不是很好阅读，本文对Karto的主流程代码进行分析，主要对局部map的维护、基于相关方法的scan-to-map匹配、闭环检测、BA图构建、计数法更新栅格地图等几个方面进行讲解。不涉及Jacobian的构造，以及图优化细节。希望可以帮到有需要的同学，错误的地方请您批评指正。

## 目录（知乎）
- [Karto源码解析(一)：工程运行](https://zhuanlan.zhihu.com/p/350852337)
- [Karto源码解析(二)：代码解析](https://zhuanlan.zhihu.com/p/352388229)

## 一、主流程

先放一个自己画的Karto流程简图。

![Image](https://github.com/smilefacehh/Karto-Note/blob/main/karto_slam.png)

每当新来一帧数据，都会执行一次完整的主流程。

1. 首先判断当前帧与前一帧是否间隔足够大，不够大则直接丢弃，可以理解为只处理关键帧；

2. 接下来要对当前帧位姿进行校正，采用scan-to-map的配准方法。map是选择当前帧时空维度上相邻的帧组成的集合，我们称之为局部map。在当前帧位姿搜索空间中，对当前帧激光点进行变换，与局部map的激光点集合进行匹配，计算响应值，选择响应值最高的匹配位姿，作为当前帧位姿；

3. 前面两步相当于SLAM的前端，得到了当前帧的优化位姿。接着，构建BA优化图，添加当前帧为新的顶点，然后添加一些边。一共要添加四种边：

 - 与前一帧建立边；
 - 与局部map中最近的一帧建立边；
 - 与当前帧的广搜邻近帧的chain，建立边，后边再解释是啥；
 - 与闭环候选帧中最接近的一帧建立边；

4. 闭环检测。在历史帧中查找当前帧时间上相隔较远，空间上相隔较近的帧，作为候选闭环帧，与当前帧进行scan-to-map匹配。如果匹配较好，认为找到闭环了；

5. 如果找到闭环，就执行BA图优化，更新所有顶点位姿。

## 二、局部map

局部map是指与当前帧在时间、空间维度上邻近的帧集合。用于与当前帧进行scan-to-map匹配，优化当前帧位姿。局部map中的激光点需要过滤一下，不能直接全部拿来匹配。因为即便局部帧与当前帧很近，但仍有可能局部帧与当前帧处在一个物体的两面，比如穿过一堵墙，那么局部帧看到墙这一面的点，当前帧已经到了墙另一面，这些点是看不到的。这样的点需要剔除，如何辨识这些点，可以看看FindValidPoints方法。

相关函数包括：AddRunningScan、AddScans、AddScan、FindValidPoints。

## 三、基于相关方法的scan-to-map匹配

scan-to-map匹配是激光雷达中一个基础的概念，scan是当前帧，map是当前帧时空维度上临近的帧组成的集合。map已经经过不断地优化，有一定精度了，而当前帧只有一个初始的估计位姿。至于怎么来的，可以是用前一帧位姿+运动计算，也可以是有外部传感器IMU、轮式里程计等得到的，总之它是不准确的。那么我们可以将当前帧的激光点，随意变换一下，再与map的激光点进行匹配。如果匹配的特别好，我们认为这个随意变换的pose就是我们需要的校正位姿。当然，随意变换是有个限制的，不能变换的特别特别大吧，所以定义了一个变化范围，称之为搜索空间。包括平移（x、y）、旋转（theta）维度上的搜索。

伪代码大致是这样的：

```
scan, map
bestPose = (x,y,theta)
bestScore = Score(scan, map)
for dx in (-dx, dx):
    for dy in (-dy, dy):
        for dtheta in (-dtheta, dtheta):
            T = (x + dx, y + dy, theta + dtheta)
            scan = T * scan
            score = Score(scan, map)
            if score > bestScore:
                bestScore = score
                bestPose = T
```

如何判断scan与map匹配的好不好呢？很简单，scan的激光点变换之后的位置，有map点接盘，不管是什么点，只要有个点就算匹配上了。统计一下匹配上的比例，用它来评判匹配好不好。为啥有效呢？可以想想。

另外要提一点，搜索空间中的角度问题。如果每平移一个（dx, dy）之后，再计算旋转dtheta角度对应的激光点位置，那计算DX*DY*DTheta次。能不能精简一下呢，可以的。我们知道旋转这个东西跟在哪没有关系，所以一开始就先计算好旋转一个dtheta角之后，激光点的偏移量，注意是偏移量。后面对于任何平移（dx, dy）之后的位置，加上这个偏移量就可以了。整个计算量就是DX*DY + DTheta次，少了很多。

scan-to-map匹配在很多地方都会用到。当前帧与局部map匹配，当前帧与闭环候选帧集合匹配，当前帧与广搜临近帧的chain匹配，都会调用这个方法。

相关函数包括：MatchScan、CorrelateScan、ComputeOffsets、GetResponse、ComputePositionalCovariance、ComputeAngularCovariance。

## 四、闭环检测

闭环就是机器人走了一圈（也有可能是好多圈），我们要在历史帧中找当前帧位置的帧，构建闭环。闭环帧跟当前帧的位置要很接近，时间呢又要求隔得比较远，不然就是相邻帧了。找到这些候选帧，称之为闭环map，再与当前帧进行scan-to-map匹配，优化当前帧的位姿。

相关函数包括：TryCloseLoop、FindPossibleLoopClosure。

## 五、BA图构建

这里属于SLAM后端部分了，BA图的顶点为历史帧，边则有很多情况。每当新来一帧，需要从图中挑出来以下这些顶点，与之建立一条边。

 - 与前一帧建立边；
 - 与局部map中最接近的一帧建立边；
 - 与闭环map中最接近的一帧建立边，可以理解为闭环边；
 - 这个复杂一点，首先在BA图中当前帧顶点位置处，广搜得到一些距离较近的帧，然后对于每一帧截取前后时间段范围内的帧，作为该邻近帧的chain。再在chain里面找一帧最近的帧，与当前帧建立边。目的都是尽可能多的把距离相近的帧关联起来，不然图过于稀疏。

相关函数包括：AddEdges、LinkNearChains、LinkChainToScan、FindNearChains、FindNearLinkedScans。

## 六、计数法更新栅格地图

每当新增一帧激光点数据，如何融入已有的地图中呢？我们为栅格定义了三种状态

```
typedef enum
{
  GridStates_Unknown = 0,
  GridStates_Occupied = 100,
  GridStates_Free = 255
} GridStates;
```

从激光发射器位置发射一束激光，打在物体表面，得到一个激光点数据。那么激光点位置记为占用，激光束穿过的的位置都是free，物体后面的区域是未知的。每当新增一帧激光点数据，对激光束路径上的点累加一次free或者occupied计数，用occupied/free来判断是否占用。

相关函数包括：CreateFromScans、AddScan、RayTrace、UpdateCell。

如有错误请您批评指正，希望内容对您有帮助，更多细节可以查看代码注释~

:) 如果对您有帮助，欢迎star~
