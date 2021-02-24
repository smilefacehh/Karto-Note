# Karto-Note
Karto Slam源码详细注释

Karto是一种2D激光SLAM解决方案，它是基于稀疏图优化的方法，带闭环检测。

Karto的代码个人觉得不是很好阅读，本文对Karto的主流程代码进行分析，主要对局部map的维护、基于相关方法的scan-to-map匹配、闭环检测、BA图构建、计数法更新栅格地图等几个方面进行讲解。不涉及Jacobian的构造，以及图优化细节。希望可以帮到有需要的同学，错误的地方请您批评指正。

## 目录（知乎）
- [Karto源码解析(一)：工程运行](https://zhuanlan.zhihu.com/p/350852337)
- [Karto源码解析(二)：代码解析](https://zhuanlan.zhihu.com/p/352388229)


如有错误请您批评指正，希望内容对您有帮助。完整内容可前往知乎文章查看~
