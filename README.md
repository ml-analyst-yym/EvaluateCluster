# 作者第一次使用GitHub放包，可能有些问题，在该文件中一一解决

# Notes：

该函数用于确定Seurat聚类时的最佳聚类数（Seurat V5）

函数包含4种算法：rogue、lisi、mrtree、CH_index（均可自行上网搜索用途）。

推荐在FindCluster时resolution=seq(0.1,1,0.1),可以设置更多分辨率，但运行速度会减慢。

其中rogue、mrtree运行较慢,均可在函数输入时选择是否运行（默认全部运行）。

rogue: 20个分辨率单核约2.5小时。##该函数主流用于评估cluster纯净度。

PS：Windows下我想用多核运行，奈何这个函数多核跑总是报错所以定死了单核。Linux为多核运行。

mrtree：20个分辨率8核约20分钟。##mrtree是CHOIR的依赖包，主要使用多分辨率下的兰德系数（AMRI）判断聚类合理性。

PS：CHOIR我也想用的，但可能是我的问题，教程学了整整一天才学明白，成功运行后，会在检查结果的步骤卡死，导致完全没有结果输出。

参数data为SeuratObject，cores在Windows下默认为1，不可更改，否则报错。服务器中可自行设定。

要求：输入的Seurat对象已经运行完RunUMAP。

最后：该包是多种检验Seurat聚类效果R包的整合包，仅CH_index为作者手搓。操作简单、源码易懂，出现BUG概不负责。

# 苦口婆心：

批次效应可能会影响聚类得分，因此在运行该R包前，应自行检查数据是否存在批次效应并去除。

一般单细胞只会在亚聚类的步骤用到该包，因此聚类数应该是合理的，所以尽量不要横跨0.1-2的分辨率，没有什么道理。而且理论上聚类数越多，总得分应该更高。所以0.2的得分和2.0的得分应该不能一同衡量的。我只推荐约5-10个分辨率的跨度下衡量得分的高低。因为我认为综合了四种算法的综合得分理应能反应相似分辨率下的微小聚类差异。

如果有用户希望一步完成图谱构建，可以直接将分辨率拉到2-3之间。同时，Seurat的FindCluster是支持0.01步长的，可以检测高分辨率下亚群的更微小的差异。

如果用户不信任我提供的综合得分，我也将四种算法计算的每种分辨率下的得分全部给出，可以自行挑选。理论上，聚类效果越好，lisi的得分越低，其他三种算法的得分越高。因CH系数在样本数较高时分数和其他三种算法得分不在同一水平，已对其进行标准化。同时，该R包涉及的每一个算法都是非常优秀的算法，可以单独使用！

四个算法中的三个都是集成了别人R包中的算法，可以单独使用，不用在意综合得分！


#感谢您的使用！

# 下载R包

devtools::install_github("ml-analyst-yym/EvaluateCluster")

# 该R包有三个依赖包

devtools::install_github("pengminshi/mrtree")#mrtree包（多分辨率下的兰德系数）

devtools::install_github("PaulingLiu/ROGUE")#ROGUE包（cluster纯净度）

devtools::install_github("immunogenomics/lisi")#lisi包（局部逆辛普森指数）

# 快速使用

EvaluateClusters(data)#只会计算CH系数

EvaluateClusters(data,rogue=T,mrtree=T,lisi=T)#如果安装了前面三个依赖包，会综合四种方法，最终返回一个列表，包含四种方法的得分和综合得分。

EvaluateClusters(data,plot=F)#不返回图像（图像可以查看得分差异，帮助选择较为合适的分辨率）
