# 作者第一次使用GitHub放包，可能有些问题，在该文件中一一解决

# Notes：

该函数用于确定Seurat聚类时的最佳聚类数（Seurat V5）

推荐在FindCluster时resolution=seq(0.1,1,0.1),可以设置更多分辨率，但运行速度会减慢。

该包同名函数整合了现有的算法，使用 ROGUE 计算cluster纯度，确定分辨率下限， AUC 用于确定过度聚类阈值，确定分辨率上限。
clisi 与 mrtree 共同确定排名前三的分辨率。最后，使用簇的内部指标（包含 CH指数、 DB指数 、 轮廓系数）的秩排序选出最为推荐的分辨率。

参数data为SeuratObject，cores在Windows下默认为1，不可更改，否则报错。服务器中可自行设定。

如果在ROGUE结果中存在NA值，请提高rogue_span参数，具体解释参考loess的span参数。

要求：输入的Seurat对象已经运行完RunUMAP、JoinLayers，且包含多个分辨率（>=5）。

最后：该包是多种检验Seurat聚类效果R包的整合包。操作简单、源码易懂，出现BUG概不负责。

# 苦口婆心：

批次效应可能会影响聚类得分，因此在运行该R包前，应自行检查数据是否存在批次效应并去除。

一般单细胞只会在亚聚类的步骤用到该包，因此聚类数应该是合理的，所以尽量不要横跨0.1-2的分辨率，没有什么道理。我只推荐约5-10个分辨率（0.1）的跨度下衡量得分的高低。因为我认为相似分辨率下的微小聚类差异能通过不同算法反应。

如果有用户希望一步完成图谱构建，可以直接将分辨率拉到2-3之间。同时，Seurat的FindCluster是支持0.01步长的，可以检测高分辨率下亚群的更微小的差异。

如果用户不信任我提供的推荐分辨率，我也将各算法计算的每种分辨率下的得分全部给出，可以自行挑选。

# 对各算法的一些解释：

1、ROGUE：计算cluster纯度，越高越好。借用ROGUE包，参考文献中选择为>0.9（单个簇），但是对于整体的分辨率，文章中有一例重聚类后为0.861，因此我将其设定为>0.85

2、AUC：计算各簇满足条件的Marker数量，当出现第一次存在某簇为0，则说明该分辨率及以后为过度聚类。借用FindAllMarkers函数，参考文献中选择阈值为AUC>0.6，其他参数默认 #不推荐更改

3、lisi：计算cluster是否正确分类的指标，已转换成越高越好。借用lisi包，参考文献没有给出多类下的具体指标（说是最高得分为簇的数量，但是针对不同分辨率要取平均，得分会显著减少）

4、mrtree：计算校正后多分辨率的兰德系数（AMRI），是ARI的变体，越高越好。借用mrtree包，参考文献没有给出具体衡量指标

因为lisi与mrtree官方没有给出具体衡量指标，但对lisi取倒数后，二者得分范围均小于1。因此使用二者均值筛选Top5

5、簇内部指标：计算数学层面上的各簇内部指标，主要衡量簇内聚集度与簇间分离度。借用clusterCrit包，每种算法都有官方参考

因内部指标纳入算法较多，且取值范围、衡量方面、数值方向并不统一。因此使用秩排序作为内部指标的综合评价，从而给出推荐分辨率

# 该R包涉及的每一个算法都是非常优秀的算法，可以单独使用，不用在意结果的推荐！

# 感谢使用！

# 下载R包

devtools::install_github("ml-analyst-yym/EvaluateCluster")

# 该R包有三个依赖包

devtools::install_github("pengminshi/mrtree") # mrtree包（校正后多分辨率下的兰德系数）

devtools::install_github("PaulingLiu/ROGUE") # ROGUE包（cluster纯度）

devtools::install_github("immunogenomics/lisi") # lisi包（局部逆辛普森指数）

# 快速使用

res_eva <- EvaluateCluster(data, cores = 1, rogue_span=0.6) # cores请自行设定

res_auc <- CalculateAUC(data, cores = 1, min_pct = 0.05, logfc_threshold = 0.1, auc_cutoff = 0.6) # 不推荐更改auc_cutoff！

x <- Choose_res(result_eva, res_auc, rogue_threshold = 0.85) # rogue_threshold默认为0.85，参考文献中推荐为0.9

x$best_resolution
