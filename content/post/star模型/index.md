&emsp;

参见 [连享会推文-空间计量专题](https://www.lianxh.cn/blogs/29.html)

> **作者：** 袁子晴 (香港大学)     
> **邮箱：** <yzq0612@foxmail.com>

---

**目录**

[TOC]

---

## 1. 问题背景

**地理学第一定律**（Tobler's first law of geography）阐述了“所有事物都与其他事物相关，但是近处的事物比远处的事物更相关”。该定律是研究空间自相关性的基础。空间计量模型通过引入空间权重矩阵来研究单元之间的关联方式和关联程度，以此反映地理学第一定律。

但是大量的经验事实和实证研究表明**时间维度**也是一个很重要的因素。我们可以随着时间的推移在不同时点上收集空间横截面数据，从而得到两种类型的时空数据：**空间面板数据** (Spatial Panels) 和**混合空间截面数据** [Spatial Data (cross-sections) Pooled over Time]。

空间面板数据是在一段时间内对空间单元 (Spatial Units) 进行**追踪**调查所得到的，即使在不同的时间点上观测到的空间单元都是一样的，这意味着空间单元必须可以被重复观测。例如，在研究环境污染的影响因素时，我们可以收集每个国家在每年的 $CO_2$ 排放量，得到**平衡的时空面板**。

但是，有些情况下一个具体的空间单位在一段时间内仅能被观察一次，比如房屋交易、事故、犯罪、公司的成立/关闭等。这种情况下不同的时点上观测到的空间单元是不同的，从而得到混合空间截面数据。下图直观地展示了在不同离散的时点上收集到的不同空间单元的观测值。

空间自相关的特殊性在于其多向作用 (Multidirectional Effect)，某一空间单元受到周围其它空间单元的影响，同时它也对周围其它空间单元产生影响，下图中的红色双向箭头描述了**空间维度的多向效应**。时间的自相关是单向性的，昨天的观测值可以影响到今天的观测值，但反之则不可，也就是时间惯性在方向上的单一性，下图中灰色单向箭头描述了**时间维度的单向效应** (Unidirectional Effect)。

如果研究的数据集是混合截面数据，或者是由于缺失值导致的非平衡空间面板数据，这种情况下可以把上述数据集看做成一个混合的大截面，此时时空权重矩阵的维度不再是空间单元的个数，而是观测值的数目。

![](https://fig-lianxh.oss-cn-shenzhen.aliyuncs.com/%E8%A2%81%E5%AD%90%E6%99%B4_STAR_Fig01.png)

> 图片来源：Dubé, J., & Legros, D. (2013); Dubé, J., & Legros, D. (2018)

## 2. 时空自回归模型

### 2.1 STAR 模型设定

因此运用混合空间截面数据研究时空自回归模型 ( Spatio-Temporal Autoregressive Regression, STAR ) 时，既要考虑到空间维度的多向效应也要兼顾时间维度的单向效应：
$$
y_{i t}=\rho y_{j t}+\psi y_{j t-1}+\epsilon_{i t}
\quad(1)
$$
其中，空间滞后项为：
$$
y_{j t}=\sum_{j=1}^{N} w_{i j}^{\star} y_{j t}
\quad(2)
$$
时间滞后项为
$$
y_{j t-1}=L \cdot y_{j t}
\quad(3)
$$

以矩阵形式可以写成：

$$
\mathbf{y}=\rho \underline{\mathbf{S}} \mathbf{y}+\psi \underline{\mathbf{W}} \mathbf{y}+\epsilon
\quad(4)
$$
其中  $\underline{\mathbf{S}}$  是衡量**当期**不同观测值之间关系的空间权重矩阵（空间维度的多向效应）； $\underline{\mathbf{W}}$  是衡量  $(t-1)$ 时期与  $t$  时期观测值之间关系的时空权重矩阵（时间维度的单向效应 ），这两个矩阵的维度均为 $\left(N_{T} \times N_{T}\right)$  。系数  $\rho$  和 $\psi$ 分别衡量同期的空间依赖程度——**空间溢出效应** ( Spatial Spillover Effect) 和与滞后一期的空间依赖程度——**动态空间效应** (Dynamic Spatial Effect)。需要注意的一点是时间滞后项不包含自身的滞后 $y_{i t-1}$ ,这是由混合空间截面数据的特点决定的。

### 2.2 时空权重矩阵 (Spatio-Temporal Weights Matrix )

一般地，将空间权重 $s_{i j}$  和 时间权重  $t_{i j}$ 相结合可以得到时空权重：


$$
\mathbf{W}=\mathbf{S} \odot \mathbf{T}=\left(\begin{array}{ccccc}
0 & s_{12} \times t_{12} & s_{13} \times t_{13} & \cdots & s_{1 N} \times t_{1 N} \\
s_{21} \times t_{21} & 0 & s_{23} \times t_{23} & \cdots & s_{2 N} \times t_{2 N} \\
s_{31} \times t_{31} & s_{32} \times t_{32} & & \cdots & s_{3 N} \times t_{3 N} \\
\vdots & \vdots & \vdots & \ddots & \vdots \\
s_{N 1} \times t_{N 1} & s_{N 2} \times t_{N 2} & s_{N 3} \times t_{N 3} & \cdots & s_{N N} \times t_{N N}
\end{array}\right)
\quad(5)
$$


- **第一种情况**：假设所有的观测值都是同时发生， 即时间权重  $t_{i j}$ 均为 1时，时空权重矩阵退化为空间权重矩阵 $(\mathbf{W}=\mathbf{S})$
- **第二种情况**：假设空间自相关只存在于同一时期内，不考虑跨期的空间自相关，此时的时空权重矩阵除了对角线外元素均为0，其中 $\mathbf{S}_{tt}$ 是描述在时点 $t$ 上观测值 $i$ 和 $j$ 之间的空间关联的空间权重矩阵

$$
\underline{\mathbf{S}}=\left(\begin{array}{ccccc}
\mathbf{S}_{11} & 0 & 0 & \cdots & 0 \\
0 & \mathbf{S}_{22} & 0 & \cdots & 0 \\
0 & 0 & \mathbf{S}_{33} & \cdots & 0 \\
\vdots & \vdots & \vdots & \ddots & \vdots \\
0 & 0 & 0 & \cdots & \mathbf{S}_{T T}
\end{array}\right)
\quad(6)
$$


- **第三种情况**：只考虑到过去和现在的观测之间的时间惯性关系，此时的时空权重矩阵是下三角矩阵，其中  矩阵 $\mathbf{W}_{qp}$  是刻画观测值 $i$  在时点 $q$  与观测值 $j$  在之前的时点 $p$ 之间的空间权重矩阵。

$$
\underline{\mathbf{W}}=\left(\begin{array}{ccccc}
0 & 0 & 0 & \cdots & 0 \\
\mathbf{W}_{21} & 0 & 0 & \cdots & 0 \\
\mathbf{W}_{31} & \mathbf{W}_{32} & 0 & \cdots & 0 \\
\vdots & \vdots & \vdots & \ddots & \vdots \\
\mathbf{W}_{T 1} & \mathbf{W}_{T 2} & \mathbf{W}_{T 3} & \cdots & 0
\end{array}\right)
\quad(7)
$$



- **第四种情况**：与第三种情况想法，考虑到预期 (anticipation) 效应，此时的时空权重矩阵是上三角矩阵。

$$
\overline{\mathbf{W}}=\left(\begin{array}{ccccc}
0 & \mathbf{W}_{12} & \mathbf{W}_{13} & \cdots & \mathbf{W}_{1 T} \\
0 & 0 & \mathbf{W}_{23} & \cdots & \mathbf{W}_{2 T} \\
0 & 0 & 0 & \cdots & \mathbf{W}_{3 T} \\
\vdots & \vdots & \vdots & \ddots & \vdots \\
0 & 0 & 0 & \cdots & 0
\end{array}\right)
\quad(8)
$$

上面之所以列出几种特殊情况下退化的时空权重矩阵，是为了可以将一般化的时空权重矩阵拆分如下：
$$
\mathbf{W}=\left(\begin{array}{ccccc}
\mathbf{S}_{11} & \mathbf{W}_{12} & \mathbf{W}_{13} & \cdots & \mathbf{W}_{1 T} \\
\mathbf{W}_{21} & \mathbf{S}_{22} & \mathbf{W}_{23} & \cdots & \mathbf{W}_{2 T} \\
\mathbf{W}_{31} & \mathbf{W}_{32} & \mathbf{S}_{33} & \cdots & \mathbf{W}_{3 T} \\
\vdots & \vdots & \vdots & \ddots & \vdots \\
\mathbf{W}_{T 1} & \mathbf{W}_{T 2} & \mathbf{W}_{T 3} & \cdots & \mathbf{S}_{T T}
\end{array}\right)
\quad(9)
$$

$$
\mathbf{W} = \underline{\mathbf{S}} + \underline{\mathbf{W}}+ \overline{\mathbf{W}}
\quad(10)
$$


## 3. Stata 实例

以下例子来源于 [Stata do文档和示例数据下载](https://www.researchgate.net/project/Estimating-spatio-temporal-autoregressive-model-using-micro-data-pooled-over-time/update/5deef1a53843b09383955570)，是基于1991年至1996年在一个虚构城市发生的房屋交易案例，数据集中包含每笔交易的房屋坐标 (X, Y) 以及交易时间，基于随着时间和空间维度变化的面板微观数据，估计时空自回归模型  (STAR, spatio-temporal autoregressive model )，Stata 中的外部命令 `geo2xy ` 可以使用地图投影变换将地理经纬度转换为坐标。

### 3.1 数据预处理

第一步：进行数据预处理，包括主要解释变量和被解释变量的对数化处理和生成年份和组别的虚拟变量：

```stata
**Loading the data base
use "Transaction-STExample.dta"

**Generating temporal (continuous) variable
quietly generate date = 365*(year - 1990) + 31*(month - 1) + day
sort date   /*Chronologically order data*/

**Generating coordinates in km (instead of meters)
** (eliminating scale problem for building weights)
quietly generate xc = XMTM/1000
quietly generate yc = YMTM/1000

**Generating dependent variable
quietly generate log_price    = log(saleprice)

/*Generating independent variables*/

**Generating group of independent dummy variables
tab quality, gen(QIndex)
tab location, gen(Sfixed)

**Generating time fixed effects (annual dummy variables) 
quietly generate d91 = (year==1991)
quietly generate d92 = (year==1992)
quietly generate d93 = (year==1993)
quietly generate d94 = (year==1994)
quietly generate d95 = (year==1995)

**Applying mathematical transformation on independent variables
quietly generate log_livearea = log(livearea)
quietly generate log_age = log(age)
```

### 3.2 生成时空权重矩阵

虽然 Stata 的 `spmatrix create` 或  `spmatrix` 命令 (详见连享会推文 [空间权重矩阵的构建](https://www.lianxh.cn/news/919826ca0da88.html) ) 可以直接根据地理坐标创建空间权重矩阵，但是这些模块通常依赖于单一的空间关系范式：距离的倒数。此外这个命令是基于循环函数的，较为耗费时间和内存。所以 Jean Dubé 和 Diègo Legros  利用 Stata 中的 Mata 模块 ( 详见连享会推文 [Stata - Mata系列 (一)：Mata 入门](https://www.lianxh.cn/news/e23df70afde87.html))，允许研究者在 Stata 中自定义地创建并导出生成的矩阵，从而用于空间自相关指数的计算或空间和时空自回归模型的估计。事实上，在时空建模中能够灵活地生成自己所需要的空间权重矩阵特别有用。

```stata
**--------------------------------------------------
**Build spatio-temporally lagged variables
**--------------------------------------------------
quietly generate WTprice = .

mata 
	/*1.导入用于计算的列向量 Import the main vectors used for calculations*/
	N  = st_nobs()
	yn = st_data(.,"log_price")    /*Vector of dependent variables*/
	XC = st_data(.,"xc")           /*Vector of X coordinates*/
	YC = st_data(.,"yc")           /*Vector of Y coordinates*/
	TC = st_data(.,"date")         /*Vector of temporal coordinates*/
	/*2.生成空间和时间距离矩阵 Generate the spatial and temporal matrices*/
	XX = (XC*J(N,1,1)') - (J(N,1,1)*XC')  /*Distance between X coordinates*/
	YY = (YC*J(N,1,1)') - (J(N,1,1)*YC')  /*Distance between Y coordinates*/
	ZZ = (TC*J(N,1,1)') - (J(N,1,1)*TC')  /*Eucledian distance between transactions*/
	DD = sqrt((XX:*XX) + (YY:*YY))        /*Distance matrix*/
	SO = exp(-DD)-I(N)         /*Spatial weight matrix (negative exponential transformation)*/
	TT = exp(-abs(ZZ))         /*Temporal weight matrix (negative exponential transformation)*/
	/*3.对于时空关系加以限制 Impose restrictions on spatial and temporal relations*/
	DK = ((DD*J(N,1,1)):/N)*J(N,1,1)'      /*Estimate & attribute to each observtion the mean distance of all observations*/
	T1 = (ZZ:<=65):*(ZZ:>25)   /*Define the temporal "past" window (in this case between 25 and 65 days)*/
	/*4.生成空间和时空权重矩阵 Generating spatial and spatio-temporal weights*/ 
	ST = (DD:<=DK):*SO         /*Truncated spatial weights matrix*/
	TU = T1:*TT                /*Temporal weights within temporal "past" window*/
	Wb = ST:*TU                /*Lower triangular spatio-temporal weight matrix*/
	/*5.行标准化 Standardizing spatio-temporal weight matrix*/
	SL  = (Wb*J(N,1,1))*J(N,1,1)' /*Row sums*/
	WSb = Wb:/(SL + (SL:==0))     /*Row standardization*/
	/*6.生成时空滞后项 Generating the spatiotemporally dynamic depedent variable*/
	Wyn = WSb*yn			     /*Creating spatio-temporally lagged dependent variable (exogenous variable)*/
	st_store(.,("WTprice"),(Wyn))
end

**Eliminating observations having no spatio-temporal neighbours
drop if WTprice==0
clear mata
```

> **代码解析**：
>
> 1. **导入用于计算的列向量**：$N$ 是观测值总数，$yn$ 是被解释变量对数化的房价，$XC$ 和 $YC$ 分别是地理位置坐标，$TC$ 代表交易时间；
> 2. **生成空间和时间距离矩阵**：其中 $J(N,1,1)$ 代表 $N\times 1$ 且元素均为 1 的矩阵，$J(N,1,1)'$ 为该矩阵的转置，在 Mata 模块中当使用冒号 ( `:` ) 的时候，包括加 (`+`)、减 (`-`)、乘 (`*`)、除(`/`)在内的代数运算符可以被用于**元素对元素**的运算。根据二维平面上两点 $\mathbf{a} (x_i,y_i)$ 与 $\mathbf{b}(x_j,y_j)$ 间的欧氏距离 (Euclidean Distance) 公式：$d_{ij} = \sqrt{(x_i-x_j)^2 + (y_i - y_j)^2}$ 计算得到距离矩阵 $DD$ ，然后对空间距离和时间距离进行负指数转换，得到空间权重矩阵为 $SO$ 和时间权重矩阵 $TT$；
> 3. **对于时空关系加以限制**：$DK$ 计算出了所有空间距离的均值；$T1$ 定义了滞后一期的时间窗口，即时间距离（两次交易时间间隔）在25 天和 65 天之间，在区间 (25,65] 内属于 $t-1$ 期；
> 4. **生成空间和时空权重矩阵**：只保留空间权重矩阵中小于或等于上述距离均值的元素，其余为0，得到 $ST$ ；时间权重矩阵中只保留  $t-1$ 期元素，其余为0，得到 $TU$
> 5. **行标准化**：首先求出行和，然后每个元素除以所在行的行和 ( 若行和为0，则除以1 ) ，最终得到行标准化后的时空权重矩阵 $WSb$; 
> 6. **生成被解释变量的时空滞后项**：生成 $Wyn$ 即公式 (4) 中的 $\underline{\mathbf{W}} \mathbf{y}$ ，存储为变量 **WTprice** 并且将没有与任何点有时空相邻的观测值剔除。



```stata
**--------------------------------------------------
**Building spatial block-diagonal weights matrix
**--------------------------------------------------
mata 
	/*1. 导入用于计算的列向量 Import the main vectors used for calculations*/
	N  = st_nobs()
	yn = st_data(.,"log_price")    /*Vector of dependent variables*/
	XC = st_data(.,"xc")           /*Vector of X coordinates*/
	YC = st_data(.,"yc")           /*Vector of Y coordinates*/
	TC = st_data(.,"date")         /*Vector of temporal coordinates*/
	/* 2. 生成空间和时间距离矩阵 Generate the spatial and temporal matrices*/
	XX = (XC*J(N,1,1)') - (J(N,1,1)*XC')  /*Distance between X coordinates*/
	YY = (YC*J(N,1,1)') - (J(N,1,1)*YC')  /*Distance abetween Y coordinates*/
	ZZ = (TC*J(N,1,1)') - (J(N,1,1)*TC')  /*Eucledian distance between transactions*/
	DD = sqrt((XX:*XX) + (YY:*YY))        /*Distance matrix*/
	SO = exp(-DD)-I(N)         /*Spatial weight matrix (negative exponential transformation)*/
	TT = exp(-abs(ZZ))         /*Temporal weight matrix (negative exponential transformation)*/
	/* 3. 对于时空关系加以限制 Impose restrictions on spatial and temporal relations*/
	DK = ((DD*J(N,1,1)):/N)*J(N,1,1)'      /*Estimate & attribute to each observtion the mean distance of all observations*/
	TO = (ZZ:<=25):*(ZZ:>=-15) /*Define the "present" temporal window or contemporaneous period*/
	/*4. 构建当期的空间权重矩阵 Creating contemporaneous spatial weight matrix*/
	TO = TO:*TT
	Sb = SO:*TO
	/*5. 行标准化 Row-standardizing the contemporaneous weight matrix*/
	SL  = (Sb*J(N,1,1))*J(N,1,1)' /*Row sums*/
	SSb = Sb:/SL                  /*Row standardization*/
	/*6. 导出 Exporting the contemporaneous weight matrix*/
	st_matrix("Sb",Sb)
end
```

> **代码解析**：
>
> **注意**：这段代码的目的是为了生成公式 (4) 中的 $\underline{\mathbf{S}} $ ，从而衡量**当期**不同观测值之间关系的空间权重矩阵（空间维度的多向效应），也就是不需要考虑跨期的空间相关性，所以只有第 3 步与上面不同。
>
> **对于时空关系加以限制**：$DK$ 计算出了所有空间距离的均值；$T0$ 定义了**当期**的时间窗口，即时间距离（两次交易时间间隔）在 -15 天和 25 天之间，在区间 [-15,25] 内属于 $t$ 期；
>
> 最后空间权重矩阵与时间权重矩阵对应元素相乘得到时空权重矩阵 $Sb$，即公式 (4) 中的 $\underline{\mathbf{S}}$

### 3.3 STAR模型估计

将上述生成的时空权重矩阵 $Sb$ 导入 Stata 中并进行行标准化，命名为 $WST$，使用全局宏 **varX** 表示一组主要解释变量，**varT** 表示一组时间虚拟变量，**varS** 表示一组地点虚拟变量，被解释变量为 **log_price**。

```stata
**--------------------------------------------------
**Importing spatial weights matrix in Stata
**--------------------------------------------------
spmat putmatrix WST Sb, normalize(row) id(ID) replace
spmat note WST : "Spatial - block diagonal - row-standardized weights matrix"
spmat graph WST, blocks(10)  /*Graphing the spatial relations of the weights matrix*/

**--------------------------------------------------
**Generating macro for independent variables
**--------------------------------------------------

global varX log_livearea log_age basement garage panoramicview QIndex2 QIndex3 bathroom
global varT d91-d95
global varS Sfixed1 Sfixed2 Sfixed4 Sfixed5

```



Stata 官方 `spregress` 命令可以估计 STAR 模型，具体地， 采用广义空间最小二乘法对 STAR 模型进行估计 (generalized spatial two-stage least squares, GS2SLS) 估计，估计结果如下，时空滞后项 **WTprice** 前面的系数估计值为 0.0773753 且统计上高度显著，即公式 (4) 中的 $\hat{\psi}$ ，度量了动态空间溢出效应，说明该效应显著为正。输出结果中的空间自回归系数 lambda 为0.1023655 且显著为正，即公式 (4) 中的 $\hat{\rho}$ ，度量了空间依赖关系。



```stata
. spreg gs2sls log_price $varX $varT $varS WTprice, id(ID) dlmat(WST)

Spatial autoregressive model                      Number of obs   =     7084
(GS2SLS estimates)

-------------------------------------------------------------------------------
    log_price |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
--------------+----------------------------------------------------------------
log_price     |
 log_livearea |   1.128009    .008789   128.34   0.000     1.110783    1.145236
      log_age |  -.0540687   .0057479    -9.41   0.000    -.0653344    -.042803
     basement |   .0789563   .0087618     9.01   0.000     .0617834    .0961292
       garage |   .0314771   .0148944     2.11   0.035     .0022846    .0606696
panoramicview |   .1619325    .023933     6.77   0.000     .1150248    .2088403
      QIndex2 |   .0472097   .0105297     4.48   0.000     .0265718    .0678476
      QIndex3 |   .0508522   .0097422     5.22   0.000     .0317578    .0699466
     bathroom |   .1025071    .008493    12.07   0.000     .0858611    .1191531
          d91 |   .2186493   .0134815    16.22   0.000     .1922261    .2450725
          d92 |   .2192898    .012974    16.90   0.000     .1938613    .2447184
          d93 |   .1438101   .0113244    12.70   0.000     .1216147    .1660055
          d94 |   .0668429    .010911     6.13   0.000     .0454577    .0882282
          d95 |   .0516173   .0104344     4.95   0.000     .0311662    .0720684
      Sfixed1 |   .0867594   .0255836     3.39   0.001     .0366164    .1369023
      Sfixed2 |  -.0756727   .0258302    -2.93   0.003    -.1262989   -.0250466
      Sfixed4 |  -.1624837     .00875   -18.57   0.000    -.1796334   -.1453339
      Sfixed5 |  -.2889483   .0118157   -24.45   0.000    -.3121067     -.26579
      WTprice |   .0773753   .0091718     8.44   0.000      .059399    .0953516
        _cons |   5.293063   .1625815    32.56   0.000     4.974409    5.611717
--------------+----------------------------------------------------------------
lambda        |
        _cons |   .1023655   .0121711     8.41   0.000     .0785105    .1262205
-------------------------------------------------------------------------------
```



## 4. 参考文献

- Dubé, J., & Legros, D. (2013). Dealing with spatial data pooled over time in statistical models. *Letters in Spatial and Resource Sciences*, *6*(1), 1-18.[[pdf]](https://sci-hub.se/https://doi.org/10.1007/s12076-012-0082-3)

- Dubé, J., & Legros, D. (2014). *Spatial econometrics using microdata*. John Wiley & Sons. [[pdf]](https://sci-hub.se/10.1002/9781119008651)

- Dubé, J., & Legros, D. (2018). Decomposing and Interpreting Spatial Effects in Spatio-Temporal Analysis: Evidences for Spatial Data Pooled Over Time. In *GeoComputational Analysis and Modeling of Regional Systems* (pp. 373-394). Springer, Cham. [[pdf]](https://sci-hub.se/10.1007/978-3-319-59511-5_19)

- Thanos, S., Dubé, J., & Legros, D. (2016). Putting time into space: the temporal coherence of spatial applications in the housing market. *Regional Science and Urban Economics*, *58*, 78-88.[[pdf]](https://sci-hub.se/https://doi.org/10.1016/j.regsciurbeco.2016.03.001)

- [连享会推文-空间计量专题](https://www.lianxh.cn/blogs/29.html)

-  [Stata - Mata系列 (一)：Mata 入门](https://www.lianxh.cn/news/e23df70afde87.html)

  

