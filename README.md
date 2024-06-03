# FerroDispCalc

Calculator of polarization displacements and polarization in ferroelectric materials.

1. 代码分成三部分：构造临近列表、计算、后处理

    * 对于构造临近列表：使用python完成，读取单个结构；可以用joblib并行，然后输出文件
    * 对于计算：用python or cpp完成，python用joblib并行，CPP用openmp并行。
    * 后处理：支持自动分层，氧空位，晶界可视化，丰富的可视化（python or cpp）

2. 支持的体系：钙钛矿（固溶体、空位、晶界）、二氧化铪