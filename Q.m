function y=Q(x)
%% 输入：x需要计算的序列
%% 输出：标准误差函数
y=erfc(x/sqrt(2))/2;