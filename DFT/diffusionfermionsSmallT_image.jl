using Images

# 创建一个随机的 1D 数组，模拟一张灰度图像
img = rand(256)

# 将数组转换为 2D 矩阵，大小为 16x16
img_mat = reshape(img, 16, 16)

# 定义高斯核的标准差
sigma = 2.0

# 构造高斯核
kernel = Gaussian((5, 5), sigma)

# 使用高斯核进行卷积
blurred_mat = conv2(img_mat, kernel, padtype=Symmetric)

# 将模糊后的图像矩阵转换回一维数组形式
blurred_img = reshape(blurred_mat, 256)

# 显示原始图像和模糊后的图像
using Plots
plot(heatmap(img_mat), heatmap(blurred_mat))
