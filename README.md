# Kalman-Filter-Example
# 理论公式
详细理论可参考DR_CAN关于卡尔曼滤波器的视频讲解。[https://www.bilibili.com/video/BV1dV411B7ME](https://www.bilibili.com/video/BV1dV411B7ME/?spm_id_from=333.788&vd_source=1d326cfdbea18a88b4c83bff15116a78)  
卡尔曼滤波公式分为**预测**和**更新**两部分。  
## 预测公式为：

    x_hat_minus[k]= A * x_hat[k-1] + B * u[k-1]

    P_minus[k] = A * P[k-1] * A.T + Q

## 更新公式为：

    K[k] = P_minus[k] * H.T * (H * P_minus[k] * H.T + R).I

    x_hat[k] = x_hat_minus[k] + K[k] * (z[k] - H * x_hat_minus[k])

    P[k] = (np.eye(1) - K[k] * H) * P_minus[k]

## 其中:  
    x_hat_minus--先验估计状态  
    P_minus--先验估计协方差矩阵  
    x_hat--后验估计状态  
    P--后验估计协方差矩阵  
    A--状态转移矩阵  
    B--控制矩阵  
    u--控制量  
    Q--过程噪声协方差矩阵  
    K--卡尔曼增益  
    H--观测矩阵  
    R--观测噪声协方差矩阵  
    z--观测量。