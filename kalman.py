# -*- coding: utf-8 -*-
"""
Created on 2023年10月25日

作者: jason

邮箱: zhengjie9510@qq.com

代码说明: 
"""
import numpy as np


def kalman_filter(x: np.ndarray, z: np.ndarray, u: np.ndarray, A: np.ndarray, B: np.ndarray, H: np.ndarray, Q: np.ndarray, R: np.ndarray, x_hat_first: np.ndarray, P_first: np.ndarray) -> np.ndarray:
    """
    卡尔曼滤波器

    Args:
        x: 状态量

        z: 观测值

        u: 控制量

        A: 状态转移矩阵

        B: 控制矩阵

        H: 观测矩阵

        Q: 系统噪声协方差矩阵

        R: 测量噪声协方差矩阵

        x_hat_first: 初始估计

        P_first: 初始估计协方差

    Returns:
        x_hat: 后验估计

    公式:
        预测:
            x_hat_minus[k]= A * x_hat[k-1] + B * u[k-1]

            P_minus[k] = A * P[k-1] * A.T + Q

        更新:
            K[k] = P_minus[k] * H.T * (H * P_minus[k] * H.T + R).I

            x_hat[k] = x_hat_minus[k] + K[k] * (z[k] - H * x_hat_minus[k])

            P[k] = (np.eye(1) - K[k] * H) * P_minus[k]
    """
    rows, columns = x.shape
    x_hat_minus = np.zeros((rows, columns))  # 先验估计
    x_hat = np.empty((rows, columns))  # 后验估计
    P_minus = np.empty((rows, columns, columns))  # 先验估计协方差
    P = np.empty((rows, columns, columns))  # 后验估计协方差
    K = np.empty((rows, columns, columns))  # 卡尔曼增益

    x_hat[0] = x_hat_first  # 初始估计
    P[0] = P_first  # 初始估计协方差

    for k in range(1, rows):
        # 预测
        x_hat_minus[k] = np.dot(A, x_hat[k-1]) + np.dot(B, u[k-1])
        P_minus[k] = np.dot(A, np.dot(P[k-1], np.transpose(A))) + Q
        # 更新
        K[k] = np.dot(P_minus[k], np.dot(np.transpose(H), np.linalg.inv(np.dot(H, np.dot(P_minus[k], np.transpose(H))) + R)))
        x_hat[k] = x_hat_minus[k] + np.dot(K[k], (z[k] - np.dot(H, x_hat_minus[k])))
        P[k] = np.dot((np.eye(columns) - np.dot(K[k], H)), P_minus[k])

    return x_hat
