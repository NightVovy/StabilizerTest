import numpy as np
from itertools import combinations

# 定义部分纠缠态 |ψ⟩ = cosθ |00⟩ + sinθ |11⟩
theta = np.pi / 5  # 示例角度，可以修改
cos_theta = np.cos(theta)
sin_theta = np.sin(theta)
sin_2theta = np.sin(2 * theta)
cos_2theta = np.cos(2 * theta)

# 输出四个值
print("cos(theta) =", cos_theta)
print("sin(theta) =", sin_theta)
print("sin(2*theta) =", sin_2theta)
print("cos(2*theta) =", cos_2theta)

# 定义所有21个算子及其作用结果（|00⟩, |01⟩, |10⟩, |11⟩ 的系数）
operators = {
    "XA": [0, sin_theta, cos_theta, 0],  # X_A |ψ⟩
    "XA_sin": [0, sin_2theta * sin_theta, sin_2theta * cos_theta, 0],  # sin2θ * X_A |ψ⟩
    "XA_cos": [0, cos_2theta * sin_theta, cos_2theta * cos_theta, 0],  # cos2θ * X_A |ψ⟩
    "XB": [0, cos_theta, sin_theta, 0],  # X_B |ψ⟩
    "XB_sin": [0, sin_2theta * cos_theta, sin_2theta * sin_theta, 0],  # sin2θ * X_B |ψ⟩
    "XB_cos": [0, cos_2theta * cos_theta, cos_2theta * sin_theta, 0],  # cos2θ * X_B |ψ⟩
    "ZA": [cos_theta, 0, 0, -sin_theta],  # Z_A |ψ⟩
    "ZA_sin": [sin_2theta * cos_theta, 0, 0, -sin_2theta * sin_theta],  # sin2θ * Z_A |ψ⟩
    "ZA_cos": [cos_2theta * cos_theta, 0, 0, -cos_2theta * sin_theta],  # cos2θ * Z_A |ψ⟩
    "ZB": [cos_theta, 0, 0, -sin_theta],  # Z_B |ψ⟩
    "ZB_sin": [sin_2theta * cos_theta, 0, 0, -sin_2theta * sin_theta],  # sin2θ * Z_B |ψ⟩
    "ZB_cos": [cos_2theta * cos_theta, 0, 0, -cos_2theta * sin_theta],  # cos2θ * Z_B |ψ⟩
    "XAXB": [sin_theta, 0, 0, cos_theta],  # X_A X_B |ψ⟩
    "XAXB_sin": [sin_2theta * sin_theta, 0, 0, sin_2theta * cos_theta],  # sin2θ * X_A X_B |ψ⟩
    "XAXB_cos": [cos_2theta * sin_theta, 0, 0, cos_2theta * cos_theta],  # cos2θ * X_A X_B |ψ⟩
    "XAZB": [0, -sin_theta, cos_theta, 0],  # X_A Z_B |ψ⟩
    "XAZB_sin": [0, -sin_2theta * sin_theta, sin_2theta * cos_theta, 0],  # sin2θ * X_A Z_B |ψ⟩
    "XAZB_cos": [0, -cos_2theta * sin_theta, cos_2theta * cos_theta, 0],  # cos2θ * X_A Z_B |ψ⟩
    "ZAXB": [0, cos_theta, -sin_theta, 0],  # Z_A X_B |ψ⟩
    "ZAXB_sin": [0, sin_2theta * cos_theta, -sin_2theta * sin_theta, 0],  # sin2θ * Z_A X_B |ψ⟩
    "ZAXB_cos": [0, cos_2theta * cos_theta, -cos_2theta * sin_theta, 0],  # cos2θ * Z_A X_B |ψ⟩
}

# 遍历所有3项组合
found = False
for combo in combinations(operators.keys(), 3):
    S1, S2, S3 = operators[combo[0]], operators[combo[1]], operators[combo[2]]

    # 构造线性方程组：a*S1 + b*S2 + c*S3 = [cosθ, 0, 0, sinθ]
    A = np.array([
        [S1[0], S2[0], S3[0]],  # |00⟩ 系数
        [S1[1], S2[1], S3[1]],  # |01⟩ 系数
        [S1[2], S2[2], S3[2]],  # |10⟩ 系数
        [S1[3], S2[3], S3[3]],  # |11⟩ 系数
    ])
    b = np.array([cos_theta, 0, 0, sin_theta])  # 目标

    # 解线性方程组 Ax = b
    try:
        x, residuals, _, _ = np.linalg.lstsq(A, b, rcond=None)
        if np.allclose(A @ x, b, atol=1e-8):  # 检查是否满足条件
            print(f"Found stabilizer combination: {combo}")
            print(f"Coefficients (a, b, c): {x}")
            print("Verification:")
            print("a*S1 + b*S2 + c*S3 =", x[0] * np.array(S1) + x[1] * np.array(S2) + x[2] * np.array(S3))
            print("Target:             ", b)
            print("---")
            found = True
    except:
        continue

if not found:
    print("No valid stabilizer combination found.")