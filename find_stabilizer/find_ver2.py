import numpy as np
from itertools import combinations, product

# 初始化参数
theta = np.pi / 5
cos_theta = np.cos(theta)
sin_theta = np.sin(theta)
cos_2theta = np.cos(2 * theta)
sin_2theta = np.sin(2 * theta)

# 输出四个值
print("cos(theta) =", cos_theta)
print("sin(theta) =", sin_theta)
print("sin(2*theta) =", sin_2theta)
print("cos(2*theta) =", cos_2theta)

# 算子定义
operators = {
    "XA": [0, sin_theta, cos_theta, 0],
    "XB": [0, cos_theta, sin_theta, 0],
    "ZA": [cos_theta, 0, 0, -sin_theta],
    "ZB": [cos_theta, 0, 0, -sin_theta],
    "XAXB": [sin_theta, 0, 0, cos_theta],
    "XAZB": [0, -sin_theta, cos_theta, 0],
    "ZAXB": [0, cos_theta, -sin_theta, 0]
}

# 系数定义
coeff_options = [1, cos_theta, sin_theta, cos_2theta, sin_2theta,
                 -1, -cos_theta, -sin_theta, -cos_2theta, -sin_2theta, 0]

# 系数名称映射
coeff_names = {
    1: "1", -1: "-1",
    cos_theta: "cosθ", -cos_theta: "-cosθ",
    sin_theta: "sinθ", -sin_theta: "-sinθ",
    cos_2theta: "cos2θ", -cos_2theta: "-cos2θ",
    sin_2theta: "sin2θ", -sin_2theta: "-sin2θ",
    0: "0"
}

# 存储已找到的解（用非零算子组合作为键）
solutions = {}

# 遍历所有组合
for ops in combinations(operators.keys(), 3):
    for a, b, c in product(coeff_options, repeat=3):
        if [a, b, c].count(0) > 1:
            continue

        # 计算结果
        res = [
            a * operators[ops[0]][0] + b * operators[ops[1]][0] + c * operators[ops[2]][0],
            a * operators[ops[0]][1] + b * operators[ops[1]][1] + c * operators[ops[2]][1],
            a * operators[ops[0]][2] + b * operators[ops[1]][2] + c * operators[ops[2]][2],
            a * operators[ops[0]][3] + b * operators[ops[1]][3] + c * operators[ops[2]][3]
        ]

        # 检查是否满足条件
        target = [cos_theta, 0, 0, sin_theta]
        if np.allclose(res, target, atol=1e-8):
            # 构建解的唯一标识（忽略系数为0的算子）
            active_ops = []
            active_coeffs = []
            for op, coeff in zip(ops, [a, b, c]):
                if not np.isclose(coeff, 0, atol=1e-8):
                    active_ops.append(op)
                    active_coeffs.append(coeff)

            # 按字母顺序排序以确保唯一性
            active_ops_tuple = tuple(sorted(active_ops))
            active_coeffs_tuple = tuple(active_coeffs)

            # 如果是新解则存储并输出
            if active_ops_tuple not in solutions:
                solutions[active_ops_tuple] = active_coeffs_tuple
                print(f"Found stabilizer combination: {ops}")
                print(f"Coefficients (a, b, c): {a}, {b}, {c}")
                print(f"Coefficient types: {coeff_names[a]}, {coeff_names[b]}, {coeff_names[c]}")
                print(f"Active operators: {active_ops} with coefficients {active_coeffs}")
                print("Verification:")
                print("Combination result:", res)
                print("Target:           ", target)
                print("---")

if not solutions:
    print("No valid stabilizer combination found.")
else:
    print(f"Found {len(solutions)} unique solutions.")