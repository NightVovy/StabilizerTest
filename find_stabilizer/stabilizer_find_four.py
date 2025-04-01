import numpy as np
from itertools import combinations, product

# 定义部分纠缠态 |ψ⟩ = cosθ |00⟩ + sinθ |11⟩
theta = np.pi / 5  # θ = π/5
cos_theta = np.cos(theta)
sin_theta = np.sin(theta)
sin_2theta = np.sin(2 * theta)
cos_2theta = np.cos(2 * theta)

# 输出四个值
print("cos(theta) =", cos_theta)
print("sin(theta) =", sin_theta)
print("sin(2*theta) =", sin_2theta)
print("cos(2*theta) =", cos_2theta)

# 定义8个算子及其作用结果（|00⟩, |01⟩, |10⟩, |11⟩ 的系数）
operators = {
    "XA": [0, sin_theta, cos_theta, 0],  # X_A |ψ⟩
    "XB": [0, cos_theta, sin_theta, 0],  # X_B |ψ⟩
    "ZA": [cos_theta, 0, 0, -sin_theta],  # Z_A |ψ⟩
    "ZB": [cos_theta, 0, 0, -sin_theta],  # Z_B |ψ⟩
    "XAXB": [sin_theta, 0, 0, cos_theta],  # X_A X_B |ψ⟩
    "XAZB": [0, -sin_theta, cos_theta, 0],  # X_A Z_B |ψ⟩
    "ZAXB": [0, cos_theta, -sin_theta, 0],  # Z_A X_B |ψ⟩
    "YAYB": [-sin_theta, 0, 0, -cos_theta]  # Y_A Y_B |ψ⟩
}

# 定义系数名称映射
coeff_names = {
    1: "1",
    cos_theta: "cosθ",
    sin_theta: "sinθ",
    cos_2theta: "cos2θ",
    sin_2theta: "sin2θ",
    -1: "-1",
    -cos_theta: "-cosθ",
    -sin_theta: "-sinθ",
    -cos_2theta: "-cos2θ",
    -sin_2theta: "-sin2θ",
    0: "0"
}

# 定义所有可能的系数（共11种）
coefficients = [
    1, cos_theta, sin_theta, cos_2theta, sin_2theta,
    -1, -cos_theta, -sin_theta, -cos_2theta, -sin_2theta,
    0
]

# 存储唯一解
unique_solutions = set()

# 遍历所有4个算子的组合（C(8,4)=70种）
found = False
for op_combo in combinations(operators.keys(), 4):
    S1, S2, S3, S4 = [operators[op] for op in op_combo]

    # 遍历所有可能的系数组合（a, b, c, d），最多2个为0
    for a, b, c, d in product(coefficients, repeat=4):
        if [a, b, c, d].count(0) > 2:  # 跳过多于2个系数为0的情况
            continue

        # 计算 a*S1 + b*S2 + c*S3 + d*S4
        result = [
            a * S1[0] + b * S2[0] + c * S3[0] + d * S4[0],  # |00⟩
            a * S1[1] + b * S2[1] + c * S3[1] + d * S4[1],  # |01⟩
            a * S1[2] + b * S2[2] + c * S3[2] + d * S4[2],  # |10⟩
            a * S1[3] + b * S2[3] + c * S3[3] + d * S4[3],  # |11⟩
        ]

        # 检查是否满足稳定子条件
        target = [cos_theta, 0, 0, sin_theta]
        if np.allclose(result, target, atol=1e-8):
            # 构建解的唯一标识（忽略系数为0的算子，并排序）
            active_ops_coeffs = [(op, coeff) for op, coeff in zip(op_combo, [a, b, c, d])
                                 if not np.isclose(coeff, 0, atol=1e-8)]
            # 按算子名称和系数排序
            solution_key = tuple(sorted(active_ops_coeffs, key=lambda x: x[0]))

            if solution_key not in unique_solutions:
                unique_solutions.add(solution_key)
                active_ops = [op for op, _ in solution_key]
                active_coeffs = [coeff for _, coeff in solution_key]

                print(f"Found stabilizer combination: {op_combo}")
                print(f"Coefficients (a, b, c, d): {a}, {b}, {c}, {d}")
                print(f"Coefficient types: {coeff_names[a]}, {coeff_names[b]}, {coeff_names[c]}, {coeff_names[d]}")
                print(f"Active operators: {active_ops} with coefficients {active_coeffs}")
                print("Verification:")
                print("Combination result:", np.round(result, 8))
                print("Target:            ", np.round(target, 8))
                print("---")
                found = True

if not found:
    print("No valid stabilizer combination found.")
else:
    print(f"Total unique solutions found: {len(unique_solutions)}")