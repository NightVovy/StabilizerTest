import numpy as np
from itertools import combinations, product

# 定义部分纠缠态 |ψ⟩ = cosθ |00⟩ + sinθ |11⟩
theta = np.pi / 5  # θ = π/5
cos_theta = np.cos(theta)
sin_theta = np.sin(theta)
sin_2theta = np.sin(2 * theta)
cos_2theta = np.cos(2 * theta)

# 定义10个算子及其作用结果（|00⟩, |01⟩, |10⟩, |11⟩ 的系数）
operators = {
    "I": [cos_theta, 0, 0, sin_theta],      # I |ψ⟩
    "XA": [0, sin_theta, cos_theta, 0],     # X_A |ψ⟩
    "XB": [0, cos_theta, sin_theta, 0],     # X_B |ψ⟩
    "ZA": [cos_theta, 0, 0, -sin_theta],    # Z_A |ψ⟩
    "ZB": [cos_theta, 0, 0, -sin_theta],    # Z_B |ψ⟩
    "XAXB": [sin_theta, 0, 0, cos_theta],   # X_A X_B |ψ⟩
    "XAZB": [0, -sin_theta, cos_theta, 0],  # X_A Z_B |ψ⟩
    "ZAXB": [0, cos_theta, -sin_theta, 0],  # Z_A X_B |ψ⟩
    "YAYB": [-sin_theta, 0, 0, -cos_theta], # Y_A Y_B |ψ⟩
    "ZAZB": [cos_theta, 0, 0, sin_theta]    # Z_A Z_B |ψ⟩
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

# 遍历所有5个算子的组合（C(10,5)=252种）
found = False
for op_combo in combinations(operators.keys(), 3):
    S1, S2, S3 = [operators[op] for op in op_combo]

    # 遍历所有可能的系数组合（a, b, c），最多1个为0
    for a, b, c in product(coefficients, repeat=3):  # 改为3个系数
        if [a, b, c].count(0) > 1:  # 跳过多于1个系数为0的情况
            continue

        # 计算 a*S1 + b*S2 + c*S3
        result = [
                a * S1[0] + b * S2[0] + c * S3[0],  # |00⟩
                a * S1[1] + b * S2[1] + c * S3[1],  # |01⟩
                a * S1[2] + b * S2[2] + c * S3[2],  # |10⟩
                a * S1[3] + b * S2[3] + c * S3[3],  # |11⟩
        ]

        # 检查是否满足nullifier条件（结果为0向量）
        target = [0, 0, 0, 0]
        if np.allclose(result, target, atol=1e-8):
            # 构建解的唯一标识（忽略系数为0的算子，并排序）
            active_ops_coeffs = [(op, coeff) for op, coeff in zip(op_combo, [a, b, c])
                                 if not np.isclose(coeff, 0, atol=1e-8)]
            # 按算子名称和系数排序
            solution_key = tuple(sorted(active_ops_coeffs, key=lambda x: x[0]))

            if solution_key not in unique_solutions:
                unique_solutions.add(solution_key)
                active_ops = [op for op, _ in solution_key]
                active_coeffs = [coeff for _, coeff in solution_key]

                print(f"Found nullifier combination: {op_combo}")
                print(f"Coefficients (a, b, c): {a}, {b}, {c}")
                print(f"Coefficient types: {coeff_names[a]}, {coeff_names[b]}, {coeff_names[c]}")
                print(f"Active operators: {active_ops} with coefficients {active_coeffs}")
                print("Verification:")
                print("Combination result:", np.round(result, 8))
                print("Target:            ", np.round(target, 8))
                print("---")
                found = True

if not found:
    print("No valid nullifier combination found.")
else:
    print(f"Total unique solutions found: {len(unique_solutions)}")

    # 新增部分：输出所有解的汇总
    print("\n=== 所有Nullifier组合的汇总 ===")
    for i, solution in enumerate(sorted(unique_solutions), 1):
        # 解构解决方案
        ops = [item[0] for item in solution]
        coeffs = [item[1] for item in solution]
        coeff_names_list = [coeff_names.get(c, str(c)) for c in coeffs]

        # 格式化输出
        print(f"\n解 {i}:")
        print("算子组合:", " + ".join(f"{coeff_names_list[i]}*{ops[i]}"
                                      for i in range(len(ops))))
        print("详细系数:")
        for op, coeff in zip(ops, coeffs):
            print(f"  {op}: {coeff:.6f} ({coeff_names.get(coeff, str(coeff))})")

    print("\n=== 汇总结束 ===")