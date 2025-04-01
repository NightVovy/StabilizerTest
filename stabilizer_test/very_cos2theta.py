import numpy as np
from qiskit.quantum_info import Operator, Statevector

# 定义Pauli矩阵
I = np.array([[1, 0], [0, 1]])  # 单位矩阵
X = np.array([[0, 1], [1, 0]])  # Pauli X
Y = np.array([[0, -1j], [1j, 0]])  # Pauli Y
Z = np.array([[1, 0], [0, -1]])  # Pauli Z

# 定义双量子比特Pauli算子
XX = np.kron(X, X)
ZZ = np.kron(Z, Z)


# 定义态 |ψ⟩ = cosθ|00⟩ + sinθ|11⟩
def create_state(theta):
    """创建态 cosθ|00⟩ + sinθ|11⟩"""
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)
    return np.array([cos_theta, 0, 0, sin_theta])


# 创建组合算子 sin(2θ)XX + cos(2θ)ZZ
def create_combined_operator(theta):
    """创建算子 sin(2θ)XX + cos(2θ)ZZ"""
    sin_2theta = np.sin(2 * theta)
    cos_2theta = np.cos(2 * theta)
    return sin_2theta * XX + cos_2theta * ZZ


# 验证稳定子条件
def verify_stabilizers(theta):
    """
    验证:
    1. ZZ|ψ⟩ = +1|ψ⟩
    2. [sin(2θ)XX + cos(2θ)ZZ]|ψ⟩ = +1|ψ⟩
    返回验证结果和计算的本征值
    """
    # 创建态和算子
    psi = create_state(theta)
    combined_op = create_combined_operator(theta)

    # 验证ZZ
    ZZ_psi = np.dot(ZZ, psi)
    ZZ_eigenvalue = np.vdot(psi, ZZ_psi)
    is_ZZ_stabilizer = np.isclose(ZZ_eigenvalue, 1.0, atol=1e-10)

    # 验证组合算子
    combined_psi = np.dot(combined_op, psi)
    combined_eigenvalue = np.vdot(psi, combined_psi)

    # 检查是否为比例关系
    # 计算两个态之间的保真度
    fidelity = np.abs(np.vdot(psi, combined_psi / np.linalg.norm(combined_psi)))
    is_combined_stabilizer = np.isclose(fidelity, 1.0, atol=1e-10)

    return {
        'ZZ': {
            'is_stabilizer': is_ZZ_stabilizer,
            'eigenvalue': ZZ_eigenvalue
        },
        'combined': {
            'is_stabilizer': is_combined_stabilizer,
            'eigenvalue': combined_eigenvalue,
            'fidelity': fidelity
        }
    }


# 测试不同的theta值
test_angles = [0, np.pi / 6, np.pi / 4, np.pi / 3, np.pi / 2]
print("验证稳定子:")
print("{:<10} {:<25} {:<15} {:<25} {:<15}".format(
    "θ", "ZZ本征值", "ZZ是否稳定子", "组合算子本征值", "组合算子是否稳定子"
))

for theta in test_angles:
    results = verify_stabilizers(theta)
    print("{:<10.4f} {:<25.10f} {:<15} {:<25.10f} {:<15}".format(
        theta,
        results['ZZ']['eigenvalue'].real,  # 本征值应为实数
        str(results['ZZ']['is_stabilizer']),
        results['combined']['eigenvalue'].real,
        str(results['combined']['is_stabilizer'])
    ))

# 详细分析一个特定角度
theta = np.pi / 3
print(f"\n详细分析 θ = {theta:.4f}:")
psi = create_state(theta)
combined_op = create_combined_operator(theta)

print("\n态 |ψ⟩:")
print(psi)

print("\n组合算子 sin(2θ)XX + cos(2θ)ZZ:")
print(combined_op)

# 应用算子
result = np.dot(combined_op, psi)
print("\n应用算子后的结果:")
print(result)

print("\n原始态与结果态的比例关系:")
# 计算比例因子
scale_factor = result[0] / psi[0] if psi[0] != 0 else result[3] / psi[3]
print(f"比例因子: {scale_factor}")

# 使用Qiskit验证
print("\n使用Qiskit验证:")
psi_sv = Statevector(psi)
combined_op_qiskit = Operator(combined_op)

# 计算期望值
expectation_value = psi_sv.expectation_value(combined_op_qiskit)
print(f"组合算子的期望值: {expectation_value:.10f}")
print(f"态是否为组合算子的本征态: {np.isclose(abs(expectation_value), 1.0, atol=1e-10)}")