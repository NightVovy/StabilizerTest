import numpy as np
from qiskit.quantum_info import Operator, Statevector

# 定义Pauli矩阵
I = np.array([[1, 0], [0, 1]])  # 单位矩阵
X = np.array([[0, 1], [1, 0]])  # Pauli X
Y = np.array([[0, -1j], [1j, 0]])  # Pauli Y
Z = np.array([[1, 0], [0, -1]])  # Pauli Z

# 定义ZZ算子 (张量积)
ZZ = np.kron(Z, Z)

# 定义态 |ψ⟩ = cosθ|00⟩ + sinθ|11⟩
def create_state(theta):
    """创建态 cosθ|00⟩ + sinθ|11⟩"""
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)
    return np.array([cos_theta, 0, 0, sin_theta])

# 验证稳定子条件
def verify_stabilizer(theta):
    """
    验证 ZZ|ψ⟩ = +1|ψ⟩
    返回验证结果和计算的本征值
    """
    # 创建态
    psi = create_state(theta)
    
    # 应用ZZ算子
    ZZ_psi = np.dot(ZZ, psi)
    
    # 检查是否为比例关系
    # 由于可能存在全局相位，我们比较两个态是否平行
    # 通过计算内积 |⟨ψ|ZZ|ψ⟩| 应该等于1
    eigenvalue = np.vdot(psi, ZZ_psi)  # ⟨ψ|ZZ|ψ⟩
    
    # 检查是否近似等于1 (考虑浮点误差)
    is_stabilizer = np.isclose(eigenvalue, 1.0, atol=1e-10)
    
    return is_stabilizer, eigenvalue

# 测试不同的theta值
test_angles = [0, np.pi/6, np.pi/4, np.pi/3, np.pi/2]
print("验证稳定子 ZZ 对态 cosθ|00⟩ + sinθ|11⟩ 的作用:")
print("{:<10} {:<25} {:<15}".format("θ", "本征值", "是否稳定子"))

for theta in test_angles:
    is_stab, eigval = verify_stabilizer(theta)
    print("{:<10.4f} {:<25} {:<15}".format(
        theta, 
        f"{eigval:.10f}", 
        str(is_stab)
    ))

# 使用Qiskit验证
print("\n使用Qiskit验证:")
theta = np.pi/3  # 示例角度
psi = Statevector(create_state(theta))
zz_op = Operator(np.kron(Z, Z))

# 计算期望值
expectation_value = psi.expectation_value(zz_op)
print(f"对于θ={theta:.4f}, ZZ的期望值: {expectation_value:.10f}")
print(f"态是否为ZZ的本征态: {np.isclose(abs(expectation_value), 1.0, atol=1e-10)}")