/stabilizer_test
  尝试验证sin2thetaXX+cos2thetaZZ是不是量子态costheta00 + sintheta11的稳定子
    ZZ|ψ⟩ = ZZ(cosθ|00⟩ + sinθ|11⟩) = cosθ|00⟩ + sinθ|11⟩ = +1|ψ⟩
    XX|ψ⟩ = XX(cosθ|00⟩ + sinθ|11⟩) = cosθ|11⟩ + sinθ|00⟩
    [sin(2θ)XX + cos(2θ)ZZ]|ψ⟩ = sin(2θ)(cosθ|11⟩ + sinθ|00⟩) + cos(2θ)(cosθ|00⟩ + sinθ|11⟩)
                           = [sin(2θ)sinθ + cos(2θ)cosθ]|00⟩ + [sin(2θ)cosθ + cos(2θ)sinθ]|11⟩
                           = cos(2θ - θ)|00⟩ + sin(2θ + θ)|11⟩  (使用三角恒等式)
                           = cosθ|00⟩ + sin3θ|11⟩
                           = +1|ψ⟩? 只有当 θ 为 π/4 的整数倍时，组合算子才是严格的稳定子


/find_stabilizer
    因为暴力解太麻烦，所以改成代码实现。
    根本上XA（等7+1项，单独的是ZAZB）作用于部分纠缠态上就是个系数问题，组合后00部分系数是costheta,11部分是sintheta,
其他部分系数是0，这样就组合出了稳定子。
    stabilizer_find_three.py: 用3个组合，但是直接用了sin2thetaXAZA之类的。
    find_ver2.py:只用了7项而不是7x3x2项，把sin2theta之类的系数放到组合里面了。
    stabilizer_find_four.py:4个组合，然后加入了YAYB=-XAZAXBZB. 最多2个0.