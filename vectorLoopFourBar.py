import math
import positionSynthesis

def basic(a, b, c, d, theta2):
    theta2 = theta2 * 2 * math.pi / 360
    K1 = d / a
    K2 = d / c
    K3 = (a ** 2 - b ** 2 + c ** 2 + d ** 2) / (2 * a * c)
    K4 = d / b
    K5 = (c ** 2 - d ** 2 - a ** 2 - b ** 2) / (2 * a * b)

    A = math.cos(theta2) - K1 - K2 * math.cos(theta2) + K3
    B = -2 * math.sin(theta2)
    C = K1 - (K2 + 1) * math.cos(theta2) + K3
    D = math.cos(theta2) - K1 + K4 * math.cos(theta2) + K5
    E = -2 * math.sin(theta2)
    F = K1 + (K4 - 1) * math.cos(theta2) + K5

    theta4_1 = 2 * math.atan((-B + math.sqrt(B ** 2 - 4 * A * C)) / (2 * A)) * 360 / (2 * math.pi)
    theta4_2 = 2 * math.atan((-B - math.sqrt(B ** 2 - 4 * A * C)) / (2 * A)) * 360 / (2 * math.pi)
    theta3_1 = 2 * math.atan((-E + math.sqrt(E ** 2 - 4 * D * F)) / (2 * D)) * 360 / (2 * math.pi)
    theta3_2 = 2 * math.atan((-E - math.sqrt(E ** 2 - 4 * D * F)) / (2 * D)) * 360 / (2 * math.pi)

    print(f"Theta4 = {theta4_1} (crossed), {theta4_2} (open)")
    print(f"Theta3 = {theta3_1} (crossed), {theta3_2} (open)")
    print()

def basicAcceleration(a, b, c, d, theta2, omega2, alpha2, p, delta3):
    theta2 = math.radians(theta2)
    delta3 = math.radians(delta3)
    K1 = d / a
    K2 = d / c
    K3 = (a ** 2 - b ** 2 + c ** 2 + d ** 2) / (2 * a * c)
    K4 = d / b
    K5 = (c ** 2 - d ** 2 - a ** 2 - b ** 2) / (2 * a * b)

    A = math.cos(theta2) - K1 - K2 * math.cos(theta2) + K3
    B = -2 * math.sin(theta2)
    C = K1 - (K2 + 1) * math.cos(theta2) + K3
    D = math.cos(theta2) - K1 + K4 * math.cos(theta2) + K5
    E = -2 * math.sin(theta2)
    F = K1 + (K4 - 1) * math.cos(theta2) + K5

    theta4_1 = math.radians(2 * math.atan((-B + math.sqrt(B ** 2 - 4 * A * C)) / (2 * A)) * 360 / (2 * math.pi))
    theta4_2 = math.radians(2 * math.atan((-B - math.sqrt(B ** 2 - 4 * A * C)) / (2 * A)) * 360 / (2 * math.pi))
    theta3_1 = math.radians(2 * math.atan((-E + math.sqrt(E ** 2 - 4 * D * F)) / (2 * D)) * 360 / (2 * math.pi))
    theta3_2 = math.radians(2 * math.atan((-E - math.sqrt(E ** 2 - 4 * D * F)) / (2 * D)) * 360 / (2 * math.pi))

    print("OPEN")
    print(f"Theta3 = {round(math.degrees(theta3_2), 3)}, Theta4 = {round(math.degrees(theta4_2), 3)}")
    omega3 = a * omega2 * math.sin(theta4_2 - theta2) / (b * math.sin(theta3_2 - theta4_2))
    omega4 = a * omega2 * math.sin(theta2 - theta3_2) / (c * math.sin(theta4_2 - theta3_2))
    VaReal = a * omega2 * (-math.sin(theta2))
    VaImag = a * omega2 * math.cos(theta2)
    Va = math.sqrt(VaReal ** 2 + VaImag ** 2)
    VaAngle = positionSynthesis.invTan(VaReal, VaImag)
    VbReal = a * omega4 * (-math.sin(theta4_2))
    VbImag = a * omega4 * math.cos(theta4_2)
    Vb = math.sqrt(VbReal ** 2 + VbImag ** 2)
    VbAngle = positionSynthesis.invTan(VbReal, VbImag)
    print(f"Omega3 = {round(omega3, 3)}, Omega4 = {round(omega4, 3)}")
    print(f"Va = {round(Va, 3)}, Va angle = {round(math.degrees(VaAngle), 3)}")
    print(f"Vb = {round(Vb, 3)}, Vb angle = {round(math.degrees(VbAngle), 3)}")

    VpReal = VaReal + p * omega3 * (-math.sin(theta3_2 + delta3))
    VpImag = VaImag + p * omega3 * math.cos(theta3_2 + delta3)
    Vp = math.sqrt(VpReal ** 2 + VpImag ** 2)
    VpAngle = positionSynthesis.invTan(VpReal, VpImag)
    print(f"Vp = {round(Vp, 3)}, Vp angle = {round(math.degrees(VpAngle), 3)}")

    basicAccelerationHelper(a, b, c, theta2, theta3_2, theta4_2, omega2, omega3, omega4, alpha2, p, delta3)
    
    print("\nCROSSED")
    print(f"Theta3 = {round(math.degrees(theta3_1), 3)}, Theta4 = {round(math.degrees(theta4_1), 3)}")
    omega3 = a * omega2 * math.sin(theta4_1 - theta2) / (b * math.sin(theta3_1 - theta4_1))
    omega4 = a * omega2 * math.sin(theta2 - theta3_1) / (c * math.sin(theta4_1 - theta3_1))
    VaReal = a * omega2 * (-math.sin(theta2))
    VaImag = a * omega2 * math.cos(theta2)
    Va = math.sqrt(VaReal ** 2 + VaImag ** 2)
    VaAngle = positionSynthesis.invTan(VaReal, VaImag)
    VbReal = a * omega4 * (-math.sin(theta4_1))
    VbImag = a * omega4 * math.cos(theta4_1)
    Vb = math.sqrt(VbReal ** 2 + VbImag ** 2)
    VbAngle = positionSynthesis.invTan(VbReal, VbImag)
    print(f"Omega3 = {round(omega3, 3)}, Omega4 = {round(omega4, 3)}")
    print(f"Va = {round(Va, 3)}, Va angle = {round(math.degrees(VaAngle), 3)}")
    print(f"Vb = {round(Vb, 3)}, Vb angle = {round(math.degrees(VbAngle), 3)}")

    VpReal = VaReal + p * omega3 * (-math.sin(theta3_1 + delta3))
    VpImag = VaImag + p * omega3 * math.cos(theta3_1 + delta3)
    Vp = math.sqrt(VpReal ** 2 + VpImag ** 2)
    VpAngle = positionSynthesis.invTan(VpReal, VpImag)
    print(f"Vp = {round(Vp, 3)}, Vp angle = {round(math.degrees(VpAngle), 3)}")

    basicAccelerationHelper(a, b, c, theta2, theta3_1, theta4_1, omega2, omega3, omega4, alpha2, p, delta3)

def basicAccelerationHelper(a, b, c, theta2, theta3, theta4, omega2, omega3, omega4, alpha2, p, delta3):
    A = c * math.sin(theta4)
    B = b * math.sin(theta3)
    C = a * alpha2 * math.sin(theta2) + a * omega2 ** 2 * math.cos(theta2) + b * omega3 ** 2 * math.cos(theta3) - c * omega4 ** 2 * math.cos(theta4)
    D = c * math.cos(theta4)
    E = b * math.cos(theta3)
    F = a * alpha2 * math.cos(theta2) - a * omega2 ** 2 * math.sin(theta2) - b * omega3 ** 2 * math.sin(theta3) + c * omega4 ** 2 * math.sin(theta4)

    alpha3 = (C * D - A * F) / (A * E - B * D)
    alpha4 = (C * E - B * F) / (A * E - B * D)
    print(f"Alpha3 = {round(alpha3, 3)}, Alpha4 = {round(alpha4, 3)}")

    AaReal = -a * alpha2 * math.sin(theta2) - a * omega2 ** 2 * math.cos(theta2)
    AaImag = a * alpha2 * math.cos(theta2) - a * omega2 ** 2 * math.sin(theta2)
    Aa = math.sqrt(AaReal ** 2 + AaImag ** 2)
    AaAngle = positionSynthesis.invTan(AaReal, AaImag)
    AbReal = -c * alpha4 * math.sin(theta4) - c * omega4 ** 2 * math.cos(theta4)
    AbImag = c * alpha4 * math.cos(theta4) - c * omega4 ** 2 * math.sin(theta4)
    Ab = math.sqrt(AbReal ** 2 + AbImag ** 2)
    AbAngle = positionSynthesis.invTan(AbReal, AbImag)
    ApReal = AaReal + p * alpha3 * (-math.sin(theta3 + delta3)) - p * omega3 ** 2 * math.cos(theta3 + delta3)
    ApImag = AaImag + p * alpha3 * math.cos(theta3 + delta3) - p * omega3 ** 2 * math.sin(theta3 + delta3)
    Ap = math.sqrt(ApReal ** 2 + ApImag ** 2)
    ApAngle = positionSynthesis.invTan(ApReal, ApImag)
    print(f"Aa = {round(Aa, 3)}, Aa angle = {round(math.degrees(AaAngle), 3)}")
    print(f"Ab = {round(Ab, 3)}, Ab angle = {round(math.degrees(AbAngle), 3)}")
    print(f"Ap = {round(Ap, 3)}, Ap angle = {round(math.degrees(ApAngle), 3)}")

def crankSlider(a, b, c, theta2, omega2, alpha2):
    theta2 = math.radians(theta2)
    theta3_1 = (math.asin((a * math.sin(theta2) - c) / b))
    theta3_2 = math.asin(-(a * math.sin(theta2) - c) / b) + math.pi
    d1 = a * math.cos(theta2) - b * math.cos(theta3_1)
    d2 = a * math.cos(theta2) - b * math.cos(theta3_2)

    print("OPEN")
    crankSliderHelper(a, b, c, theta2, theta3_2, d2, omega2, alpha2)

    print("\nCROSSED")
    crankSliderHelper(a, b, c, theta2, theta3_1, d1, omega2, alpha2)

def crankSliderHelper(a, b, c, theta2, theta3, d, omega2, alpha2):
    print(f"Theta3 = {math.degrees(theta3)}, d = {d}")
    omega3 = a * math.cos(theta2) * omega2 / (b * math.cos(theta3))
    dDot = -a * omega2 * math.sin(theta2) + b * omega3 * math.sin(theta3)
    print(f"Omega3 = {omega3}, d* = {dDot}")
    VaReal = a * omega2 * (-math.sin(theta2))
    VaImag = a * omega2 * math.cos(theta2)
    Va = math.sqrt(VaReal ** 2 + VaImag ** 2)
    VaAngle = positionSynthesis.invTan(VaReal, VaImag)
    print(f"Va = {Va}, Va angle = {math.degrees(VaAngle)}")
    alpha3Num = a * alpha2 * math.cos(theta2) - a * omega2 ** 2 * math.sin(theta2) + b * omega3 ** 2 * math.sin(theta3)
    alpha3Den = b * math.cos(theta3)
    alpha3 = alpha3Num / alpha3Den
    dDoubleDot = -a * alpha2 * math.sin(theta2) - a * omega2 ** 2 * math.cos(theta2) + b * alpha3 * math.sin(theta3) + b * omega3 ** 2 * math.cos(theta3)
    AaReal = -a * alpha2 * math.sin(theta2) - a * omega2 ** 2 * math.cos(theta2)
    AaImag = a * alpha2 * math.cos(theta2) - a * omega2 ** 2 * math.sin(theta2)
    Aa = math.sqrt(AaReal ** 2 + AaImag ** 2)
    AaAngle = positionSynthesis.invTan(AaReal, AaImag)
    print(f"Aa = {Aa}, Aa angle = {math.degrees(AaAngle)}")
    print(f"Alpha3 = {alpha3}, d** = {dDoubleDot}")

def invCrankSlider(a, c, d, gamma, theta2, omega2, alpha2):
    theta2 = math.radians(theta2)
    gamma = math.radians(gamma)

    P = a * math.sin(theta2) * math.sin(gamma) + (a * math.cos(theta2) - d) * math.cos(gamma)
    Q = -a * math.sin(theta2) * math.cos(gamma) + (a * math.cos(theta2) - d) * math.sin(gamma)
    R = -c * math.sin(gamma)
    S = R - Q
    T = 2 * P
    U = Q + R

    theta4_1 = 2 * math.atan((-T + math.sqrt(T ** 2 - 4 * S * U)) / (2 * S))
    theta4_2 = 2 * math.atan((-T - math.sqrt(T ** 2 - 4 * S * U)) / (2 * S))
    theta3_1 = theta4_1 + gamma
    theta3_2 = theta4_2 + gamma
    b_1 = (a * math.sin(theta2) - c * math.sin(theta4_1)) / math.sin(theta3_1)
    b_2 = (a * math.sin(theta2) - c * math.sin(theta4_2)) / math.sin(theta3_2)


    print("OPEN")
    invCrankSliderHelper(a, c, d, gamma, theta2, omega2, alpha2, theta3_1, theta4_1, b_1)

    print("\nCROSSED")
    invCrankSliderHelper(a, c, d, gamma, theta2, omega2, alpha2, theta3_2, theta4_2, b_2)

def invCrankSliderHelper(a, c, d, gamma, theta2, omega2, alpha2, theta3, theta4, b):
    print(f"Theta3 = {math.degrees(theta3)}, Theta4 = {math.degrees(theta4)}, b = {b}")

    omega4 = a * omega2 * math.cos(theta2 - theta3) / (b + c * math.cos(theta4 - theta3))
    omega3 = omega4
    bDotNum = -a * omega2 * math.sin(theta2) + omega4 * (b * math.sin(theta3) + c * math.sin(theta4))
    bDot = bDotNum / math.cos(theta3)
    VaReal = a * omega2 * (-math.sin(theta2))
    VaImag = a * omega2 * math.cos(theta2)
    Va = math.sqrt(VaReal ** 2 + VaImag ** 2)
    VaAngle = positionSynthesis.invTan(VaReal, VaImag)
    VbReal = c * omega4 * (-math.sin(theta4))
    VbImag = c * omega4 * math.cos(theta4)
    Vb = math.sqrt(VbReal ** 2 + VbImag ** 2)
    VbAngle = positionSynthesis.invTan(VbReal, VbImag)
    print(f"Va = {Va}, Va Angle = {math.degrees(VaAngle)}, Vb = {Vb}, Vb Angle = {math.degrees(VbAngle)}")
    print(f"Omega3 = Omega4 = {omega4}, b* = {bDot}")

    alpha4Num1 = a * (alpha2 * math.cos(theta3 - theta2) + omega2 ** 2 * math.sin(theta3 - theta2))
    alpha4Num2 = c * omega4 ** 2 * math.sin(theta4 - theta3) - 2 * bDot * omega3
    alpha4Den = b + c * math.cos(theta3 - theta4)
    alpha4 = (alpha4Num1 + alpha4Num2) / alpha4Den
    AaReal = -a * alpha2 * math.sin(theta2) - a * omega2 ** 2 * math.cos(theta2)
    AaImag = a * alpha2 * math.cos(theta2) - a * omega2 ** 2 * math.sin(theta2)
    Aa = math.sqrt(AaReal ** 2 + AaImag ** 2)
    AaAngle = positionSynthesis.invTan(AaReal, AaImag)
    AbReal = -c * alpha4 * math.sin(theta4) - c * omega4 ** 2 * math.cos(theta4)
    AbImag = -c * alpha4 * (-math.cos(theta4)) - c * omega4 ** 2 * math.sin(theta4)
    Ab = math.sqrt(AbReal ** 2 + AbImag ** 2)
    AbAngle = positionSynthesis.invTan(AbReal, AbImag)
    bDoubleDotNum1 = a * omega2 ** 2 * (b * math.cos(theta3 - theta2) + c * math.cos(theta4 - theta2))
    bDoubleDotNum2 = a * alpha2 * (b * math.sin(theta2 - theta3) - c * math.sin(theta4 - theta2))
    bDoubleDotNum3 = 2 * bDot * c * omega4 * math.sin(theta4 - theta3)
    bDoubleDotNum4 = omega4 ** 2 * (b ** 2 + c ** 2 + 2 * b * c * math.cos(theta4 - theta3))
    bDoubleDotDen = b + c * math.cos(theta3 - theta4)
    bDoubleDot = -(bDoubleDotNum1 + bDoubleDotNum2 + bDoubleDotNum3 - bDoubleDotNum4) / bDoubleDotDen
    print(f"Aa = {Aa}, Aa Angle = {AaAngle}, Ab = {Ab}, Ab Angle = {math.degrees(AbAngle)}")
    print(f"alpha3 = alpha4 = {alpha4}, b** = {bDoubleDot}")
