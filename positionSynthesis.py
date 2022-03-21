import math
import numpy
from scipy.optimize import fsolve


def threePositionFixedPivots(p21x, p21y, p31x, p31y, alpha2, alpha3, O2x, O2y, O4x, O4y):
    alpha2 = math.radians(alpha2)
    alpha3 = math.radians(alpha3)

    r2x = O2x + p21x
    r2y = O2y + p21y
    r3x = O2x + p31x
    r3y = O2y + p31y

    betaArray = threePositionFixedPivotsHelper(O2x, O2y, r2x, r2y, r3x, r3y, alpha2, alpha3)
    print(f'beta2 = {round(math.degrees(betaArray[0]), 3)} ({round(math.degrees(betaArray[0]) - 360, 3)}), beta3 = {round(math.degrees(betaArray[1]), 3)} ({round(math.degrees(betaArray[1]) - 360, 3)})')
    
    r2x = O4x + p21x
    r2y = O4y + p21y
    r3x = O4x + p31x
    r3y = O4y + p31y
    
    gammaArray = threePositionFixedPivotsHelper(O4x, O4y, r2x, r2y, r3x, r3y, alpha2, alpha3)
    print(
        f'gamma2 = {round(math.degrees(gammaArray[0]), 3)} ({round(math.degrees(gammaArray[0]) - 360, 3)}), gamma3 = {round(math.degrees(gammaArray[1]), 3)} ({round(math.degrees(gammaArray[1]) - 360, 3)})')

    print()
    p21 = math.sqrt(p21x ** 2 + p21y ** 2)
    p31 = math.sqrt(p31x ** 2 + p31y ** 2)
    delta2 = invTan(p21x, p21y)
    delta3 = invTan(p31x, p31y)
    threePosition(p21, p31, math.degrees(delta2), math.degrees(delta3), math.degrees(alpha2), math.degrees(alpha3), math.degrees(betaArray[0]), math.degrees(betaArray[1]), math.degrees(gammaArray[0]), math.degrees(gammaArray[1]))


def threePositionFixedPivotsHelper(Ox, Oy, r2x, r2y, r3x, r3y, alpha2, alpha3):
    zeta1 = invTan(Ox, Oy)
    zeta2 = invTan(r2x, r2y)
    zeta3 = invTan(r3x, r3y)
    r1 = math.sqrt(Ox ** 2 + Oy ** 2)
    r2 = math.sqrt(r2x ** 2 + r2y ** 2)
    r3 = math.sqrt(r3x ** 2 + r3y ** 2)

    C1 = r3 * math.cos(alpha2 + zeta3) - r2 * math.cos(alpha3 + zeta2)
    C2 = r3 * math.sin(alpha2 + zeta3) - r2 * math.sin(alpha3 + zeta2)
    C3 = r1 * math.cos(alpha3 + zeta1) - r3 * math.cos(zeta3)
    C4 = -r1 * math.sin(alpha3 + zeta1) + r3 * math.sin(zeta3)
    C5 = r1 * math.cos(alpha2 + zeta1) - r2 * math.cos(zeta2)
    C6 = -r1 * math.sin(alpha2 + zeta1) + r2 * math.sin(zeta2)

    A1 = -C3 ** 2 - C4 ** 2
    A2 = C3 * C6 - C4 * C5
    A3 = -C4 * C6 - C3 * C5
    A4 = C2 * C3 + C1 * C4
    A5 = C4 * C5 - C3 * C6
    A6 = C1 * C3 - C2 * C4

    K1 = A2 * A4 + A3 * A6
    K2 = A3 * A4 + A5 * A6
    K3 = (A1 ** 2 - A2 ** 2 - A3 ** 2 - A4 ** 2 - A6 ** 2) / 2

    numer = K2 + math.sqrt(K1 ** 2 + K2 ** 2 - K3 ** 2)
    denom = K1 + K3
    beta3 = 2 * math.atan(numer / denom)
    if(round(beta3, 2) == round(alpha3, 2)):
        numer = K2 - math.sqrt(K1 ** 2 + K2 ** 2 - K3 ** 2)
        beta3 = 2 * math.atan(numer / denom)
    # beta3 = 2 * invTan(K1 + K3, K2 ** 2 + math.sqrt(K1 ** 2 + K2 ** 2 - K3 ** 2))
    beta2 = invTan(-(A5 * math.sin(beta3) + A3 * math.cos(beta3) + A6),
                   -(A3 * math.sin(beta3) + A2 * math.cos(beta3) + A4))

    return [beta2, beta3]

def threePosition(p21, p31, delta2, delta3, alpha2, alpha3, beta2, beta3, gamma2, gamma3):
    delta2 = math.radians(delta2)
    delta3 = math.radians(delta3)
    alpha2 = math.radians(alpha2)
    alpha3 = math.radians(alpha3)
    beta2 = math.radians(beta2)
    beta3 = math.radians(beta3)
    gamma2 = math.radians(gamma2)
    gamma3 = math.radians(gamma3)

    M3 = threePositionHelper(p21, p31, delta2, delta3, alpha2, alpha3, beta2, beta3)
    print(f'w1x = {round(M3[0], 3)}, w1y = {round(M3[1], 3)}')
    print(f'z1x = {round(M3[2], 3)}, z1y = {round(M3[3], 3)}')
    W1 = math.sqrt(M3[0] ** 2 + M3[1] ** 2)
    theta = invTan(M3[0], M3[1])
    Z1 = math.sqrt(M3[2] ** 2 + M3[3] ** 2)
    phi = invTan(M3[2], M3[3])
    print(f'W1 = {round(W1, 3)}, theta = {round(math.degrees(theta), 3)}')
    print(f'Z1 = {round(Z1, 3)}, phi = {round(math.degrees(phi), 3)}\n')

    M3 = threePositionHelper(p21, p31, delta2, delta3, alpha2, alpha3, gamma2, gamma3)
    print(f'u1x = {round(M3[0], 3)}, u1y = {round(M3[1], 3)}')
    print(f's1x = {round(M3[2], 3)}, s1y = {round(M3[3], 3)}')
    U1 = math.sqrt(M3[0] ** 2 + M3[1] ** 2)
    sigma = invTan(M3[0], M3[1])
    S1 = math.sqrt(M3[2] ** 2 + M3[3] ** 2)
    psi = invTan(M3[2], M3[3])
    print(f'U1 = {round(U1, 3)}, sigma = {round(math.degrees(sigma), 0)}')
    print(f'S1 = {round(S1, 3)}, psi = {round(math.degrees(psi), 3)}\n')

    v1x = Z1 * math.cos(phi) - S1 * math.cos(psi)
    v1y = Z1 * math.sin(phi) - S1 * math.sin(psi)
    thetaV = invTan(v1x, v1y)
    v = math.sqrt(v1x ** 2 + v1y ** 2)
    print(f'V = {round(v, 3)} @ {round(math.degrees(thetaV), 3)} ({round(math.degrees(thetaV) - 360, 2)})')

    g1x = W1 * math.cos(theta) + v * math.cos(thetaV) - U1 * math.cos(sigma)
    g1y = W1 * math.sin(theta) + v * math.sin(thetaV) - U1 * math.sin(sigma)
    thetaG = invTan(g1x, g1y)
    g = math.sqrt(g1x ** 2 + g1y ** 2)
    print(f'G = {round(g, 3)} @ {round(math.degrees(thetaG), 3)} ({round(math.degrees(thetaG) - 360, 3)})')


def threePositionHelper(p21, p31, delta2, delta3, alpha2, alpha3, beta2, beta3):
    A = math.cos(beta2) - 1
    B = math.sin(beta2)
    C = math.cos(alpha2) - 1
    D = math.sin(alpha2)
    E = p21 * math.cos(delta2)
    F = math.cos(beta3) - 1
    G = math.sin(beta3)
    H = math.cos(alpha3) - 1
    K = math.sin(alpha3)
    L = p31 * math.cos(delta3)
    M = p21 * math.sin(delta2)
    N = p31 * math.sin(delta3)

    M1 = [[A, -B, C, -D], [F, -G, H, -K], [B, A, D, C], [G, F, K, H]]
    M2 = [E, L, M, N]

    M1Inv = numpy.linalg.inv(M1)
    return(numpy.dot(M1Inv, M2))

def twoPositionsM1():
    return

def twoPositionM2(z, phi, beta2, s, psi, gamma2, p2x, p2y, ang1, ang2):
    phi = math.radians(phi)
    beta2 = math.radians(beta2)
    psi = math.radians(psi)
    gamma2 = math.radians(gamma2)
    ang1 = math.radians(ang1)
    ang2 = math.radians(ang2)

    p = math.sqrt(p2y ** 2 + p2x ** 2)
    delta2 = invTan(p2x, p2y)
    alpha2 = ang2 - ang1
    print(f'p = {p}, delta2 = {math.degrees(delta2)}, alpha2 = {math.degrees(alpha2)}')

    A = math.cos(beta2) - 1
    B = math.sin(beta2)
    C = math.cos(alpha2) - 1
    D = math.sin(alpha2)
    E = p * math.cos(delta2)
    F = p * math.sin(delta2)
    z1x = z * math.cos(phi)
    z1y = z * math.sin(phi)

    w1x = (A * (-C * z1x + D * z1y + E) + B * (-C * z1y - D * z1x + F)) / (-2 * A)
    w1y = (A * (-C * z1y - D * z1x + F) + B * (C * z1x - D * z1y - E)) / (-2 * A)
    print(f'w1x = {w1x}, w1y = {w1y}')
    w = math.sqrt(w1x ** 2 + w1y ** 2)
    theta = invTan(w1x, w1y)
    print(f'w = {w}, theta = {math.degrees(theta)} ({math.degrees(theta) - 360})')

    A = math.cos(gamma2) - 1
    B = math.sin(gamma2)
    C = math.cos(alpha2) - 1
    D = math.sin(alpha2)
    E = p * math.cos(delta2)
    F = p * math.sin(delta2)
    s1x = s * math.cos(psi)
    s1y = s * math.sin(psi)

    u1x = (A * (-C * s1x + D * s1y + E) + B * (-C * s1y - D * s1x + F)) / (-2 * A)
    u1y = (A * (-C * s1y - D * s1x + F) + B * (C * s1x - D * s1y - E)) / (-2 * A)
    print(f'u1x = {u1x}, u1y = {u1y}')
    u = math.sqrt(u1x ** 2 + u1y ** 2)
    sigma = invTan(u1x, u1y)
    print(f'u = {u}, sigma = {math.degrees(sigma)} ({math.degrees(sigma) - 360})')

    v1x = z * math.cos(phi) - s * math.cos(psi)
    v1y = z * math.sin(phi) - s * math.sin(psi)
    thetaV = invTan(v1x, v1y)
    v = math.sqrt(v1x ** 2 + v1y ** 2)
    print(f'V = {v} @ {math.degrees(thetaV)} ({math.degrees(thetaV) - 360})')
    
    g1x = w * math.cos(theta) + v * math.cos(thetaV) - u * math.cos(sigma)
    g1y = w * math.sin(theta) + v * math.sin(thetaV) - u * math.sin(sigma)
    thetaG = invTan(g1x, g1y)
    g = math.sqrt(g1x ** 2 + g1y ** 2)
    print(f'G = {g} @ {math.degrees(thetaG)} ({math.degrees(thetaG) - 360})')

    origin2x = -z * math.cos(phi) - w * math.cos(theta)
    origin2y = -z * math.sin(phi) - w * math.sin(theta)
    origin4x = -s * math.cos(psi) - u * math.cos(sigma)
    origin4y = -s * math.sin(psi) - u * math.sin(sigma)
    print(f'origin 2: {origin2x}, {origin2y}')
    print(f'origin 4: {origin4x}, {origin4y}')

    return

def invTan(x, y):
    preAngle = math.atan(y / x)
    if (x < 0) and (y < 0):
        angle = preAngle + math.pi
    elif x < 0:
        angle = math.pi + preAngle
    elif y < 0:
        angle = 2 * math.pi + preAngle
    else:
        angle = preAngle
    return angle

def myFunction(z, A, B, C, D, E, F):
    x = z[0]
    y = z[1]

    func = numpy.empty((2))
    func[0] = 0

    return
