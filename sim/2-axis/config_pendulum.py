"""
    Pendulum
    Mr =  mass of rod (kg)
    Yr =  normalized distance of rod along rod
    Mw =  mass of wheel (kg)
    Yw =  normalized distance of wheel along rod
    lp =  length of pendulum (m)
    Mp = Mw*Yw + Mr*Yr

    X-axis Constants
    Ix =  pendulum MOI (kg*m^2)
    Bx =  friction in pendulum bearing (Nms)

    Y-axis Constants
    Iy =  pendulum MOI (kg*m^2)
    By =  friction in pendulum bearing (Nms)

    Z-axis Constants
    Iz =  pendulum MOI (kg*m^2)
    Bz =  friction in pendulum bearing (Nms)
    Fspin_z =  force spinning around rod axis, term needed for stability
"""

estimate = {'Mr': 0.2,
            'Yr': 0.5,
            'Mw': 0.8,
            'Yw': 1,
            'lp': 0.33,
            'Ix': 0.05,
            'Bx': 0,
            'Iy': 0.05,
            'By': 0,
            'Iz': 0.05,
            'Bz': 0,
            'Fspin_z': 1}
