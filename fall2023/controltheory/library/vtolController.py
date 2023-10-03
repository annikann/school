import library.vtolParam as P

class vtolController:
    def __init__(self):
        self.mr = P.mr
        self.mc = P.mc
        self.kP = 0.09
        self.kD = 0.75

    def update(self, hc, state):
        h = state[2]
        hdot = state[3]
        feq = P.g*(self.mc + 2.*self.mr)
        fc = self.kP*(hc - h)  - self.kD*hdot
        f = feq + fc
        return f[0]