

def compute_macroscopic(self, e):
    """ rho, ux, uy """
    temp_x = self.f * e.x[:, None, None]
    temp_y = self.f * e.y[:, None, None]
    self.rho = self.f.sum(axis=0)
    self.u.x = (self.rho / (1e-30 + self.rho) ** 2) * temp_x.sum(axis=0)
    self.u.y = (self.rho / (1e-30 + self.rho) ** 2) * temp_y.sum(axis=0)
    # account for obstacles
    self.rho *= 1. * (self.solid <= 0)
    self.u.x *= 1. * (self.solid <= 0)
    self.u.y *= 1. * (self.solid <= 0)
    self.check_cfl()