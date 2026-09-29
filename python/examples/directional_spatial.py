import numpy as np

from green_laplace import reconstruct_directional_field


length = 2 * np.pi
points = 64
x = np.arange(points) * length / points
y = np.arange(points) * length / points
X, Y = np.meshgrid(x, y)
eta1 = 0.010 * np.cos(2 * X) + 0.006 * np.cos(3 * X + Y + 0.3)

result = reconstruct_directional_field(
    eta1,
    lx=length,
    ly=length,
    depth=1.0,
    order=3,
    eta22_rank=6,
)
print("eta22 L2 norm:", np.linalg.norm(result.components["eta22"]))
print("psi22 L2 norm:", np.linalg.norm(result.components["psi22"]))
print("eta33 L2 norm:", np.linalg.norm(result.components["eta33"]))
print("psi33 L2 norm:", np.linalg.norm(result.components["psi33"]))
