import numpy as np

from green_laplace import reconstruct_unidirectional_timeseries


sample_count = 256
dt = 0.2
time = np.arange(sample_count) * dt
eta1 = 0.010 * np.cos(2 * np.pi * 5 * np.arange(sample_count) / sample_count)
eta1 += 0.006 * np.cos(
    2 * np.pi * 8 * np.arange(sample_count) / sample_count + 0.3
)

result = reconstruct_unidirectional_timeseries(
    eta1,
    time,
    depth=1.0,
    energy_fraction=0.999,
    order=3,
)
print("retained energy:", result.audit["retained_positive_frequency_energy"])
print("eta20 L2 norm:", np.linalg.norm(result.components["eta20"]))
print("eta22 L2 norm:", np.linalg.norm(result.components["eta22"]))
print("eta33 L2 norm:", np.linalg.norm(result.components["eta33"]))
