"""
Plotting from Database
Shawn Macon
09.09.2025
"""


import numpy as np
import matplotlib.pyplot as plt

# =============================
# 1. Create a simple phantom
# =============================

N = 64   # grid size
x = np.linspace(-1, 1, N)
y = np.linspace(-1, 1, N)
xx, yy = np.meshgrid(x, y)

phantom = ((xx**2 + yy**2) < 0.3**2).astype(float)  # a simple disk


# ==========================================
# 2. Build Expanded-Ray Projection Matrix A
# ==========================================

def ray_pixel_intersection_length(x0, y0, theta, pixel_center, pixel_size=2/N):
    """
    Approximate ray–pixel intersection length using a footprint model.
    Very fast surrogate: weight = max(0, 1 - distance_to_ray / pixel_halfdiag)
    """
    px, py = pixel_center
    # Ray direction
    dx = np.cos(theta)
    dy = np.sin(theta)

    # Distance from pixel center to infinite line of ray
    dist = abs(dy*(px - x0) - dx*(py - y0))

    pixel_radius = np.sqrt(2)*(pixel_size/2)

    if dist > pixel_radius:
        return 0.0
    else:
        # Smooth footprint (expanded ray)
        return 1 - dist/pixel_radius


# Define projection geometry
n_angles = 40
n_rays = N
angles = np.linspace(0, np.pi, n_angles)

# Total number of equations
M = n_angles * n_rays

# Allocate projection matrix
A = np.zeros((M, N*N), dtype=float)
b = np.zeros(M, dtype=float)

pixel_centers = np.stack([xx.flatten(), yy.flatten()], axis=1)

row = 0
for th in angles:
    # Ray positions are shifts along perpendicular direction
    perp = np.array([-np.sin(th), np.cos(th)])
    shifts = np.linspace(-1, 1, n_rays)

    for s in shifts:
        x0, y0 = s * perp      # ray origin on detector line

        # Compute all weights for this ray
        weights = np.zeros(N*N)
        for i, (px, py) in enumerate(pixel_centers):
            weights[i] = ray_pixel_intersection_length(
                x0, y0, th, (px, py)
            )

        A[row, :] = weights
        b[row] = np.dot(weights, phantom.flatten())
        row += 1


# =============================
# 3. Run ART iterations
# =============================

def ART_reconstruction(A, b, n_iter=5, relaxation=1.0):
    M, K = A.shape
    x = np.zeros(K)
    residuals = []

    for it in range(n_iter):
        for i in range(M):
            ai = A[i]
            denom = ai @ ai
            if denom == 0:
                continue

            # Kaczmarz/ART update
            x = x + relaxation * (b[i] - ai @ x) / denom * ai

        # track residual
        residuals.append(np.linalg.norm(A @ x - b))

    return x, residuals


recon_vec, residuals = ART_reconstruction(A, b, n_iter=15)
recon = recon_vec.reshape(N, N)


# =============================
# 4. Plot results
# =============================

plt.figure(figsize=(14, 6))

plt.subplot(1, 3, 1)
plt.title("Original Phantom")
plt.imshow(phantom, cmap='gray')
plt.axis('off')

plt.subplot(1, 3, 2)
plt.title("ART Reconstruction")
plt.imshow(recon, cmap='gray')
plt.axis('off')

plt.subplot(1, 3, 3)
plt.title("Residual Norm vs Iteration")
plt.plot(residuals, marker='o')
plt.xlabel("Iteration")
plt.ylabel("||Ax - b||")

plt.tight_layout()
plt.show()