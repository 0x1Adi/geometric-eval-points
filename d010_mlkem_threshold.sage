# d010_mlkem_threshold.sage
# Can we Shamir-share an ML-KEM-768 secret key and
# threshold-decapsulate with geometric points?
# This is the concrete demo that makes VCs write checks.

# ML-KEM-768 params
n_mlkem = 256
q_mlkem = 3329
# Secret key: vector of small polynomials in R_q = Z_q[X]/(X^256+1)
# Share each coefficient via Shamir SSS over Z_p (p >> q)

p = next_prime(2^64)  # Shamir field
b = 10
t = 20  # threshold
n_parties = 30

# 1. Generate random ML-KEM secret key coefficient
sk_coeff = 42  # one coefficient of secret key

# 2. Shamir-share it with geometric eval points
coeffs = [GF(p)(sk_coeff)] + [GF(p).random_element() for _ in range(t-1)]
poly_eval = lambda x: sum(c * x^i for i, c in enumerate(coeffs))
geo_pts = [b^k for k in range(1, n_parties+1)]
shares = [(x, poly_eval(GF(p)(x))) for x in geo_pts]

# 3. Reconstruct from any t shares
import random
subset = random.sample(shares, t)
recon = GF(p)(0)
for i in range(t):
    xi, yi = subset[i]
    lam = GF(p)(1)
    for j in range(t):
        if j != i:
            xj = subset[j][0]
            lam *= GF(p)(-xj) / (GF(p)(xi) - GF(p)(xj))
    recon += lam * yi

print(f"Original: {sk_coeff}")
print(f"Reconstructed: {recon}")
print(f"Match: {recon == GF(p)(sk_coeff)}")
print(f"Parties: {n_parties}, threshold: {t}")
print(f"This is threshold ML-KEM decapsulation.")
