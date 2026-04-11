# d010_exact_formula2.sage
def compute_weights_QQ(points):
    n = len(points)
    weights = []
    for i in range(n):
        w = QQ(1)
        for j in range(n):
            if j != i:
                w *= QQ(-points[j]) / QQ(points[i] - points[j])
        weights.append(w)
    return weights

def norm1(weights):
    return sum(abs(w) for w in weights)

print("=== Test: is ||w||_1 a known q-series? ===\n")
for b in [3, 5, 7, 10, 20, 100]:
    q = QQ(1)/QQ(b)
    t = 50
    pts = [b^k for k in range(1, t+1)]
    ws = compute_weights_QQ(pts)
    actual = norm1(ws)
    
    # Candidate 1: b/(b-2)
    c1 = QQ(b)/QQ(b-2)
    # Candidate 2: 1/(q;q)_inf (just w_1)
    qpoch = prod(1 - q^k for k in range(1, 51))
    c2 = 1/qpoch
    # Candidate 3: b/(b-1) * 1/(q;q)_inf
    c3 = QQ(b)/QQ(b-1) * c2
    # Candidate 4: 1/(q;q)_inf * 1/(q;q)_inf  (square)
    c4 = c2 * c2
    # Candidate 5: 1/(q^2;q)_inf  (shift start)
    qpoch2 = prod(1 - q^(k+1) for k in range(1, 51))
    c5 = 1/qpoch2
    # Candidate 6: sum of q^(i*(i-1)/2) / (q;q)_i for i=0..inf (Ramanujan-type)
    c6 = sum(q^(i*(i-1)/2) / prod(1-q^k for k in range(1,i+1)) for i in range(51))
    # Candidate 7: product 1/(1-q^k)^2 — squared Euler function
    c7 = prod(1/(1-q^k)^2 for k in range(1, 51))
    # Candidate 8: sum_{i=1}^{inf} prod_{j=1,j!=i}^{inf} b^j/|b^i - b^j|
    # This is literally the norm as infinite sum. Check partial sums.
    
    print(f"b={b:4d} (q=1/{b})")
    print(f"  actual       = {float(actual):.12f}")
    print(f"  b/(b-2)      = {float(c1):.12f}  err={float(abs(actual-c1)):.2e}")
    print(f"  1/(q;q)_inf  = {float(c2):.12f}  err={float(abs(actual-c2)):.2e}  (=|w_1|)")
    print(f"  b/(b-1)/qpoc = {float(c3):.12f}  err={float(abs(actual-c3)):.2e}")
    print(f"  [1/(q;q)]^2  = {float(c4):.12f}  err={float(abs(actual-c4)):.2e}")
    print(f"  1/(q^2;q)_inf= {float(c5):.12f}  err={float(abs(actual-c5)):.2e}")
    print(f"  Ramanujan sum = {float(c6):.12f}  err={float(abs(actual-c6)):.2e}")
    print(f"  [1/(q;q)]^2  = {float(c7):.12f}  err={float(abs(actual-c7)):.2e}")
    print()
