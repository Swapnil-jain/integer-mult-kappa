"""O6/O7 exact checks on random integer data: outer wrapper and remainder pairing; negative controls."""
import random
rng = random.Random(7)
def F(H, D): return D[::-1], H[::-1]            # F_n(H,D) = (C D, C H)
ok = True
for n in (1, 2, 7, 64, 529 * 3 + 5):
    H = [rng.randint(-10**9, 10**9) for _ in range(n)]; D = [rng.randint(-10**9, 10**9) for _ in range(n)]
    h, d = H[:], D[:]
    d = [y - x for x, y in zip(h, d)]; h, d = F(h, d); d = [y + x for x, y in zip(h, d)]; h, d = F(h, d); d = [x - y for x, y in zip(h, d)]
    ok &= (h, d) == (D, H)
    # remainder pairing: n = m f + delta; first delta head atoms <-> last delta tail atoms reversed; the rest is F on two contiguous blocks
    m = 529; f, de = divmod(n, m)
    h2, d2 = H[:], D[:]
    for i in range(de): h2[i], d2[n - 1 - i] = D[n - 1 - i], H[i]
    a, b = F(H[de:], D[:n - de]); h2[de:] = a; d2[:n - de] = b
    ok &= (h2, d2) == F(H, D)
    # negative: wrong final step (D <- D - H) must fail
    h, d = H[:], D[:]
    d = [y - x for x, y in zip(h, d)]; h, d = F(h, d); d = [y + x for x, y in zip(h, d)]; h, d = F(h, d); d = [y - x for x, y in zip(h, d)]
    ok &= (h, d) != (D, H) or n == 0
print('O6/O7 wrapper and remainder pairing exact', 'PASS' if ok else 'FAIL')
