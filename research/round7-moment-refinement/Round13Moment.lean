import Round7

/- Finite rational certificate only. Analytic interpretation uses the atanh
   remainder and exponential Taylor-tail inequalities, which are not formalized
   here. The physical network is not formalized here. Imported Round7.bitHist
   is the complete frozen histogram, not an arbitrary supplied sum.
   All exp terms are rounded upward; this avoids enormous power denominators. -/

set_option maxRecDepth 100000

namespace Round13Moment

def scale : Nat := 10^40

def ceil40 (a : Q) : Q :=
  norm ((a.1 * scale + a.2 - 1) / a.2, scale)

-- Pair (term, sum) scaled by 10^40, with each term rounded upward.
-- Inductively term_j/scale >= a^j/j! for a >= 0.
def expUpLoop (a : Q) : Nat → Nat × Nat
  | 0 => (scale, scale)
  | j + 1 =>
    let previous := expUpLoop a j
    let den := (j+1)*a.2
    let term := (previous.1*a.1 + den - 1)/den
    (term, previous.2 + term)

-- Tail after degree8 <= (upper_term8*a/9)/(1-a/10), for a < 1/1000.
def exp8Up (a : Q) : Q :=
  let p := expUpLoop a 8
  ceil40 (qadd (p.2, scale)
    (qdiv (qdiv (qmul (p.1, scale) a) (9, 1))
      (qsub (1, 1) (qdiv a (10, 1)))))

def momentTaylorOK (h : List (Nat × Nat)) (a : Q) (m W : Nat) : Bool :=
  let step := fun (acc : Bool × Q) (p : Nat × Nat) =>
    let v := ceil40 (qmul a (lnUp (m, p.1)))
    let ok := acc.1 && (0 < p.1) && (p.1 < m) && qlt v (1, 1000)
    if ok then
      (true, qadd acc.2 (qmul (p.1*p.2, m*W) (exp8Up v)))
    else (false, acc.2)
  let result := h.foldl step (true, (0, 1))
  result.1 && qlt result.2 (1, 1)

theorem refined_bit_moment :
    momentTaylorOK bitHist (15995753310011, 250000000000000000)
      529 108516254 = true := by decide

theorem refined_kappa_assembly :
    kappaOK (15995753310011, 250000000000000000)
      (36926111, 500000000000) (1, 1000)
      (9999360210803241, 10000000000000000)
      (14695, 1) (14694, 1) (799736495947, 12500000000000000)
      576 119453132304 = true := by decide +kernel

end Round13Moment

#print axioms Round13Moment.refined_bit_moment
#print axioms Round13Moment.refined_kappa_assembly
