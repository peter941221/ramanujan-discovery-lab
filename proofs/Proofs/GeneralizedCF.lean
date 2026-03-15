import Mathlib

namespace Proofs

structure GCFData (R : Type*) where
  b0 : R
  a : Nat -> R
  b : Nat -> R

section Semiring

variable {R : Type*} [Semiring R]

def continuantNum (g : GCFData R) : Nat -> R
  | 0 => g.b0
  | 1 => g.b 0 * g.b0 + g.a 0
  | n + 2 => g.b (n + 1) * continuantNum g (n + 1) + g.a (n + 1) * continuantNum g n

def continuantDen (g : GCFData R) : Nat -> R
  | 0 => 1
  | 1 => g.b 0
  | n + 2 => g.b (n + 1) * continuantDen g (n + 1) + g.a (n + 1) * continuantDen g n

@[simp]
theorem continuantNum_zero (g : GCFData R) :
    continuantNum g 0 = g.b0 := rfl

@[simp]
theorem continuantNum_one (g : GCFData R) :
    continuantNum g 1 = g.b 0 * g.b0 + g.a 0 := rfl

@[simp]
theorem continuantNum_succ_succ (g : GCFData R) (n : Nat) :
    continuantNum g (n + 2) =
      g.b (n + 1) * continuantNum g (n + 1) + g.a (n + 1) * continuantNum g n := rfl

@[simp]
theorem continuantDen_zero (g : GCFData R) :
    continuantDen g 0 = 1 := rfl

@[simp]
theorem continuantDen_one (g : GCFData R) :
    continuantDen g 1 = g.b 0 := rfl

@[simp]
theorem continuantDen_succ_succ (g : GCFData R) (n : Nat) :
    continuantDen g (n + 2) =
      g.b (n + 1) * continuantDen g (n + 1) + g.a (n + 1) * continuantDen g n := rfl

end Semiring

section Equivalence

variable {R : Type*} [CommSemiring R]

/-!
Equivalence transforms for generalized continued fractions.

Given stage scales `rₙ` (here indexed as `r k := r_{k+1}`), define a transformed continued fraction by
- `b₀` unchanged
- `b_{k+1}' = r_{k+1} * b_{k+1}` (Lean index: `b' k = r k * b k`)
- `a₁' = r₁ * a₁` (Lean index: `a' 0 = r 0 * a 0`)
- `a_{k+1}' = r_k * r_{k+1} * a_{k+1}` for `k ≥ 1` (Lean index: `a' (k+1) = r k * r (k+1) * a (k+1)`).

This is the coefficient-level rule used in the Python pipeline and preserves all finite convergents.
-/

def scaleProd (r : Nat → R) : Nat → R
  | 0 => 1
  | n + 1 => scaleProd r n * r n

@[simp] theorem scaleProd_zero (r : Nat → R) : scaleProd r 0 = 1 := rfl
@[simp] theorem scaleProd_succ (r : Nat → R) (n : Nat) : scaleProd r (n + 1) = scaleProd r n * r n := rfl

def equivalenceTransform (g : GCFData R) (r : Nat → R) : GCFData R where
  b0 := g.b0
  b := fun k => r k * g.b k
  a := fun k =>
    match k with
    | 0 => r 0 * g.a 0
    | k + 1 => r k * r (k + 1) * g.a (k + 1)

@[simp] theorem equivalenceTransform_b (g : GCFData R) (r : Nat → R) (k : Nat) :
    (equivalenceTransform g r).b k = r k * g.b k := rfl

@[simp] theorem equivalenceTransform_a_zero (g : GCFData R) (r : Nat → R) :
    (equivalenceTransform g r).a 0 = r 0 * g.a 0 := rfl

@[simp] theorem equivalenceTransform_a_succ (g : GCFData R) (r : Nat → R) (k : Nat) :
    (equivalenceTransform g r).a (k + 1) = r k * r (k + 1) * g.a (k + 1) := rfl

theorem continuantNum_equivalenceTransform (g : GCFData R) (r : Nat → R) :
    ∀ n : Nat, continuantNum (equivalenceTransform g r) n = scaleProd r n * continuantNum g n := by
  intro n
  induction' n using Nat.twoStepInduction with n ih0 ih1
  · simp [continuantNum, equivalenceTransform, scaleProd]
  · simp [continuantNum, equivalenceTransform, scaleProd, mul_add, mul_assoc, mul_left_comm, mul_comm]
  · -- step
    -- Unfold the recurrence at stage `n+2`, rewrite subcalls via IH, then normalize.
    simp [continuantNum, scaleProd, ih0, ih1, mul_add, mul_assoc, mul_left_comm, mul_comm]

theorem continuantDen_equivalenceTransform (g : GCFData R) (r : Nat → R) :
    ∀ n : Nat, continuantDen (equivalenceTransform g r) n = scaleProd r n * continuantDen g n := by
  intro n
  induction' n using Nat.twoStepInduction with n ih0 ih1
  · simp [continuantDen, scaleProd]
  · simp [continuantDen, scaleProd]
  ·
    simp [continuantDen, scaleProd, ih0, ih1, mul_add, mul_assoc, mul_left_comm, mul_comm]

end Equivalence

section DivisionRing

variable {R : Type*} [DivisionRing R]

def convergent (g : GCFData R) (n : Nat) : R :=
  continuantNum g n / continuantDen g n

@[simp]
theorem convergent_zero (g : GCFData R) :
    convergent g 0 = g.b0 := by
  simp [convergent]

@[simp]
theorem convergent_one (g : GCFData R) :
    convergent g 1 = (g.b 0 * g.b0 + g.a 0) / g.b 0 := by
  simp [convergent]

end DivisionRing

section EquivalenceConvergent

variable {R : Type*} [Field R]

theorem scaleProd_ne_zero (r : Nat → R) (hr : ∀ k : Nat, r k ≠ 0) :
    ∀ n : Nat, scaleProd r n ≠ 0 := by
  intro n
  induction n with
  | zero =>
      simp [scaleProd]
  | succ n ih =>
      exact mul_ne_zero ih (hr n)

theorem convergent_equivalenceTransform (g : GCFData R) (r : Nat → R) (hr : ∀ k : Nat, r k ≠ 0) :
    ∀ n : Nat, convergent (equivalenceTransform g r) n = convergent g n := by
  intro n
  have hs : scaleProd r n ≠ 0 := scaleProd_ne_zero (r := r) hr n
  -- Use the continuant scaling lemmas and cancel the common factor.
  simp [convergent, continuantNum_equivalenceTransform, continuantDen_equivalenceTransform, hs, mul_div_mul_left]

end EquivalenceConvergent

end Proofs
