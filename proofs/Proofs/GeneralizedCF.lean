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

end Proofs
