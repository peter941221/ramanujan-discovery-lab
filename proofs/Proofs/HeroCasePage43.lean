import Mathlib
import Mathlib.Algebra.Polynomial.Laurent

namespace Proofs
namespace HeroCase

set_option linter.unnecessarySimpa false

open scoped LaurentPolynomial
open LaurentPolynomial

noncomputable section

/-!
Bounded family exclusions for the hero case `cb60fd71d1d7`.

This file formalizes the "page-43 monomial substitution" no-hit results used in the
research notes.  We work in Laurent polynomials `Rat[T;T⁻¹]` so integer shifts are
represented directly.

The target is the step-reduced reciprocal coefficients

* `aₙ = T^n + T^(2n)`
* `bₙ = 1 + T^n`.

The page-43 families are represented in the same coefficient language, and we prove that
they cannot match the target even at stages `n = 1,2`.  This is stronger than the bounded
search box used in the Python pipeline (`A,B,L ∈ [-3,3]`, `stages=3`), hence it implies
zero hits in that box.
-/

abbrev LPoly := Rat[T;T⁻¹]

def targetA (n : Nat) : LPoly :=
  T (Int.ofNat n) + T (Int.ofNat (2 * n))

def targetB (n : Nat) : LPoly :=
  C (1 : Rat) + T (Int.ofNat n)

namespace Page43

/-!
Family definitions (matching `research.py:page43_monomial_parameter_search`).

We only need the low-stage shapes to state exclusion theorems.
-/

def f4A (A L : ℤ) (alpha gamma : Rat) (n : Nat) : LPoly :=
  C alpha * T (A + 1) + C gamma * T (L + Int.ofNat n)

def f4B (A B : ℤ) (alpha beta : Rat) (n : Nat) : LPoly :=
  C (1 : Rat) + C (-alpha) * T (A + 1) + C beta * T (B + Int.ofNat n)

def f2A (A B L : ℤ) (alpha beta gamma : Rat) (n : Nat) : LPoly :=
  C gamma * T (L + Int.ofNat n) + C (-(alpha * beta)) * T (A + B + Int.ofNat (2 * n))

def f2B (A B : ℤ) (alpha beta : Rat) (n : Nat) : LPoly :=
  C (1 : Rat) + C beta * T (B + Int.ofNat n) + C alpha * T (A + Int.ofNat n + 1)

private lemma C_mul_T_eq_single (c : Rat) (e : ℤ) :
    (C c * T e : LPoly) = Finsupp.single e c := by
  simp [LaurentPolynomial.single_eq_C_mul_T]

private lemma C_mul_C_apply_zero (a b : Rat) : ((C a * C b : LPoly) 0) = a * b := by
  have h := congrArg (fun p : LPoly => p 0) ((C : Rat →+* LPoly).map_mul a b)
  -- `h : (C (a*b)) 0 = (C a * C b) 0`
  have h0 : ((C a * C b : LPoly) 0) = (C (a * b) : LPoly) 0 := h.symm
  calc
    ((C a * C b : LPoly) 0) = (C (a * b) : LPoly) 0 := h0
    _ = a * b := by
      simpa using (LaurentPolynomial.C_apply (R := Rat) (t := (a * b)) (n := (0 : ℤ)))

private lemma targetA_support_one : (targetA 1).support = ({(1 : ℤ), (2 : ℤ)} : Finset ℤ) := by
  classical
  -- `T k` is `single k 1`.
  simpa [targetA, LaurentPolynomial.T, -LaurentPolynomial.single_eq_C_mul_T] using
    (Finsupp.support_single_add_single (f₁ := (1 : ℤ)) (f₂ := (2 : ℤ))
      (g₁ := (1 : Rat)) (g₂ := (1 : Rat)) (by decide) (by decide) (by decide))

private lemma targetA_support_two : (targetA 2).support = ({(2 : ℤ), (4 : ℤ)} : Finset ℤ) := by
  classical
  simpa [targetA, LaurentPolynomial.T, -LaurentPolynomial.single_eq_C_mul_T] using
    (Finsupp.support_single_add_single (f₁ := (2 : ℤ)) (f₂ := (4 : ℤ))
      (g₁ := (1 : Rat)) (g₂ := (1 : Rat)) (by decide) (by decide) (by decide))

private lemma targetB_support_one : (targetB 1).support = ({(0 : ℤ), (1 : ℤ)} : Finset ℤ) := by
  classical
  -- `C 1` is `single 0 1`.
  simpa [targetB, LaurentPolynomial.C, LaurentPolynomial.T, -LaurentPolynomial.single_eq_C_mul_T] using
    (Finsupp.support_single_add_single (f₁ := (0 : ℤ)) (f₂ := (1 : ℤ))
      (g₁ := (1 : Rat)) (g₂ := (1 : Rat)) (by decide) (by decide) (by decide))

private lemma support_f4A_subset (A L : ℤ) (alpha gamma : Rat) (n : Nat) :
    (f4A A L alpha gamma n).support ⊆ ({A + 1, L + Int.ofNat n} : Finset ℤ) := by
  classical
  have h := (Finsupp.support_add (g₁ := (C alpha * T (A + 1) : LPoly))
    (g₂ := (C gamma * T (L + Int.ofNat n) : LPoly)))
  refine h.trans ?_
  refine Finset.union_subset_iff.2 ?_
  constructor
  · exact (LaurentPolynomial.support_C_mul_T (R := Rat) (a := alpha) (n := A + 1)).trans (by simp)
  · exact (LaurentPolynomial.support_C_mul_T (R := Rat) (a := gamma) (n := L + Int.ofNat n)).trans (by simp)

private lemma support_f2A_subset (A B L : ℤ) (alpha beta gamma : Rat) (n : Nat) :
    (f2A A B L alpha beta gamma n).support ⊆ ({L + Int.ofNat n, A + B + Int.ofNat (2 * n)} : Finset ℤ) := by
  classical
  have h := (Finsupp.support_add (g₁ := (C gamma * T (L + Int.ofNat n) : LPoly))
    (g₂ := (C (-(alpha * beta)) * T (A + B + Int.ofNat (2 * n)) : LPoly)))
  refine h.trans ?_
  refine Finset.union_subset_iff.2 ?_
  constructor
  · exact (LaurentPolynomial.support_C_mul_T (R := Rat) (a := gamma) (n := L + Int.ofNat n)).trans (by simp)
  ·
    exact
      (LaurentPolynomial.support_C_mul_T (R := Rat) (a := (-(alpha * beta)))
            (n := A + B + Int.ofNat (2 * n))).trans
        (by simp)

private lemma support_f2B_subset (A B : ℤ) (alpha beta : Rat) (n : Nat) :
    (f2B A B alpha beta n).support ⊆ ({0, B + Int.ofNat n, A + Int.ofNat n + 1} : Finset ℤ) := by
  classical
  -- Reassociate as `(C 1 + beta-term) + alpha-term`, then bound supports by unions of singletons.
  have h2 :=
    (Finsupp.support_add
      (g₁ := (C (1 : Rat) + C beta * T (B + Int.ofNat n) : LPoly))
      (g₂ := (C alpha * T (A + Int.ofNat n + 1) : LPoly)))
  have h1 :=
    (Finsupp.support_add
      (g₁ := (C (1 : Rat) : LPoly))
      (g₂ := (C beta * T (B + Int.ofNat n) : LPoly)))

  refine h2.trans ?_
  refine Finset.union_subset_iff.2 ?_
  constructor
  · -- Left part: `C 1 + beta-term`.
    refine h1.trans ?_
    refine Finset.union_subset_iff.2 ?_
    constructor
    · -- `C 1` is supported at exponent `0`.
      -- First rewrite `C 1` as a single at exponent `0`, then widen to the target finset.
      have h0 : (C (1 : Rat) : LPoly).support ⊆ ({(0 : ℤ)} : Finset ℤ) := by
        simpa [LaurentPolynomial.C, -LaurentPolynomial.single_eq_C_mul_T] using
          (Finsupp.support_single_subset (a := (0 : ℤ)) (b := (1 : Rat)))
      exact h0.trans (by simp)
    · -- `beta * T^(B+n)` is supported at exponent `B+n`.
      exact
        (LaurentPolynomial.support_C_mul_T (R := Rat) (a := beta) (n := B + Int.ofNat n)).trans
          (by simp)
  · -- Right part: `alpha * T^(A+n+1)` is supported at exponent `A+n+1`.
    exact
      (LaurentPolynomial.support_C_mul_T (R := Rat) (a := alpha) (n := A + Int.ofNat n + 1)).trans
        (by simp)

/-!
Stage-2 exclusions.
-/

/-- No `f4` monomial substitution can match the target at stages `n = 1,2`. -/
theorem no_f4_match_stage2 (A B L : ℤ) :
    ¬ ∃ alpha beta gamma : Rat,
        f4A A L alpha gamma 1 = targetA 1 ∧
        f4A A L alpha gamma 2 = targetA 2 ∧
        f4B A B alpha beta 1 = targetB 1 ∧
        f4B A B alpha beta 2 = targetB 2 := by
  classical
  intro h
  rcases h with ⟨alpha, beta, gamma, hA1, hA2, -, -⟩

  have hsA1 : (f4A A L alpha gamma 1).support = ({(1 : ℤ), (2 : ℤ)} : Finset ℤ) := by
    simpa [targetA_support_one] using congrArg (fun p : LPoly => p.support) hA1

  have hsA2 : (f4A A L alpha gamma 2).support = ({(2 : ℤ), (4 : ℤ)} : Finset ℤ) := by
    simpa [targetA_support_two] using congrArg (fun p : LPoly => p.support) hA2

  have h12 : ({(1 : ℤ), (2 : ℤ)} : Finset ℤ) ⊆ ({A + 1, L + 1} : Finset ℤ) := by
    -- stage 1 support is `{1,2}` and is bounded by `{A+1, L+1}`.
    have hsub := support_f4A_subset (A := A) (L := L) (alpha := alpha) (gamma := gamma) (n := 1)
    simpa [hsA1] using hsub

  have h24 : ({(2 : ℤ), (4 : ℤ)} : Finset ℤ) ⊆ ({A + 1, L + 2} : Finset ℤ) := by
    have hsub := support_f4A_subset (A := A) (L := L) (alpha := alpha) (gamma := gamma) (n := 2)
    simpa [hsA2] using hsub

  -- Extract `A+1 ∈ {1,2}` from `1 ∈ {A+1,L+1}` and `2 ∈ {A+1,L+1}` (and similarly for `{2,4}`).
  have hA_in_12 : A + 1 = (1 : ℤ) ∨ A + 1 = (2 : ℤ) := by
    have h1_mem : (1 : ℤ) ∈ ({A + 1, L + 1} : Finset ℤ) := h12 (by simp)
    have h2_mem : (2 : ℤ) ∈ ({A + 1, L + 1} : Finset ℤ) := h12 (by simp)
    have h1_or : (1 : ℤ) = A + 1 ∨ (1 : ℤ) = L + 1 := by
      simpa [Finset.mem_insert, Finset.mem_singleton] using h1_mem
    have h2_or : (2 : ℤ) = A + 1 ∨ (2 : ℤ) = L + 1 := by
      simpa [Finset.mem_insert, Finset.mem_singleton] using h2_mem
    rcases h1_or with hA1 | hL1
    · left; exact hA1.symm
    · rcases h2_or with hA2 | hL2
      · right; exact hA2.symm
      · exfalso; linarith [hL1, hL2]

  have hA_in_24 : A + 1 = (2 : ℤ) ∨ A + 1 = (4 : ℤ) := by
    have h2_mem : (2 : ℤ) ∈ ({A + 1, L + 2} : Finset ℤ) := h24 (by simp)
    have h4_mem : (4 : ℤ) ∈ ({A + 1, L + 2} : Finset ℤ) := h24 (by simp)
    have h2_or : (2 : ℤ) = A + 1 ∨ (2 : ℤ) = L + 2 := by
      simpa [Finset.mem_insert, Finset.mem_singleton] using h2_mem
    have h4_or : (4 : ℤ) = A + 1 ∨ (4 : ℤ) = L + 2 := by
      simpa [Finset.mem_insert, Finset.mem_singleton] using h4_mem
    rcases h2_or with hA2 | hL2
    · left; exact hA2.symm
    · rcases h4_or with hA4 | hL4
      · right; exact hA4.symm
      · exfalso; linarith [hL2, hL4]

  have hA_eq2 : A + 1 = (2 : ℤ) := by
    rcases hA_in_12 with hA1 | hA2
    · -- A+1=1 contradicts A+1∈{2,4}.
      rcases hA_in_24 with hA2' | hA4'
      · linarith [hA1, hA2']
      · linarith [hA1, hA4']
    · exact hA2

  -- Stage 1: with A+1=2, the only way to cover exponent 1 is L+1=1.
  have h1_mem : (1 : ℤ) ∈ ({A + 1, L + 1} : Finset ℤ) := h12 (by simp)
  have h1_or : (1 : ℤ) = A + 1 ∨ (1 : ℤ) = L + 1 := by
    simpa [Finset.mem_insert, Finset.mem_singleton] using h1_mem
  have hL_eq0 : L = 0 := by
    rcases h1_or with hA | hL
    · linarith [hA, hA_eq2]
    · linarith

  -- Stage 2 then forces `L+2 = 4`, contradicting `L=0`.
  have h4_mem : (4 : ℤ) ∈ ({A + 1, L + 2} : Finset ℤ) := h24 (by simp)
  have h4_or : (4 : ℤ) = A + 1 ∨ (4 : ℤ) = L + 2 := by
    simpa [Finset.mem_insert, Finset.mem_singleton] using h4_mem
  rcases h4_or with hA | hL2
  · linarith [hA, hA_eq2]
  · have : L = 2 := by linarith
    linarith [this, hL_eq0]

/-- No `f2` monomial substitution can match the target at stages `n = 1,2`. -/
theorem no_f2_match_stage2 (A B L : ℤ) :
    ¬ ∃ alpha beta gamma : Rat,
        f2A A B L alpha beta gamma 1 = targetA 1 ∧
        f2A A B L alpha beta gamma 2 = targetA 2 ∧
        f2B A B alpha beta 1 = targetB 1 ∧
        f2B A B alpha beta 2 = targetB 2 := by
  classical
  intro h
  rcases h with ⟨alpha, beta, gamma, hA1, hA2, hB1Eq, -⟩

  have hsA1 : (f2A A B L alpha beta gamma 1).support = ({(1 : ℤ), (2 : ℤ)} : Finset ℤ) := by
    simpa [targetA_support_one] using congrArg (fun p : LPoly => p.support) hA1

  have hsA2 : (f2A A B L alpha beta gamma 2).support = ({(2 : ℤ), (4 : ℤ)} : Finset ℤ) := by
    simpa [targetA_support_two] using congrArg (fun p : LPoly => p.support) hA2

  have h12 : ({(1 : ℤ), (2 : ℤ)} : Finset ℤ) ⊆ ({L + 1, A + B + 2} : Finset ℤ) := by
    have hsub := support_f2A_subset (A := A) (B := B) (L := L) (alpha := alpha) (beta := beta) (gamma := gamma) (n := 1)
    simpa [hsA1] using hsub

  have h24 : ({(2 : ℤ), (4 : ℤ)} : Finset ℤ) ⊆ ({L + 2, A + B + 4} : Finset ℤ) := by
    have hsub := support_f2A_subset (A := A) (B := B) (L := L) (alpha := alpha) (beta := beta) (gamma := gamma) (n := 2)
    simpa [hsA2] using hsub

  -- `L+1` must be `1` or `2`.
  have h1_mem : (1 : ℤ) ∈ ({L + 1, A + B + 2} : Finset ℤ) := h12 (by simp)
  have h2_mem : (2 : ℤ) ∈ ({L + 1, A + B + 2} : Finset ℤ) := h12 (by simp)
  have h1_or : (1 : ℤ) = L + 1 ∨ (1 : ℤ) = A + B + 2 := by
    simpa [Finset.mem_insert, Finset.mem_singleton] using h1_mem
  have h2_or : (2 : ℤ) = L + 1 ∨ (2 : ℤ) = A + B + 2 := by
    simpa [Finset.mem_insert, Finset.mem_singleton] using h2_mem

  -- `L+2` must be `2` or `4` (otherwise `{2,4}` can't fit in a 2-element set).
  have h2_mem' : (2 : ℤ) ∈ ({L + 2, A + B + 4} : Finset ℤ) := h24 (by simp)
  have h4_mem' : (4 : ℤ) ∈ ({L + 2, A + B + 4} : Finset ℤ) := h24 (by simp)
  have h2_or' : (2 : ℤ) = L + 2 ∨ (2 : ℤ) = A + B + 4 := by
    simpa [Finset.mem_insert, Finset.mem_singleton] using h2_mem'
  have h4_or' : (4 : ℤ) = L + 2 ∨ (4 : ℤ) = A + B + 4 := by
    simpa [Finset.mem_insert, Finset.mem_singleton] using h4_mem'

  have hL1 : L + 1 = (1 : ℤ) ∨ L + 1 = (2 : ℤ) := by
    rcases h1_or with hL | hAB
    · left; exact hL.symm
    · rcases h2_or with hL | hAB2
      · right; exact hL.symm
      · exfalso; linarith [hAB, hAB2]

  have hL2 : L + 2 = (2 : ℤ) ∨ L + 2 = (4 : ℤ) := by
    rcases h2_or' with hL | hAB
    · left; exact hL.symm
    · rcases h4_or' with hL | hAB4
      · right; exact hL.symm
      · exfalso; linarith [hAB, hAB4]

  have hL_eq0 : L = 0 := by
    -- If `L+2 = 4` then `L+1 = 3`, contradicting `L+1 ∈ {1,2}`.
    rcases hL2 with hL2_eq2 | hL2_eq4
    · linarith
    · have : L + 1 = (3 : ℤ) := by linarith [hL2_eq4]
      rcases hL1 with hL1_eq1 | hL1_eq2
      · linarith [this, hL1_eq1]
      · linarith [this, hL1_eq2]

  have hAB_eq0 : A + B = 0 := by
    -- With L=0, the stage-1 exponents are `{1, A+B+2}` and must cover `{1,2}`.
    rcases h2_or with h2L | h2AB
    · -- 2 = L+1 contradicts L=0.
      linarith [hL_eq0, h2L]
    · linarith

  -- Now use the stage-1 denominator match to derive a contradiction.
  have hsB1 : (f2B A B alpha beta 1).support = ({(0 : ℤ), (1 : ℤ)} : Finset ℤ) := by
    simpa [targetB_support_one] using congrArg (fun p : LPoly => p.support) hB1Eq

  have h01 :
      ({(0 : ℤ), (1 : ℤ)} : Finset ℤ) ⊆ ({0, B + Int.ofNat 1, A + Int.ofNat 1 + 1} : Finset ℤ) := by
    have hsub := support_f2B_subset (A := A) (B := B) (alpha := alpha) (beta := beta) (n := 1)
    -- Rewrite the left-hand side using the known support equality.
    -- Use `rw` instead of `simp` to avoid normalizing the finset subset into arithmetic disjunctions.
    have hsub' := hsub
    rw [hsB1] at hsub'
    exact hsub'

  have h1_memB : (1 : ℤ) ∈ ({0, B + Int.ofNat 1, A + Int.ofNat 1 + 1} : Finset ℤ) := h01 (by simp)
  have hBA_cases : B = 0 ∨ A + 1 = 0 := by
    -- Membership of `1` means `1 = B+1` or `1 = A+2` (the `1=0` case is impossible).
    have : (1 : ℤ) = 0 ∨ (1 : ℤ) = B + Int.ofNat 1 ∨ (1 : ℤ) = A + Int.ofNat 1 + 1 := by
      simpa [Finset.mem_insert, Finset.mem_singleton] using h1_memB
    rcases this with h10 | hrest
    · exfalso; exact (by decide : (1 : ℤ) ≠ 0) h10
    rcases hrest with hB1 | hA2
    ·
      -- `1 = B+1` forces `B = 0`.
      left
      have h' : B + (1 : ℤ) = (1 : ℤ) := by simpa using hB1.symm
      have : B = 0 := by
        have := congrArg (fun x : ℤ => x - 1) h'
        simpa using this
      exact this
    ·
      -- `1 = A+2` forces `A+1 = 0`.
      right
      have h' : A + Int.ofNat 1 + (1 : ℤ) = (1 : ℤ) := by simpa using hA2.symm
      have := congrArg (fun x : ℤ => x - 1) h'
      -- simplify `(A+1+1)-1 = A+1` and `1-1 = 0`.
      simpa using this

  rcases hBA_cases with hB0 | hA1eq0
  · -- B = 0, so A = 0 (since A+B=0), but then the `alpha` term sits at exponent 2.
    have hA0 : A = 0 := by linarith [hAB_eq0, hB0]
    -- Coefficient at exponent 2 on the numerator stage-1 match is exactly `-(alpha*beta)` and equals 1.
    have hCoeffA2 : (f2A A B L alpha beta gamma 1) 2 = (targetA 1) 2 := by
      simpa using congrArg (fun p : LPoly => p 2) hA1
    have hMul : alpha * beta = -1 := by
      -- Unfold at exponent 2 using `L=0` and `A+B=0`.
      have hCoeff : -((C alpha * C beta : LPoly) 0) = (1 : Rat) := by
        simpa [f2A, targetA, hL_eq0, hAB_eq0, C_mul_T_eq_single, LaurentPolynomial.T,
          -LaurentPolynomial.single_eq_C_mul_T] using hCoeffA2
      have : -(alpha * beta) = (1 : Rat) := by
        simpa [C_mul_C_apply_zero] using hCoeff
      linarith
    -- Coefficient at exponent 2 in `b₁` must be 0 on the target side, hence forces `alpha = 0`.
    have hCoeffB2 : (f2B A B alpha beta 1) 2 = (targetB 1) 2 := by
      simpa using congrArg (fun p : LPoly => p 2) hB1Eq
    have hAlpha0 : alpha = 0 := by
      simpa [f2B, targetB, hA0, hB0, C_mul_T_eq_single, LaurentPolynomial.C, LaurentPolynomial.T,
        -LaurentPolynomial.single_eq_C_mul_T] using hCoeffB2
    -- Contradiction: `alpha = 0` would make `alpha*beta = 0`, but we already have `alpha*beta = -1`.
    have hneq : (0 : Rat) ≠ (-1 : Rat) := by norm_num
    have : (0 : Rat) = (-1 : Rat) := by
      simpa [hAlpha0] using hMul
    exact hneq this
  · -- A+1=0 implies A = -1; then B = 1 from A+B=0, so the `beta` term sits at exponent 2.
    have hAneg1 : A = -1 := by linarith [hA1eq0]
    have hB1' : B = 1 := by linarith [hAB_eq0, hAneg1]
    have hCoeffA2 : (f2A A B L alpha beta gamma 1) 2 = (targetA 1) 2 := by
      simpa using congrArg (fun p : LPoly => p 2) hA1
    have hMul : alpha * beta = -1 := by
      have hCoeff : -((C alpha * C beta : LPoly) 0) = (1 : Rat) := by
        simpa [f2A, targetA, hL_eq0, hAB_eq0, C_mul_T_eq_single, LaurentPolynomial.T,
          -LaurentPolynomial.single_eq_C_mul_T] using hCoeffA2
      have : -(alpha * beta) = (1 : Rat) := by
        simpa [C_mul_C_apply_zero] using hCoeff
      linarith
    have hCoeffB2 : (f2B A B alpha beta 1) 2 = (targetB 1) 2 := by
      simpa using congrArg (fun p : LPoly => p 2) hB1Eq
    have hBeta0 : beta = 0 := by
      simpa [f2B, targetB, hAneg1, hB1', C_mul_T_eq_single, LaurentPolynomial.C, LaurentPolynomial.T,
        -LaurentPolynomial.single_eq_C_mul_T] using hCoeffB2
    have hneq : (0 : Rat) ≠ (-1 : Rat) := by norm_num
    have : (0 : Rat) = (-1 : Rat) := by
      simpa [hBeta0] using hMul
    exact hneq this

end Page43

end

end HeroCase
end Proofs
