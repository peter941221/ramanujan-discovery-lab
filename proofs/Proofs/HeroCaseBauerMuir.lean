import Proofs.HeroCaseObjects

namespace Proofs
namespace HeroCase

/-!
Finite-step Bauer–Muir transform class exclusion (bounded, computation-checked).

The research pipeline defines a tiny fixed family of low-complexity Bauer–Muir modifiers `wₙ`
and exhaustively checks 1/2/3-step chains against exact rational sample points.

This file mirrors that computation in Lean over `Rat` and proves there are no hits for:

* RR reciprocal source → hero target
* cubic reciprocal source → hero target

at depth `4`, checked at `t = 1/10, 1/7, 1/5`.
-/

namespace BauerMuir

open Proofs

def samples : List Rat :=
  [(1 : Rat) / 10, (1 : Rat) / 7, (1 : Rat) / 5]

structure Trunc where
  b0 : Rat
  a : Array Rat
  b : Array Rat
deriving DecidableEq

def truncOfData (g : GCFData Rat) (depth : Nat) : Trunc :=
  let aTerms :=
    Array.ofFn (fun i : Fin (depth + 1) =>
      if i.val = 0 then
        0
      else
        g.a (i.val - 1))
  let bTerms :=
    Array.ofFn (fun i : Fin (depth + 1) =>
      if i.val = 0 then
        0
      else
        g.b (i.val - 1))
  { b0 := g.b0, a := aTerms, b := bTerms }

def patterns : Array (Nat → Rat → Rat) :=
  let base : List (Nat → Rat → Rat) := [fun _ _ => 0]
  let scales : List Rat := [(-2 : Rat), (-1 : Rat), (1 : Rat), (2 : Rat)]
  let others : List (Nat → Rat → Rat) :=
    List.flatMap
      (fun s =>
        [ (fun n q => s * (q ^ n - 1)),
          (fun n q => s * (q ^ (2 * n) - 1)),
          (fun n q => s * (q ^ n - q ^ (2 * n))) ])
      scales
  Array.mk (base ++ others)

def wTerms (wFn : Nat → Rat → Rat) (t : Rat) (depth : Nat) : Array Rat :=
  Array.ofFn (fun i : Fin (depth + 1) => wFn i.val t)

def bauerMuirTransformTrunc (x : Trunc) (w : Array Rat) : Option Trunc :=
  let depth := x.a.size - 1
  if w.size != depth + 1 then
    none
  else
    let lam : Array Rat :=
      Array.ofFn (fun i : Fin (depth + 1) =>
        match i.val with
        | 0 => 0
        | n + 1 =>
            let a_n := x.a[n + 1]!
            let b_n := x.b[n + 1]!
            let w_prev := w[n]!
            let w_n := w[n + 1]!
            a_n - w_prev * (b_n + w_n))

    -- Avoid dividing by 0 in later stages.
    let badDenom := (List.range (depth + 1)).drop 2 |>.any (fun n =>
      -- check `lam[n-1] = 0` for n>=2
      decide (lam[n - 1]! = 0))
    if badDenom then
      none
    else
      let b0' := x.b0 + w[0]!
      let a' : Array Rat :=
        Array.ofFn (fun i : Fin (depth + 1) =>
          match i.val with
          | 0 => 0
          | 1 => lam[1]!
          | n + 2 =>
              let denom := lam[n + 1]!
              x.a[n + 1]! * lam[n + 2]! / denom)
      let b' : Array Rat :=
        Array.ofFn (fun i : Fin (depth + 1) =>
          match i.val with
          | 0 => 0
          | 1 => x.b[1]! + w[1]!
          | n + 2 =>
              let denom := lam[n + 1]!
              x.b[n + 2]! + w[n + 2]! - w[n]! * lam[n + 2]! / denom)
      some { b0 := b0', a := a', b := b' }

def truncMatches (x y : Trunc) : Bool :=
  decide (x = y)

def applyChainAt (source target : Rat → GCFData Rat) (depth : Nat) (t : Rat) (chain : List Nat) : Bool :=
  let targetTrunc := truncOfData (target t) depth
  let rec go (current : Trunc) : List Nat → Bool
    | [] => truncMatches current targetTrunc
    | idx :: rest =>
        let wFn := patterns[idx]!
        let w := wTerms wFn t depth
        match bauerMuirTransformTrunc current w with
        | none => false
        | some next => go next rest
  go (truncOfData (source t) depth) chain

def chainMatches (source target : Rat → GCFData Rat) (depth : Nat) (chain : List Nat) : Bool :=
  samples.all (fun t => applyChainAt source target depth t chain)

def allChains (count : Nat) : Nat → List (List Nat)
  | 0 => [[]]
  | steps + 1 =>
      let idxs := List.range count
      let tails := allChains count steps
      List.flatMap (fun i => tails.map (fun tail => i :: tail)) idxs

def hitChains (source target : Rat → GCFData Rat) (depth steps : Nat) : List (List Nat) :=
  (allChains patterns.size steps).filter (fun chain => chainMatches source target depth chain)

def rrHits (steps : Nat) : List (List Nat) :=
  hitChains (source := rrReciprocalDataGeneric) (target := heroDataGeneric) (depth := 4) (steps := steps)

def cubicHits (steps : Nat) : List (List Nat) :=
  hitChains (source := cubicReciprocalDataGeneric) (target := heroDataGeneric) (depth := 4) (steps := steps)

theorem rrHits_step1_eq_nil : rrHits 1 = [] := by
  native_decide

theorem rrHits_step2_eq_nil : rrHits 2 = [] := by
  native_decide

theorem rrHits_step3_eq_nil : rrHits 3 = [] := by
  native_decide

theorem cubicHits_step1_eq_nil : cubicHits 1 = [] := by
  native_decide

theorem cubicHits_step2_eq_nil : cubicHits 2 = [] := by
  native_decide

theorem cubicHits_step3_eq_nil : cubicHits 3 = [] := by
  native_decide

end BauerMuir

end HeroCase
end Proofs
