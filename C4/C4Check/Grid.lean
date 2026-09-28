import C4Check.Search

/-!
# The two charts and their initial grids

Direct chart: `(a, b, c) ∈ [1/2, 7/4] × [1/8, 1] × [0, 7/4]`, cut into `40 × 28 × 56` cubes of
side `1/32`; the dependent coordinate `x = cos θ` is clamped to `[-1, 1]`.

Blow-up chart: `(b, al, ga) ∈ [0, 1/8] × [-2, 1] × [0, 3]`, cut into `8 × 96 × 96` boxes of
sides `1/64 × 1/32 × 1/32`; the dependent coordinate is `xi`.

The cells are numbered `n = (i · J + j) · K + k`.  `allCells P n₀ n₁ hs` checks the cells
`n₀ ≤ n < n₁` with the hints `hs`; the certificates in `C4Cert` evaluate it with
`decide +kernel`, and chain the results into `Cover P 0 N` with `allCells_cover` and
`Cover.trans`.
-/

namespace C4

def dirMode : Mode where
  kl A B C := (dirK (some A) (some B) (some C), dirL (some A) (some B) (some C))
  gP R z0 z1 z2 z3 :=
    (@dirG _ (jetOps R) (some z0) (some z1) (some z2) (some z3),
     @dirP _ (jetOps R) (some z0) (some z1) (some z2) (some z3))
  trJ R z0 z1 z2 z3 y0 y1 :=
    @trS _ (jetOps R) (@dirCert _ (jetOps R) (some z0) (some z1) (some z2) (some z3))
      (some (J.cst y0)) (some (J.cst y1))
  clamp := true

def chMode : Mode where
  kl A B C := (chK (some A) (some B) (some C), chL (some A) (some B) (some C))
  gP R z0 z1 z2 z3 :=
    (@chG _ (jetOps R) (some z0) (some z1) (some z2) (some z3),
     @chP _ (jetOps R) (some z0) (some z1) (some z2) (some z3))
  trJ R z0 z1 z2 z3 y0 y1 :=
    @trS _ (jetOps R) (@chCert _ (jetOps R) (some z0) (some z1) (some z2) (some z3))
      (some (J.cst y0)) (some (J.cst y1))
  clamp := false

/-- `[n/32, (n+1)/32]` for `0 ≤ n` -/
def iv32 (n : Nat) : I := ⟨Nat.sub Bb (Nat.mul n (2 ^ 91)), Nat.add Bb (Nat.mul (Nat.succ n) (2 ^ 91))⟩

/-- `[n/64, (n+1)/64]` for `0 ≤ n` -/
def iv64 (n : Nat) : I := ⟨Nat.sub Bb (Nat.mul n (2 ^ 90)), Nat.add Bb (Nat.mul (Nat.succ n) (2 ^ 90))⟩

/-- `[(j - 64)/32, (j - 63)/32]` for `j < 96` -/
def ivAl (j : Nat) : I :=
  ⟨Nat.sub (Nat.add Bb (Nat.mul 64 (2 ^ 91))) (Nat.mul j (2 ^ 91)),
    Nat.sub (Nat.add Bb (Nat.mul j (2 ^ 91))) (Nat.mul 63 (2 ^ 91))⟩

/-- `[(16+i)/32, (17+i)/32] × [(4+j)/32, (5+j)/32] × [k/32, (k+1)/32]` -/
def dirBox (i j k : Nat) : Box := ⟨iv32 (16 + i), iv32 (4 + j), iv32 k⟩

/-- `[i/64, (i+1)/64] × [(j-64)/32, (j-63)/32] × [k/32, (k+1)/32]` -/
def chBox (i j k : Nat) : Box := ⟨iv64 i, ivAl j, iv32 k⟩

/-- depth bound of the recursion -/
def depth : Nat := 100

/-- cell `n` of the direct grid (`n = (28 i + j) 56 + k`) is certified by the hint `h` -/
def dirCell (n h : Nat) : Bool :=
  (checkBoxH dirMode depth (dirBox (n / 1568) (n / 56 % 28) (n % 56)) h).isSome

/-- cell `n` of the blow-up grid (`n = (96 i + j) 96 + k`) is certified by the hint `h` -/
def chCell (n h : Nat) : Bool :=
  (checkBoxH chMode depth (chBox (n / 9216) (n / 96 % 96) (n % 96)) h).isSome

/-- `P n h_n` for `n₀ ≤ n < n₁`, where `hs = [h_{n₀}, h_{n₀+1}, …]` -/
def allCells (P : Nat → Nat → Bool) : Nat → Nat → List Nat → Bool
  | n, e, [] => Nat.ble e n
  | n, e, h :: hs => sel (P n h) (allCells P (Nat.succ n) e hs) false

/-! ## assembling the certificates -/

/-- every cell `a ≤ n < b` has a hint that certifies it -/
def Cover (P : Nat → Nat → Bool) (a b : Nat) : Prop := ∀ n, a ≤ n → n < b → ∃ h, P n h = true

theorem allCells_cover {P : Nat → Nat → Bool} :
    ∀ {a b : Nat} {hs : List Nat}, allCells P a b hs = true → Cover P a b
  | a, b, [], h => by
    have hb : b ≤ a := Nat.le_of_ble_eq_true h
    intro n h1 h2
    exact absurd (Nat.lt_of_lt_of_le h2 hb) (Nat.not_lt.mpr h1)
  | a, b, x :: hs, h => by
    cases hP : P a x with
    | false =>
      have h' : allCells P a b (x :: hs) = false := by
        show Bool.casesOn (motive := fun _ => Bool) (P a x) false (allCells P (Nat.succ a) b hs)
          = false
        rw [hP]
      rw [h'] at h
      exact absurd h Bool.false_ne_true
    | true =>
      have h' : allCells P a b (x :: hs) = allCells P (Nat.succ a) b hs := by
        show Bool.casesOn (motive := fun _ => Bool) (P a x) false (allCells P (Nat.succ a) b hs)
          = allCells P (Nat.succ a) b hs
        rw [hP]
      rw [h'] at h
      have ih := allCells_cover h
      intro n h1 h2
      cases Nat.eq_or_lt_of_le h1 with
      | inl e => exact ⟨x, e ▸ hP⟩
      | inr e => exact ih n e h2

theorem Cover.trans {P : Nat → Nat → Bool} {a b c : Nat} (h1 : Cover P a b) (h2 : Cover P b c) :
    Cover P a c := fun n ha hc =>
  if hb : n < b then h1 n ha hb else h2 n (Nat.not_lt.mp hb) hc

end C4
