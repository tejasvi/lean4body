import C4.SearchSound
import C4Cert.Cover

/-!
# The two caps

The certificates (`C4Cert`) check by `decide +kernel` that `checkBoxH` succeeds, with the recorded
hints, on every cell of the two grids of `C4Check/Grid.lean` (`dirCover`, `chCover`);
`checkBoxH_sound` turns each success into a statement about the real functions, and the cells
cover the two regions.

* `dir_cap`: on `[1/2, 7/4] × [1/8, 1] × [0, 7/4]`, every feasible zero of `dirP` with `|x| ≤ 1`
  has a `y` with `tr S(y) < 3/4`;
* `ch_cap`: the same on `[0, 1/8] × [-2, 1] × [0, 3]` for the blow-up chart.
-/

namespace C4

noncomputable section

/-! ## grid cells -/

theorem mem_mk {v : ℝ} {L H : ℕ} (h1 : -v ≤ val L) (h2 : v ≤ val H) : I.mem v ⟨L, H⟩ :=
  ⟨h1, h2⟩

theorem val_subBb {m : ℕ} (hm : m ≤ Bb) : val (Nat.sub Bb m) = -(m : ℝ) / 2 ^ 96 := by
  unfold val
  rw [cast_sub_of_le hm, cBb]
  ring

theorem val_addsubBb {b c : ℕ} (h : c ≤ Nat.add Bb b) :
    val (Nat.sub (Nat.add Bb b) c) = ((b : ℝ) - c) / 2 ^ 96 := by
  unfold val
  rw [cast_sub_of_le h, cast_add, cBb]
  ring

theorem mul_pow_le {n k : ℕ} (hn : n ≤ 1000) (hk : k ≤ 91) : Nat.mul n (2 ^ k) ≤ Bb := by
  change n * 2 ^ k ≤ 2 ^ 256
  calc n * 2 ^ k ≤ 1000 * 2 ^ 91 :=
        Nat.mul_le_mul hn (Nat.pow_le_pow_right (by norm_num) hk)
    _ ≤ 2 ^ 256 := by norm_num

theorem mem_iv32 {t : ℝ} {n : ℕ} (hn : n ≤ 1000) (h1 : (n : ℝ) ≤ 32 * t)
    (h2 : 32 * t ≤ n + 1) : I.mem t (iv32 n) := by
  rw [iv32]
  refine mem_mk ?_ ?_
  · rw [val_subBb (mul_pow_le hn (le_refl 91)), cast_mul, Nat.cast_pow, Nat.cast_ofNat]
    have e : -((n : ℝ) * 2 ^ 91) / 2 ^ 96 = -(n / 32) := by ring
    rw [e]
    linarith
  · rw [val_addB', cast_mul, Nat.cast_succ, Nat.cast_pow, Nat.cast_ofNat]
    have e : ((n : ℝ) + 1) * 2 ^ 91 / 2 ^ 96 = (n + 1) / 32 := by ring
    rw [e]
    linarith

theorem mem_iv64 {t : ℝ} {n : ℕ} (hn : n ≤ 1000) (h1 : (n : ℝ) ≤ 64 * t)
    (h2 : 64 * t ≤ n + 1) : I.mem t (iv64 n) := by
  rw [iv64]
  refine mem_mk ?_ ?_
  · rw [val_subBb (mul_pow_le hn (by norm_num)), cast_mul, Nat.cast_pow, Nat.cast_ofNat]
    have e : -((n : ℝ) * 2 ^ 90) / 2 ^ 96 = -(n / 64) := by ring
    rw [e]
    linarith
  · rw [val_addB', cast_mul, Nat.cast_succ, Nat.cast_pow, Nat.cast_ofNat]
    have e : ((n : ℝ) + 1) * 2 ^ 90 / 2 ^ 96 = (n + 1) / 64 := by ring
    rw [e]
    linarith

theorem mem_ivAl {t : ℝ} {j : ℕ} (hj : j ≤ 1000) (h1 : (j : ℝ) ≤ 32 * t + 64)
    (h2 : 32 * t + 64 ≤ j + 1) : I.mem t (ivAl j) := by
  rw [ivAl]
  refine mem_mk ?_ ?_
  · rw [val_addsubBb (le_trans (mul_pow_le hj (le_refl 91)) (Nat.le_add_right _ _)), cast_mul,
      cast_mul, Nat.cast_pow, Nat.cast_ofNat, Nat.cast_ofNat]
    have e : ((64 : ℝ) * 2 ^ 91 - j * 2 ^ 91) / 2 ^ 96 = (64 - j) / 32 := by ring
    rw [e]
    linarith
  · rw [val_addsubBb (le_trans (mul_pow_le (by norm_num) (le_refl 91)) (Nat.le_add_right _ _)),
      cast_mul, cast_mul, Nat.cast_pow, Nat.cast_ofNat, Nat.cast_ofNat]
    have e : ((j : ℝ) * 2 ^ 91 - 63 * 2 ^ 91) / 2 ^ 96 = (j - 63) / 32 := by ring
    rw [e]
    linarith

/-- a point of `[0, N]` lies in one of the cells `[i, i + 1]`, `i < N` -/
theorem cell_of {t : ℝ} {N : ℕ} (hN : 0 < N) (h0 : 0 ≤ t) (h1 : t ≤ N) :
    ∃ i : ℕ, i < N ∧ (i : ℝ) ≤ t ∧ t ≤ i + 1 := by
  refine ⟨min ⌊t⌋₊ (N - 1), by omega, ?_, ?_⟩
  · calc ((min ⌊t⌋₊ (N - 1) : ℕ) : ℝ) ≤ (⌊t⌋₊ : ℝ) := by exact_mod_cast min_le_left _ _
      _ ≤ t := Nat.floor_le h0
  · rcases le_total ⌊t⌋₊ (N - 1) with h | h
    · rw [min_eq_left h]
      exact (Nat.lt_floor_add_one t).le
    · rw [min_eq_right h]
      have : ((N - 1 : ℕ) : ℝ) = (N : ℝ) - 1 := by rw [Nat.cast_sub (by omega), Nat.cast_one]
      rw [this]
      linarith

/-! ## the direct grid -/

theorem dirBox_mem {i j k : ℕ} {a b c : ℝ} (hi : i < 40) (hj : j < 28) (hk : k < 56)
    (ha1 : (i : ℝ) ≤ 32 * a - 16) (ha2 : 32 * a - 16 ≤ i + 1)
    (hb1 : (j : ℝ) ≤ 32 * b - 4) (hb2 : 32 * b - 4 ≤ j + 1)
    (hc1 : (k : ℝ) ≤ 32 * c) (hc2 : 32 * c ≤ k + 1) : (dirBox i j k).mem a b c :=
  ⟨mem_iv32 (by omega) (by push_cast; linarith) (by push_cast; linarith),
    mem_iv32 (by omega) (by push_cast; linarith) (by push_cast; linarith),
    mem_iv32 (by omega) hc1 hc2⟩

/-- the cap of the direct chart, given the cover -/
theorem dir_cap_of (hC : Cover dirCell 0 62720) {a b c x : ℝ} (ha1 : 1 / 2 ≤ a) (ha2 : a ≤ 7 / 4)
    (hb1 : 1 / 8 ≤ b) (hb2 : b ≤ 1) (hc1 : 0 ≤ c) (hc2 : c ≤ 7 / 4) (hx1 : -1 ≤ x) (hx2 : x ≤ 1)
    (hG : ∀ i, 0 ≤ (dirGR a b c x).get i) (hP : dirPR a b c x = 0) :
    ∃ y0 y1, trSR (dirCertR a b c x) y0 y1 < 3 / 4 := by
  obtain ⟨i, hi, hi1, hi2⟩ := cell_of (N := 40) (t := 32 * a - 16) (by norm_num) (by linarith)
    (by push_cast; linarith)
  obtain ⟨j, hj, hj1, hj2⟩ := cell_of (N := 28) (t := 32 * b - 4) (by norm_num) (by linarith)
    (by push_cast; linarith)
  obtain ⟨k, hk, hk1, hk2⟩ := cell_of (N := 56) (t := 32 * c) (by norm_num) (by linarith)
    (by push_cast; linarith)
  obtain ⟨h, hh⟩ := hC ((28 * i + j) * 56 + k) (Nat.zero_le _) (by omega)
  have e1 : ((28 * i + j) * 56 + k) / 1568 = i := by omega
  have e2 : ((28 * i + j) * 56 + k) / 56 % 28 = j := by omega
  have e3 : ((28 * i + j) * 56 + k) % 56 = k := by omega
  rw [dirCell, e1, e2, e3] at hh
  obtain ⟨h', hh'⟩ := Option.isSome_iff_exists.mp hh
  exact checkBoxH_sound dirMode_sound depth _ h h' hh' a b c x
    (dirBox_mem hi hj hk hi1 hi2 hj1 hj2 hk1 hk2) (fun _ => ⟨hx1, hx2⟩) hG hP

/-- the cap of the direct chart -/
theorem dir_cap {a b c x : ℝ} (ha1 : 1 / 2 ≤ a) (ha2 : a ≤ 7 / 4) (hb1 : 1 / 8 ≤ b) (hb2 : b ≤ 1)
    (hc1 : 0 ≤ c) (hc2 : c ≤ 7 / 4) (hx1 : -1 ≤ x) (hx2 : x ≤ 1)
    (hG : ∀ i, 0 ≤ (dirGR a b c x).get i) (hP : dirPR a b c x = 0) :
    ∃ y0 y1, trSR (dirCertR a b c x) y0 y1 < 3 / 4 :=
  dir_cap_of dirCover ha1 ha2 hb1 hb2 hc1 hc2 hx1 hx2 hG hP

/-! ## the blow-up grid -/

theorem chBox_mem {i j k : ℕ} {b al ga : ℝ} (hi : i < 8) (hj : j < 96) (hk : k < 96)
    (h1 : (i : ℝ) ≤ 64 * b) (h2 : 64 * b ≤ i + 1)
    (h3 : (j : ℝ) ≤ 32 * al + 64) (h4 : 32 * al + 64 ≤ j + 1)
    (h5 : (k : ℝ) ≤ 32 * ga) (h6 : 32 * ga ≤ k + 1) : (chBox i j k).mem b al ga :=
  ⟨mem_iv64 (by omega) h1 h2, mem_ivAl (by omega) h3 h4, mem_iv32 (by omega) h5 h6⟩

/-- the cap of the blow-up chart, given the cover -/
theorem ch_cap_of (hC : Cover chCell 0 73728) {b al ga xi : ℝ} (hb1 : 0 ≤ b) (hb2 : b ≤ 1 / 8)
    (hal1 : -2 ≤ al) (hal2 : al ≤ 1) (hga1 : 0 ≤ ga) (hga2 : ga ≤ 3)
    (hG : ∀ i, 0 ≤ (chGR b al ga xi).get i) (hP : chPR b al ga xi = 0) :
    ∃ Y0 Y1, trSR (chCertR b al ga xi) Y0 Y1 < 3 / 4 := by
  obtain ⟨i, hi, hi1, hi2⟩ := cell_of (N := 8) (t := 64 * b) (by norm_num) (by linarith)
    (by push_cast; linarith)
  obtain ⟨j, hj, hj1, hj2⟩ := cell_of (N := 96) (t := 32 * al + 64) (by norm_num) (by linarith)
    (by push_cast; linarith)
  obtain ⟨k, hk, hk1, hk2⟩ := cell_of (N := 96) (t := 32 * ga) (by norm_num) (by linarith)
    (by push_cast; linarith)
  obtain ⟨h, hh⟩ := hC ((96 * i + j) * 96 + k) (Nat.zero_le _) (by omega)
  have e1 : ((96 * i + j) * 96 + k) / 9216 = i := by omega
  have e2 : ((96 * i + j) * 96 + k) / 96 % 96 = j := by omega
  have e3 : ((96 * i + j) * 96 + k) % 96 = k := by omega
  rw [chCell, e1, e2, e3] at hh
  obtain ⟨h', hh'⟩ := Option.isSome_iff_exists.mp hh
  exact checkBoxH_sound chMode_sound depth _ h h' hh' b al ga xi
    (chBox_mem hi hj hk hi1 hi2 hj1 hj2 hk1 hk2) (fun h => absurd h Bool.false_ne_true) hG hP

/-- the cap of the blow-up chart -/
theorem ch_cap {b al ga xi : ℝ} (hb1 : 0 ≤ b) (hb2 : b ≤ 1 / 8) (hal1 : -2 ≤ al) (hal2 : al ≤ 1)
    (hga1 : 0 ≤ ga) (hga2 : ga ≤ 3) (hG : ∀ i, 0 ≤ (chGR b al ga xi).get i)
    (hP : chPR b al ga xi = 0) :
    ∃ Y0 Y1, trSR (chCertR b al ga xi) Y0 Y1 < 3 / 4 :=
  ch_cap_of chCover hb1 hb2 hal1 hal2 hga1 hga2 hG hP

end

end C4
