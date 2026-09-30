# Flirt Training Pipeline — fixed plan (v2)

Goal: analyse each message to teach an LLM to flirt like the operator.
Derived from deep-dive review; the seven defects it fixes are recorded at the end.

## Teaching mechanism: PREFERENCE PAIRS
The artifact format serves preference-style teaching (chosen vs rejected on the
same cue). `flirt_pairs.jsonl` is the primary signal; positives/negatives are
its raw material. A pair is valid only when chosen and rejected respond to the
SAME cue — cross-case pairs are marked cross_case.

## Two tracks, always separate
- attraction/ — the flirt itself. Gold case: gabby (fully annotated, operator-
  corrected). Katherine/Myah/Tina next per evidence class.
- execution/ — confirmation, logistics, timing, repair. Gold case: jisu (plan
  transmitted, date completed). The gabby date-day sequence is a pure negative
  exemplar on this track: the strongest attraction thread in the corpus still
  died at execution. Teaching them fused teaches failure as style.
Every exemplar carries track. Every pair is track-tagged.

## Unit of analysis: the MOVE
A move = the minimal unit that changes the relationship state, usually 1-3
bubbles, with before-state and after-state. Single bubbles like "Back and bis"
are inert alone — the pair with the backyard line is the move. Never annotate
or teach a fragment as if it were the move.

## Procedure traces
`procedure_traces.jsonl` preserves the live candidate→reject→refine loops from
the gabby transcript (5 traces, 20 steps). These are direct evidence of the
search PROCEDURE, richer than outcome labels. The 10-step algorithm in
flirt_model.md was reverse-engineered from outcomes; the traces show the
procedure happening. Both are taught.

## Evidence classes drive depth
(a) full transcripts → full schema (gabby: done).
(b) conversation sequences embedded in the audit (tina, parts of sophia/jasmine)
    → full schema possible today.
(c) solo sent lists → mechanism-level only, cue_class null, observed_return
    UNKNOWN where undocumented. `audit_minimal_pass.jsonl` covers all 200 audit
    rows at this floor — "each message" is literally true at the minimum grade.
NEVER compose an observed return. Extraction only.

## Hold-out identity test
`evals/test_identity_holdout.py` — 25 operator-sent lines from cases outside the
annotation corpus must not trip any critic kill. If a kill rejects his own
verified lines, the kill is overfit to gabby-shaped data. Currently: PASS.

## Graduation gate (D7)
No annotation's rules or exemplars enter the training files until the operator
has seen and confirmed/corrected them. The gabby corpus is operator-corrected
(that is what makes it true). New annotations are pending until reviewed.
KATHERINE FORK: CLOSED. Full message-by-message annotation supplied and stored at
training/cases/katherine_analysis.md (59 messages, A/B/C-graded). Phase 2 cross-case
validation COMPLETE — all 8 Gaby rules confirmed; refinements + 1 overfit kill removed
('prove it', state-dependent) written into flirt_model.md. Class-(b) threads DONE: tina, sophia, jasmine annotated at compact grade
(training/cases/), 6 exemplars added, ownership-ambiguous lines quarantined,
undocumented returns UNKNOWN. status: PENDING OPERATOR VERIFICATION per the gate.
Audit-holding cases DONE: aditi (real writeup — the vulnerability sequence as
positive model + neg-aditi contrast), abby, arielle, mahi (compact grade). 6
exemplars added. All PENDING OPERATOR VERIFICATION.
The audit's verified-sent corpus is now FULLY ANNOTATED at some depth — every
case in the provenance audit has a case file or an audit-grade entry.
Next: operator review pass of all pending annotations (the graduation gate),
or further material beyond the audit.

## Fixed-defects register
D1 teaching mechanism unspecified → preference pairs, built.
D2 fused tracks → two tracks, tagged; jisu gold identified.
D3 message-unit → move-unit, before/after state.
D4 outcomes-only → procedure traces preserved.
D5 no hold-out → identity test, 25 lines, PASS.
D6 no user gate → graduation gate; corpus marked operator-corrected vs pending.
D7 scope honesty → minimal pass over all 200 audit rows; depth bounded by class.
