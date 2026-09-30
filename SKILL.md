---
name: dc-extended
description: |
  A state machine for dating-app threads with executable gates, a structured
  read layer, a first-hit critic, a heat-aware retrieval writer, a nightly
  queue, and an append-only ledger. The model reads her bubble and writes
  the sentence. The code owns everything else.
---

# dc-extended

A state machine with executable gates, a structured read layer, a heat-aware
retrieval writer, a first-hit critic, a nightly queue, and an append-only ledger.

## The loop

1. **Check in nightly.** Build the queue. Every live thread has one next action
   with a deadline, or a coded close.
2. **One turn at a time.** Call `run_turn` — load, read, classify, write, gate,
   critic, write-forward, emit.
3. **Read first.** `dc_extended/read.md` — one dominant read (refuse / dry /
   escalate / tease / offer / logistics / question / withdraw / topic / invest /
   none) with heat, she_extended, confidence, and reversal. Stored on the case,
   not printed.
4. **Confirm sends.** `sent_confirmed` is set only when the operator says it
   left the phone.
5. **Record outcomes.** Her response class lands in the ledger against the
   send's family.

## The twelve gates

Nine gates run in `dc_extended/gates.py`:

wrong_person · wrong_platform · cross_thread · dead_clock · restack ·
dual_plan · stale_goal · unconfirmed_send · stated_boundary

Three gates run in `dc_extended/switchboard.py` before or around the gate layer:

- `paste_without_send` — an audit claim stops the turn before any gate runs.
- `refusal_close` — a refuse read closes the thread before any gate runs.
- `staleness` — a stale case (>48h) stops the turn before any gate runs.

`no_write` is reserved for a host that cannot persist state. It never fires
in this runtime because `save_case` is atomic and always called.

Rules of the gate layer:
- Gates are code. Everything else is craft.
- A new gate enters only after two logged incidents of a new family.
- Refusal kills a family permanently. Dry outcomes go to the rate table,
  never to bans.
- A missed window is dead. Never propose the adjacent day as a synonym.

## The read layer

`dc_extended/read.md` — eleven kinds, one dominant read per turn. Heat is 0, 1,
or 2. Not a maybe. A refusal stops the action before any claim work. Those are
evidence, not nerves.

## The writer

Retrieval first: proven lines from the operator's own corpus, ranked by outcome
class, family rate, heat fit, and amplitude fit. Parameterized with live case
facts. Screened through `dc_extended/critic.md` — first hit kills, next candidate
is the replacement. Fallback: minimal edit of his draft. Never cold invention —
no corpus match and no draft returns `DATA_MISSING`, not a guess.

## The variety layer

`dc_extended/variety.py` — seeded randomness, in selection only. Gates, reads,
and critic kills stay deterministic. Randomness picks among critic-surviving
candidates by weighted sampling (higher-scored candidates win more often, not
always), and corpus frames may carry a list of phrasings instead of one string.

Every seed is derived from case + turn and logged in the skeleton. Same seed +
same input = same output. Any send can be replayed exactly.

`run_turn(..., rng_seed=None, temp=1.0)` — pass `rng_seed` to force a pick,
`temp` to flatten (higher) or sharpen (lower) the weighting.

## The critic

`dc_extended/critic.md` — 23 named kills. First hit kills the string. No
replacement from the critic. The other draft is the only replacement. If none
remain, Resolve, and name the kill.

## The skeleton

`dc_extended/skeleton.md` — the full-output contract. Every `run_turn` populates
it from actual pipeline state before emitting. Code fills every field. The model
does not invent values.

```
CASE:     {name} | {stage} | {channel} | goal={goal.type}
READ:     {kind} | heat={0|1|2} | she_extended={bool} | confidence={high|split}
          reversal: {the one fact that would kill this read}
CLAIM:    {type} | family={family} | builds_on_last={bool}
GATES:    {gate_id: PASS|FAIL — detail}
CANDIDATES:
  1. {text} | source={retrieved|edited} | family_rate={0.00} | heat_fit={0|1|2}
     critic: {kill_id | PASS}
     simulate: {extends|answers|dry|silence|topic_change}
  2. ...
  3. ...
SELECT:   {winner_index} | because {reason}
SURFACE:  protected_tokens={PASS|FAIL} | reverify={PASS|FAIL}
RELEASE:  {final text}
```

The skeleton is returned in every `run_turn` result and stored on the case
chronology. It is not printed to her.

## Non-negotiables

- Emit and write-forward are one transaction. A run that emits without writing
  is a failed run.
- Analysis is not a send. Audit requests return `paste_without_send`.
- Closed threads generate nothing. Closure is a decision with a cause, not a drift.

## Output contract

Success: `{text, claim, family, source, read, heat}` — one copyable string.
Stop: `{gate | why, detail}` — one line, no replacement.

## Cases

One JSON file per thread (`cases/*.json`). Schema in `cases/_template.json`.
Ported live cases: rae, mary, jisu, camille. Write the case forward every turn
or the file is lying.

## Running a turn

```python
from dc_extended.switchboard import run_turn

result = run_turn(
    case_id="rae",           # the thread to act on
    cases_dir="cases/",      # directory of case JSON files
    corpus_path="corpus/sends.jsonl",  # the operator's send corpus
    operator_text="mean",    # what the operator wants to do
    now=None,                # datetime; defaults to datetime.now()
    her_bubble=None,         # her latest message; defaults to case's last
    operator_draft=None,     # his draft text; used as writer fallback
)
# result: {"action": "send", "text": ..., ...} or {"action": "stop", ...}
```

## Evals

`python3 evals/test_gates.py` — gate units.
`python3 evals/test_queue.py` — window/decay.
`python3 evals/test_read.py` — read layer kinds, heat, reversal.
`python3 evals/test_critic.py` — all 23 kill ids + clean-line pass.
`python3 evals/test_writer.py` — retrieval ranking, heat fit, critic screening,
fallback.
`python3 evals/replay_incidents.py` — the acceptance suite: every named loss
from the prior system's ledger, replayed. Green means none can recur.

## Phase 0 additions

- `dc_extended/intake.py` — `open_turn` parses raw trigger text (`she:`/`her last:`/
  `me:` marks) into bubble + claim; precedence by span-aware last-match, so
  "don't book thursday" locks `bit`. `bind_case` stops duplicate case files.
- `dc_extended/recovery.py` — `record_deviation(case, event, mechanism, family, gate)`
  is the single recovery path (rejection transaction + transitions + incident
  intake). Mechanism stays UNKNOWN unless named. `forced_resolution` degrades
  non-integrity holds to amplitude-0 Play/Trade; integrity stops never degrade.
  `update_reciprocity` is the per-case thread-level indicator (receipt-gated,
  clamped, two dries cap amplitude at 0).
- `variety.is_synonym` — guards the two-draft rule against synonym re-rolls.
- Case schema gains `endpoint`, `reciprocity`, `amp_cap`, `deviations`.

Acceptance: 9/9 suites green (8 existing + test_phase0, 41 checks).

## Phase 1 additions

- `dc_extended/readlayer.py` — the read layer, upgraded:
  - `read_bubble` (base read) moved here from switchboard; switchboard re-exports.
  - `read_bubble_v2` — real confidence: `split` when one bubble supports two
    kinds that select different actions (charge vs withdraw, tease vs withdraw,
    question vs withdraw). Split picks the challenge (read.md).
  - `push_window` — last 3 of her bubbles on the case; window patterns
    (`testing_investment`, `push_pull`, `fading`) annotate low-energy reads.
  - Compute contract: `run_turn(..., force_fast=False)` returns
    `{"action": "deep", "why": "split_read"}` with no line. Deep mode writes
    state (`pending_read` in the skeleton), never essays. `resolve_split(case,
    chosen)` applies the judgment; re-run with `force_fast=True`.
  - `grade_pending_read(case, response_class)` — her confirmed reply grades the
    split (extends/counters -> primary; dry/silence -> secondary; answers/
    literal/refuses -> no count). Accumulates into `read_stats.json`: the read
    layer becomes a scored classifier from data that arrives anyway.

Acceptance: 10/10 suites green (9 existing + test_phase1, 30 checks).

## Phase 2 additions

- `dc_extended/goals.py` —
  - `lock_endpoint(text, stored)` — the operator's goal THIS TURN wins
    (kernel O\*; `substituted_end` guard). Never ladder-inferred. Stored with
    source, like claim. Unknown text -> None (UNRESOLVED legal).
  - `sim_rank(endpoint)` — endpoint-weighted outcome ranking: the same reply
    is worth different amounts against different endpoints (banter boosts
    extends, date_lock boosts counters, off_app boosts answers).
  - `triage(case)` — receptive / neutral / cooling / unreceptive, derived from
    reciprocity + recent outcomes. Computed per turn, never stored.
    cooling = his one-move budget is spent; hold until she re-invests.
  - `queue_weight(case)` — triage as queue priority (receptive first).
- `run_turn(..., endpoint=None)` — endpoint kwarg or parsed from operator text;
  skeleton carries `endpoint` + `triage`.
- Cooling gate: `triage == "cooling"` and no extension -> Hold (no line),
  unless `force_fast`. Her extension lifts it.
- `build_queue` sorts by (deadline, triage weight).

Acceptance: 11/11 suites green (10 existing + test_phase2, 30 checks).

## Phase 3 additions

- `dc_extended/lexicon.py` — ONE wordlist for read layer + critic + claim map.
  Identity (`is`) is the drift-proof guarantee; test_phase3 asserts it.
- Corpus rows are now FRAMES: `frame: true`, `id`, `slots`, `forecast`
  (build-time-certified response classes), `protected` tokens, `endpoint`.
  Harvested sends carry `frame: false` and skip certification.
- `dc_extended/frames.py` — the certification oracle: `certify_frame` /
  `certify_corpus`. Checks: required fields, critic-pass under fixtures,
  particularity (referent-swapped fills must differ), filler-lexicon lint,
  forecast validity. Certification IS the test suite — continuous, zero
  tooling; a stricter critic re-fails stale frames loudly.
- Behavior seal un-deferred: every send emits `frame: {id, variant}` + the
  logged seed. Any send replays exactly from (case state, seed, frame_id).
  The specimen-banned basil variant was replaced; the frame now certifies.

Acceptance: 12/12 suites green (11 existing + test_phase3, 17 checks).

## Phase 5 additions — live data

- `dc_extended/live.py` — the live-thread harness:
  - `confirm_send(cases_dir, case_id, text=None)` — receipt discipline; the
    confirmed text must be the candidate exactly.
  - `her_reply(..., bubble, response_class)` — one confirmed reply does three
    jobs: reciprocity, split-read grading (archived by resolve_split), and
    frame outcome attribution via `last_send.frame_id`.
  - `record_frame_outcome` / `shrunk_rate` — frame statistics with shrinkage
    (k=5): n < k leans on the family prior; five rows never crown a frame.
- `run_turn(..., frame_stats_path=None)` — candidate scoring uses shrunk frame
  rates when stats exist.
- `LIVE_RUNBOOK.md` — the per-turn protocol for real threads.

Acceptance: 14/14 suites green (13 existing + test_phase5, 21 checks).
Phase 4 closed: gate-stop deviations record the state slice (banned_families,
open_plan); replay applies it verbatim; engine-recorded rows are
self-consistent by construction.

## Fix plan execution (audit remediation)

Two full audits (wiring + deep logic) found: unwired components, dead critic
kills, fixture-set fields production never computes, and read-layer false
positives. All seven fixes are executed and regression-locked in
`evals/test_fixes.py`:

1. **Read layer** — `don't`/`do not` removed from refusal tokens ("don't be
   silly" no longer closes a thread); charge/withdraw/classify keywords are
   word-boundary matched ("I know right" no longer escalates); dead window
   pattern removed.
2. **Clock layer** — `materialize_plan` stores the day-word fact;
   `derive_plan_window` recomputes the clock every turn (never stored — the
   Jisu staleness failure in miniature). Claim-level dead_clock gate; critic
   clock kill scans candidate text against the live clock. Refusal-close sets
   `closed` (the abolished `stage` field is gone from schema and fixtures).
3. **Seams** — `trigger_text` wired to `open_turn` at run_turn; intake claims
   preserve `family` from classify; `amp_cap` excludes over-cap frames;
   `forced_resolution` degrades dead_clock holds only (unconfirmed_send always
   stops — no second bid); `is_synonym` drops pool variants re-rolling killed
   text; `writer.draft`/`pick_variant`/`jitter` deleted.
4. **Gate honesty** — `run_gates` emits a real per-gate PASS/FAIL log (the
   fabricated all-PASS block is gone); unconfirmed_send = unconfirmed-no-bubble
   OR builds_on_last union; no_write removed; goal_smuggle reads the
   turn-locked endpoint first.
5. **Live integrity** — `her_reply` is idempotent (duplicate replies return
   `duplicate: True`, keyed on bubble+class+send ts); `charter_query` no longer
   re-proposes banned families.
6. **Hygiene** — `save_case` strips transient underscore keys; fixture schema
   validated by `test_schema.py`; FLEET.md refreshed.
7. **Critic kill audit** — 23 -> 19 kills. Deleted: provenance, sand,
   cold_read, no_affordance (all test-flag-only). Re-wired from real state:
   negotiates_no (banned families), missed_heat (_policy/_amplitude),
   dead_end (sim_classes from run_turn). critic.md rewritten; corpus
   re-certifies.

16/16 suites green including `test_fixes.py` and `test_schema.py`.

## Fix Plan v3 — seal the boundary, restore the signal

Post-mortem of the first two fix plans found the logic was right but the data
flow across context boundaries was unsealed. All eight fixes executed:

1. **Queue per-stage decay** — `derive_stage(case)` derives urgency from live
   fields (open_plan, channel, last_send, closed). plan_pending=48h, off_app=144h,
   matched=72h, responded=96h. The dead `stage` field and `DM_DEADLINE_HOURS`
   branch are gone. The operator's nightly view now reflects relationship depth.
2. **Stats write path sealed** — `her_reply` and `grade_pending_read` raise
   on `stats_path=None`. Production writes require explicit paths. Tests pass
   temp paths.
3. **frame_stats.json reset** — fake test data (n=10, extends=10) removed.
   Shrinkage (k=5) leans on family priors until real outcomes arrive.
4. **Boundary regression test** — `evals/test_stats_boundary.py` redirects
   stats paths to temp files, runs a full turn, asserts production files are
   byte-identical. The seal is permanent.
5. **simulate.py corpus token removed** — "the basil" no longer leaks into
   affordance detection.
6. **skeleton.md synced** — seed, triage, endpoint, frame documented.
7. **read.md** — "No apology, no restart" replaces the token-overlapping phrase.
8. **ALL_GATES split** — GATE_LAYER_GATES (10) mirrors run_gates exactly.
   SWITCHBOARD_STOPS (3) documented separately. No phantom traps.

17/17 suites green. frame_stats.json is empty and sealed.

## Flirt model — generation layer

`dc_extended/flirt_model.md` is the generation constitution, distilled from
message-level analysis of live threads (Gaby thread, dual-analysed) plus the
provenance audit (152 sent / 42 rejected lines). Training data:
`training/flirt_positives.jsonl` (24 exemplars, provenance-graded per
authorship controls — grades present: A, A/B, B; no R-graded entries included,
reconstructed wording excluded from the positive corpus — every exemplar shipped
with the cue that produced it) and
`training/flirt_negatives.jsonl` (25 exemplars, every one with its documented
rejection reason).

Rules of use:
- Detached lines are NOT training-equivalent. Never imitate a line without
  its cue.
- The style model is literal uptake → relational tilt → low-friction return.
- Run the 10-step generation algorithm before producing any line.
- Score candidates on the 6-dimension rubric (>=9/12, no hard rejection).
- The critic kills stop bad lines; the flirt model shapes good ones. Both run.

## Flirt training pipeline v2

`training/PLAN.md` governs. Teaching mechanism: preference pairs
(`training/flirt_pairs.jsonl`). Two tracks (attraction/execution), move-unit
annotation, procedure traces (`training/procedure_traces.jsonl`), minimal pass
over all 200 audit rows (`training/audit_minimal_pass.jsonl`), hold-out identity
test (`evals/test_identity_holdout.py` — PASS), graduation gate on the operator
before any annotation enters the corpus. Open fork: Katherine files.
