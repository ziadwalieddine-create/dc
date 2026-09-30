# Live runbook — Phase 5

The framework is only real when the outcomes are. Per turn:

1. `run_turn(...)` — parse, read, gate, retrieve, emit. If `action: "deep"`,
   resolve the split (model judgment), then re-run `force_fast=True`.
2. Operator sends the line, then: `confirm_send(cases_dir, case_id)` — the
   receipt. Nothing counts before this (dc: actually_sent only from him).
3. She replies. Class it (extends / answers / counters / mirrors_length /
   literal / dry / topic_change / silence / refuses), then:
   `her_reply(cases_dir, case_id, bubble, response_class)` — one call updates
   reciprocity, grades any pending split read, and records the frame outcome.
4. Nightly: `nightly_checkin(cases_dir)` — queue + honesty ratio + proposed
   gates + gate hit rates. Two incidents of a family propose a gate; the human
   admits it.

Growth rule: new frames enter only where the ledger shows use or failure.
Shrinkage (k=5) is live from the first outcome — five rows never crown a frame.
