# Critic

First hit kills the string. The critic returns the id and no replacement. The other candidate is the only replacement. If none remain, Resolve, and name the kill.

1. `not_referable` — the string is empty.
2. `claim_swap` — the string books or names a time while `claim.type` is not `plan`, or changes channel while `claim.type` is not `channel`.
3. `clock` — a weekday token in the string is dead at the live clock (`case["_now"]`).
4. `not_particular` — the line is a member of GENERIC_LINES.
5. `label_question` — dream, fear, childhood, job, origin, type, perfect day, crystal ball.
6. `interview` — more than one question.
7. `self_decode` — the joke is explained. "that was a joke", "just kidding".
8. `aimed_at_her` — the butt is her body, her intelligence, her friends, or women as a class.
9. `missed_heat` — amplitude 0 in Play/Trade while `she_extended` (run_turn sets `_policy`/`_amplitude`).
10. `negotiates_no` — the claim's family is in `case["banned"]`.
11. `dead_end` — simulation classes are only silence and topic_change, and another candidate is wider (run_turn passes `sim_classes`).
12. `length` — over 220 characters when her length is null or short.
13. `echo` — the exact string is already banned, or the specimen sentence: "You kept the basil and threw the recipe. That tracks."
14. `unnatural` — "not x but y", or a line that starts with "I hope this helps".
15. `logistics_mush` — "when are you free", "your call", "pick a day", "you pick".
16. `frame_limit` — "fwb" or "situationship" under `no_arrangement_label`, or "our relationship" under `no_relationship_speech`.
17. `known_fact` — the line is a question that already contains a FACT longer than eight characters.
18. `restart` — channel is ig or sms and the line says "nice to match", "found you on hinge", "we matched", or "hey from hinge".

Every kill reads real state or the string itself. There is no 1-10 score and no percent. Neediness is `binds_her`. A flat line on an extended flirt is `missed_heat`, not caution.

Removed in the kill audit: `provenance` (superseded by frame certification's `requires`), `sand` (no tease-removal operation exists in a frame pipeline), `cold_read` (read.md's correctable-read allowance governs), `no_affordance` (certification requires declared affordance slots). Flag-only paths `joke_only`, `generic`, `must_be_fact`, `tease_removed`, `continues_banned`, `cold_read_as_fact` are gone; the live halves of `not_referable` (empty text) and `not_particular` (GENERIC_LINES) remain.
