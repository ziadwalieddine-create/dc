# Skeleton — the full-output contract

Every `run_turn` populates this skeleton before emitting. Code fills every field.
The model does not invent values. Unresolved stays unresolved.

```
CASE:     {name} | {stage} | {channel} | goal={goal.type}
READ:     {kind} | heat={0|1|2} | she_extended={bool} | confidence={high|split}
          reversal: {the one fact that would kill this read}
CLAIM:    {type} | family={family} | builds_on_last={bool}
SEED:     {int | null}              — reproducibility seed (variety layer)
TRIAGE:   {receptive | neutral | cooling | unreceptive | closed | null}
ENDPOINT: {off_app | date_lock | banter | close | null}
FRAME:    {id, variant, card}       — certified frame reference + justification
GATES:    {gate_id: PASS|FAIL — detail}
CANDIDATES:
  1. {text} | source={retrieved|edited} | family_rate={0.00} | heat_fit={0|1|2} | amp_fit={0|1|2}
     critic: {kill_id | PASS}
     simulate: {extends|answers|dry|silence|topic_change}
  2. ...
  3. ...
SELECT:   {winner_index} | because {reason}
SURFACE:  protected_tokens={PASS|FAIL} | reverify={PASS|FAIL}
RELEASE:  {final text}
```

## Rules

- Every candidate that survives the critic gets a simulation.
- The winner is the candidate with the best simulation outcome, then heat fit,
  then amplitude fit, then family rate.
- If no candidate survives the critic, Resolve: name the kill, emit no line.
- The skeleton is stored on the case chronology, not printed to her.
