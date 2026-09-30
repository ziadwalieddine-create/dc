# Fleet

dc-extended runs a dating turn through a single release authority.

Architecture of record:
- 01-date-controller (dc_extended/controller.py) owns final selection/release.
- 02-date-state (dc_extended/events.py) owns factual truth — immutable, append-only.
- The legacy run_turn remains available as a reference/replay path during migration.

A second callable sender without the controller's authority is the failure this
file exists to prevent.

Hinge profile text is not this turn. The romance persona is not this turn.
