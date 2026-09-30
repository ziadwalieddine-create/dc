"""dc_extended.amplitude — heat causally influences the permitted/target band.

Not mechanically monotonic: heat 2 OPENS a stronger band than heat 0, all else
equal; it never FORCES a hotter answer — literal obligations, date coordination,
timing or thread state can still make restraint strategically superior.
logistics=true never imposes an amplitude penalty. Caps are reversible on later
genuine heat (a dry streak must not permanently neuter the thread).
"""
from .licences import Licences


def permitted_band(read, licences=None, thread_cold=False, literal_obligations=None):
    """Returns (floor, ceiling) of permitted response amplitude.
    heat shifts the band; obligations can argue for restraint inside it."""
    heat = (read or {}).get("heat") or 0
    lic = licences or Licences(read or {})
    floor, ceiling = 0, 0
    if heat == 0:
        ceiling = 0
    elif heat == 1:
        ceiling = 1
    else:  # heat 2: stronger band — but see below
        ceiling = 2
    # licence widens the band, never narrows it
    if lic.escalation_licensed and ceiling < 2:
        ceiling = 2
    # literal obligations NEVER reduce the ceiling (non-regression);
    # they only recommend where in the band to aim
    recommended = ceiling
    if literal_obligations and heat >= 1:
        recommended = max(1, min(ceiling, 1 + (1 if lic.escalation_licensed else 0)))
    # thread cold (dry streak) narrows only if no fresh heat evidence
    if thread_cold and heat == 0:
        ceiling = min(ceiling, 0)
    return (floor, ceiling, recommended)


def restore_band(case, read):
    """A later genuine heat event lifts a dry-streak cap (reversibility)."""
    heat = (read or {}).get("heat") or 0
    if heat >= 2:
        case.pop("amp_cap", None)
        return True
    return False
