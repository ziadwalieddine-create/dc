"""dc_extended.licences — semantic licences, not phrase blacklists.

"prove it" is neither globally allowed nor globally forbidden. A move family is
licensed by EVIDENCE in the read/state, never by wording. Paired invariant for
every licence: genuine licence permits the move; the same wording without
licence does not. (Katherine's central lesson: reciprocal challenge vs imposed
audition; qualification licensed only by her frame.)
"""
from .read_v2 import read_v2


class Licences:
    """Evidence-bound permissions derived from the read + case state."""

    def __init__(self, read, case=None):
        r = read or {}
        case = case or {}
        self.r = r
        self.playful = "playful_resistance" == r.get("refusal_scope") or                        "tease" in (r.get("relational_functions") or [])
        # she opened a deserve/test frame (katherine msg 32-style): qualification licensed
        self.qualification_licensed = bool(
            (case.get("open_frames") or {}).get("deserve_game")) or             "Are you sure you deserve" in str(case.get("last_her_challenge") or "")
        # she issued an explicit challenge to his originality/ability
        self.reciprocal_challenge_licensed = bool(
            (case.get("open_frames") or {}).get("her_challenge")) or             any(x in str(case.get("last_her_challenge") or "").lower()
                for x in ("prove", "be more original", "try me", "make me", "you think you"))
        # misread licensed: positive evidence ABOUT HIM; warmth alone is NOT evidence
        rf = r.get("relational_functions") or []
        self.misread_licensed = (r.get("heat") or 0) >= 1 and                                 any(x in rf for x in ("compliment", "tease", "attraction_cue"))
        # escalation licensed: explicit attraction cue or heat 2 with touch
        self.escalation_licensed = "attraction_cue" in (r.get("relational_functions") or [])                                    or r.get("heat") == 2

    def permits(self, move_family):
        """Decision gate. Same wording, different licence -> different verdict."""
        table = {
            "qualification": self.qualification_licensed,
            "challenge":     self.reciprocal_challenge_licensed,
            "misread":       self.misread_licensed,
            "escalation":    self.escalation_licensed,
        }
        return table.get(move_family, True)   # unlisted families: no licence required
