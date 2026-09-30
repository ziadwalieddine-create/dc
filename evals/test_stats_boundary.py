"""Fix 4: no eval writes to production stats. The boundary is sealed."""
import sys, os, json, tempfile, shutil
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import dc_extended.live as live_mod
import dc_extended.readlayer as rl_mod
from dc_extended.switchboard import run_turn
from dc_extended.live import confirm_send, her_reply
from datetime import datetime

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROD_FRAME_STATS = os.path.join(BASE, "frame_stats.json")
PROD_READ_STATS = os.path.join(BASE, "read_stats.json")
NOW = datetime(2026, 9, 18, 18, 11)


def t(name, cond, detail=""):
    print(("PASS " if cond else "FAIL ") + name + (" | " + str(detail) if detail and not cond else ""))
    if not cond:
        sys.exit(1)


def snapshot(path):
    try:
        return open(path, "rb").read()
    except FileNotFoundError:
        return None


# Snapshot production stats
fs_before = snapshot(PROD_FRAME_STATS)
rs_before = snapshot(PROD_READ_STATS)

# Redirect the module defaults to temp files
d = tempfile.mkdtemp()
live_mod._STATS = os.path.join(d, "frame_stats.json")
rl_mod._STATS_PATH = os.path.join(d, "read_stats.json")

shutil.copytree(os.path.join(BASE, "cases"), os.path.join(d, "cases"))
cp = os.path.join(d, "sends.jsonl")
shutil.copy(os.path.join(BASE, "corpus", "sends.jsonl"), cp)

# Run a full turn + confirm + reply — all writes go to the redirected temp paths
res = run_turn("rae", os.path.join(d, "cases"), cp, "mean", NOW,
               her_bubble="the mean prompt is unhinged", rng_seed=7)
t("turn sends", res["action"] == "send")
confirm_send(os.path.join(d, "cases"), "rae")
out = her_reply(os.path.join(d, "cases"), "rae", "haha yeah you wish", "extends",
                stats_path=os.path.join(d, "frame_stats.json"),
                read_stats_path=os.path.join(d, "read_stats.json"))
t("reply recorded", out["frame"] is not None, out)

# Production stats must be untouched
fs_after = snapshot(PROD_FRAME_STATS)
rs_after = snapshot(PROD_READ_STATS)
t("frame_stats.json untouched", fs_after == fs_before,
  "changed!" if fs_after != fs_before else "")
t("read_stats.json untouched", rs_after == rs_before,
  "changed!" if rs_after != rs_before else "")

# The guard fires on None
try:
    her_reply(os.path.join(d, "cases"), "rae", "x", "extends")
    t("guard fires on None stats_path", False)
except ValueError:
    t("guard fires on None stats_path", True)

from dc_extended.readlayer import grade_pending_read
try:
    grade_pending_read({}, "extends")
    t("guard fires on None stats_path (read)", False)
except ValueError:
    t("guard fires on None stats_path (read)", True)

# Temp stats got the writes
temp_fs = json.load(open(os.path.join(d, "frame_stats.json")))
t("temp frame stats populated", len(temp_fs) > 0, list(temp_fs.keys()))

print("ALL STATS BOUNDARY TESTS PASS")
