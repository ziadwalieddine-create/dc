"""Variety tests — seeded randomness, weighted sampling, variant expansion."""
import sys, os, shutil, tempfile, json, random
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from dc_extended.variety import derive_seed, weighted_sample, is_synonym
from dc_extended.switchboard import run_turn
from datetime import datetime

SRC = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "cases")
FRI_EVEN = datetime(2026, 9, 18, 18, 11)


def t(name, cond, detail=""):
    print(("PASS " if cond else "FAIL ") + name + (" | " + str(detail) if detail and not cond else ""))
    if not cond:
        sys.exit(1)


def fixture():
    d = tempfile.mkdtemp()
    shutil.copytree(SRC, os.path.join(d, "cases"))
    corpus = os.path.join(d, "sends.jsonl")
    shutil.copy(os.path.join(os.path.dirname(SRC), "corpus", "sends.jsonl"), corpus)
    return d, corpus


# 1. derive_seed is deterministic and varies by input
t("seed deterministic", derive_seed("rae", "2026-09-18T18:11:00") == derive_seed("rae", "2026-09-18T18:11:00"))
t("seed varies", derive_seed("rae", "a") != derive_seed("rae", "b"))

# 2. same seed -> same pick (reproducibility)
items = ["a", "b", "c"]
weights = [3.0, 2.0, 1.0]
r1 = random.Random(42); i1 = weighted_sample(items, weights, r1)
r2 = random.Random(42); i2 = weighted_sample(items, weights, r2)
t("same seed same pick", i1 == i2)

# 3. different seeds can differ
picks = {weighted_sample(items, weights, random.Random(s)) for s in range(50)}
t("variety over seeds", len(picks) > 1, picks)

# 4. higher weight wins more often
wins = {"a": 0, "b": 0, "c": 0}
for s in range(300):
    i = weighted_sample(items, weights, random.Random(s))
    wins[items[i]] += 1
t("weighting respected", wins["a"] > wins["b"] > wins["c"], wins)

# 5. argmax fallback when rng is None
t("rng None -> argmax", weighted_sample(items, weights, None) == 0)
t("single item", weighted_sample(["x"], [1.0], random.Random(1)) == 0)

# 8. integration: run_turn accepts rng_seed, logs seed, output is a valid variant
d, cp = fixture()
r1 = run_turn("rae", os.path.join(d, "cases"), cp, "mean", FRI_EVEN, rng_seed=99)
t("seed in skeleton", r1["skeleton"].get("seed") == 99, r1["skeleton"].get("seed"))
t("valid variant", any(r1["text"] in v for v in [
        ["Come be mean then", "Still waiting on that mean", "Mean. Now."],
        ["Hundreds 😂 Luckily only one was worth finding"],
    ]), r1["text"])

# 9. same seed replays identically
d, cp = fixture()
r2 = run_turn("rae", os.path.join(d, "cases"), cp, "mean", FRI_EVEN, rng_seed=99)
t("replay identical", r2["text"] == r1["text"], (r1["text"], r2["text"]))

# 10. different seed can differ (variants exist)
texts = set()
for s in range(30):
    d, cp = fixture()
    r = run_turn("rae", os.path.join(d, "cases"), cp, "mean", FRI_EVEN, rng_seed=s)
    texts.add(r["text"])
t("output varies by seed", len(texts) > 1, texts)

# is_synonym guard
t("synonym trailing adjective", is_synonym("the basil plant is a lot", "the basil plant is a lot tall"))
t("not synonyms", not is_synonym("basil plant died", "thursday at eight"))

print("ALL VARIETY TESTS PASS")
