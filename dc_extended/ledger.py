"""dc_extended.ledger — append-only outcome log; measured family rates."""
import json, os

RANK = {"extends":3,"counters":3,"answers":2,"mirrors_length":1,
        "dry":0,"silence":0,"refuses":-1}

def load_corpus(path):
    if not os.path.exists(path):
        return []
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                rows.append(json.loads(line))
    return rows

def record_send(path, case_id, family, text, ts):
    with open(path, "a") as f:
        f.write(json.dumps({"case":case_id,"family":family,"text":text,
                            "outcome":None,"ts":ts,"frame":False})+"\n")

def record_outcome(path, case_id, text, outcome, ts):
    rows = load_corpus(path)
    hit = False
    for r in reversed(rows):
        if r["case"]==case_id and r["text"]==text and r.get("outcome") is None:
            r["outcome"] = outcome; r["outcome_ts"] = ts; hit = True
            break
    if hit:
        with open(path,"w") as f:
            for r in rows:
                f.write(json.dumps(r)+"\n")

def family_stats(rows):
    stats = {}
    for r in rows:
        fam = r.get("family"); oc = r.get("outcome")
        if not fam or oc is None: continue
        s = stats.setdefault(fam, {"n":0,"score":0.0})
        s["n"] += 1
        s["score"] += max(RANK.get(oc,0),0) / 3.0
    return {k:{"n":v["n"],"rate":round(v["score"]/v["n"],2)} for k,v in stats.items()}
