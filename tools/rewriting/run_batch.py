"""Run prover.py on several systems in parallel with sanitized log names.
Usage: .venv/bin/python run_batch.py TAG "SYS1 SYS2 ..." -- [prover args]"""
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
tag, systems = sys.argv[1], sys.argv[2].split()
rest = sys.argv[sys.argv.index("--") + 1:] if "--" in sys.argv else []
par = int(os.environ.get("PAR", "4"))
os.makedirs(os.path.join(HERE, "logs"), exist_ok=True)


def safe(n):
    return n.replace(">", "E").replace("<", "B").replace("*", "S")


def run(name):
    fn = os.path.join(HERE, "logs", f"{tag}_{safe(name)}.log")
    with open(fn, "w") as fh:
        subprocess.run([os.path.join(HERE, ".venv/bin/python"), os.path.join(HERE, "prover.py"),
                        name] + rest, stdout=fh, stderr=subprocess.STDOUT, cwd=HERE)
    last = [l for l in open(fn) if "PROVED" in l or "result:" in l or "first step" in l
            or "FAILED" in l]
    return name, " | ".join(l.strip() for l in last)


with ThreadPoolExecutor(par) as ex:
    for name, res in ex.map(run, systems):
        print(f"{name:20s} {res}", flush=True)
