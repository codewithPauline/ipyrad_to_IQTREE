#!/usr/bin/env python3
"""Generate fictional DNA sequences for a software demonstration; no research data."""
import random
from pathlib import Path
rng = random.Random(2026)
ancestor = [rng.choice("ACGT") for _ in range(1200)]
def mutate(sequence, probability):
    return [rng.choice([x for x in "ACGT" if x != base]) if rng.random() < probability else base
            for base in sequence]
rows = []
for group in "ABCD":
    parent = mutate(ancestor, 0.07)
    for index in (1, 2):
        rows.append((f"Demo_{group}{index}", ''.join(mutate(parent, 0.025))))
path = Path(__file__).with_name("demo.phy")
path.write_text("8 1200\n" + ''.join(f"{name}  {seq}\n" for name, seq in rows))
print(path)
