"""A small syndrome-only MLE example; no code-sized solve is required."""
import argparse

import numpy as np

from armillary.decoders import MLEDecoder

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--backend", choices=("scipy", "copt"), default="scipy")
args = parser.parse_args()

H = np.array([[1, 1, 0], [0, 1, 1]], dtype=np.uint8)
priors = np.array([.1, .1, .1])
syndrome = np.array([1, 1], dtype=np.uint8)
answer = MLEDecoder(H, priors, backend=args.backend).decode(syndrome)
print("Syndrome:", syndrome.tolist())
print("Correction:", answer["correction"].tolist() if answer["accepted"] else None)
print("Diagnostics:", answer["diagnostics"])
