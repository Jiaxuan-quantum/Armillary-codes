"""Catalog identity and physical replay tests; no distance optimization."""
from collections import Counter
import unittest

from armillary.algebra import indices, legacy_logicals
from armillary.candidates import candidate_definition, list_candidates, load_candidate


class CandidateTests(unittest.TestCase):
    def test_catalog_population(self):
        rows = list_candidates()
        self.assertEqual(len(rows), 122)
        self.assertEqual(len({r["matrix_sha256"] for r in rows}), 122)
        self.assertEqual(Counter(r["exact_distance"] for r in rows),
                         {12: 18, 14: 96, 16: 8})
        for row in rows:
            self.assertEqual(row["n"], 12 * row["q"] * row["p"])
            self.assertEqual(candidate_definition(row["candidate_id"])["matrix_sha256"],
                             row["matrix_sha256"])

    def test_all_reconstructed_matrices(self):
        for row in list_candidates():
            with self.subTest(candidate=row["candidate_id"]):
                code = load_candidate(row["candidate_id"])
                self.assertTrue(code["_verification"]["slots_restored"])
                self.assertTrue(code["_verification"]["forward_reverse_agree"])
                self.assertTrue(code["_verification"]["girth_ge6"])
                self.assertEqual(len(code["_hx"]), code["n"] // 4)
                self.assertEqual(len(code["_hz"]), code["n"] // 4)

    def test_highlighted_aliases(self):
        for alias, expected in {"8749": (3636, 1822, 12),
                                "9053": (4932, 2470, 14),
                                "p181": (8688, 4348, 16)}.items():
            row = candidate_definition(alias)
            self.assertEqual((row["n"], row["k"], row["exact_distance"]), expected)

    def test_logical_basis_convention(self):
        # The four-qubit CSS code with XXXX and ZZZZ stabilizers has k=2.
        logicals = legacy_logicals([0b1111], [0b1111], 4)
        self.assertEqual(logicals, [0b0011, 0b0101])
        self.assertEqual(indices(logicals[1]), [0, 2])
        self.assertTrue(all((row & 0b1111).bit_count() % 2 == 0 for row in logicals))


if __name__ == "__main__":
    unittest.main()
