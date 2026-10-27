"""Regression checks for scientific comparison failures, not algorithm mirrors."""
import unittest

import pandas as pd

from recompute import compare_frames


class ScientificComparisonTests(unittest.TestCase):
    def setUp(self):
        self.expected = pd.DataFrame({"species_id": [1, 2], "area": [100.0, 200.0], "p_value": [0.001, 0.02]})

    def test_missing_result_column_fails(self):
        report = compare_frames(self.expected.drop(columns="p_value"), self.expected.copy(), "missing", ["species_id"])
        self.assertEqual(report["status"], "FAIL")
        self.assertEqual(report["missing_columns"], ["p_value"])

    def test_changed_species_set_fails_even_when_values_match(self):
        actual = self.expected.copy()
        actual["species_id"] = [1, 3]
        self.assertEqual(compare_frames(actual, self.expected.copy(), "identifiers", ["species_id"])["status"], "FAIL")

    def test_changed_result_fails_and_reordered_equal_rows_pass(self):
        actual = self.expected.iloc[::-1].copy()
        self.assertEqual(compare_frames(actual.copy(), self.expected.copy(), "reordered", ["species_id"])["status"], "PASS")
        actual.loc[actual["species_id"].eq(1), "p_value"] = 0.5
        self.assertEqual(compare_frames(actual, self.expected.copy(), "wrong", ["species_id"])["status"], "FAIL")


if __name__ == "__main__":
    unittest.main()
