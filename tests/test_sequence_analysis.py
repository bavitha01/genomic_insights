import unittest
from scripts.sequence_analysis import calculate_gc_content

class TestSequenceAnalysis(unittest.TestCase):
    def test_gc_content(self):
        # Mock sequence with 50% GC content
        class MockRecord:
            seq = "ATGC" * 10
        record = MockRecord()
        self.assertEqual(calculate_gc_content(record), 50.0)

if __name__ == "__main__":
    unittest.main()