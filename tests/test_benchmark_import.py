import unittest


class BenchmarkImportTests(unittest.TestCase):
    def test_benchmark_module_imports_with_metadata_schema_changes(self):
        import funmirbench.benchmark  # noqa: F401


if __name__ == "__main__":
    unittest.main()
