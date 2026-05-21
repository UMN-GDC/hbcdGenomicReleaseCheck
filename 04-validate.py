#!/usr/bin/env python3
"""Run the release-data validation tests and log the output.

Python equivalent of 02-validate.R — runs the test suite and captures
results to a log file for review.

Usage
-----
python 04-validate.py [--verbose]
"""
import argparse
import sys


def main():
    parser = argparse.ArgumentParser(
        description="Validate HBCD genomics release data"
    )
    parser.add_argument(
        "--verbose", "-v", action="store_true",
        help="Print test output to stdout as well as the log file"
    )
    args = parser.parse_args()

    import pytest

    test_file = "tests/test_release_data.py"
    log_file = "test-log.txt"

    # Build pytest arguments
    pyargs = ["-x", test_file]  # stop on first failure
    if args.verbose:
        pyargs.insert(0, "-v")

    print(f"Running: pytest {' '.join(pyargs)}")
    print(f"Logging to: {log_file}")

    # Run and capture output
    with open(log_file, "w") as f:
        exit_code = pytest.main(
            pyargs,
            plugins=[],
        )

    # Re-run with verbose to log the full output
    with open(log_file, "w") as f:
        exit_code = pytest.main(
            ["-v", test_file],
            plugins=[],
        )

    if exit_code == 0:
        print("\nAll tests passed.")
    else:
        print(f"\n{exit_code} test(s) failed — see {log_file} for details.")

    sys.exit(exit_code)


if __name__ == "__main__":
    main()
