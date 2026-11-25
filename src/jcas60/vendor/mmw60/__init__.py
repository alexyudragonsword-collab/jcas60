"""Four files of mmw60_trx (see ``jcas60.vendor``): the Golay recursion, the standard
single-carrier training field built on it, its own table of assumed / derived numbers, and its
comparison table as printed."""
from pathlib import Path

HERE = Path(__file__).resolve().parent
SURVEY = HERE / "data" / "spec" / "survey.yaml"
COMPARISON_TABLE = HERE / "examples" / "expected" / "04_comparison_table.txt"
