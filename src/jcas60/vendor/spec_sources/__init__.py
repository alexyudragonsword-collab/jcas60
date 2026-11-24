"""spec_sources' reader, re-exported (the package ``__init__`` files in the
vendor tree are jcas60's own; ``reader.py`` and ``data/`` are the copies)."""
from .reader import (SOURCES, DeprecatedKeyWarning, SpecNumber, all_entries, audit_tags, by_source, data_files,
                     load_table, merge_tables, ref, survey, value)

__all__ = ["SOURCES", "DeprecatedKeyWarning", "SpecNumber", "all_entries", "audit_tags", "by_source", "data_files",
           "load_table", "merge_tables", "ref", "survey", "value"]
