"""Centralized logging configuration for fun-miRBenvh."""

from __future__ import annotations

import logging
import sys
from typing import Final

DEFAULT_LOG_FORMAT: Final[str] = (
    "%(asctime)s | %(levelname)s | %(name)s | %(message)s"
)
DEFAULT_DATE_FORMAT: Final[str] = "%Y-%m-%d %H:%M:%S"
DEFAULT_LOG_LEVEL: Final[int] = logging.INFO


def parse_log_level(value: str | None) -> int:
    """Parse a user-provided log level string into a logging module constant."""
    if value is None:
        return DEFAULT_LOG_LEVEL

    normalized = str(value).strip().upper()
    mapping = {
        "CRITICAL": logging.CRITICAL,
        "ERROR": logging.ERROR,
        "WARNING": logging.WARNING,
        "INFO": logging.INFO,
        "DEBUG": logging.DEBUG,
    }
    if normalized not in mapping:
        valid = ", ".join(mapping)
        raise ValueError(f"Invalid log level {value!r}. Expected one of: {valid}")
    return mapping[normalized]


def setup_logging(level: int | str = DEFAULT_LOG_LEVEL) -> None:
    """Configure root logging once for the whole application.

    This should be called near the start of each CLI entry point.
    """
    if isinstance(level, str):
        level = parse_log_level(level)

    root_logger = logging.getLogger()
    root_logger.setLevel(level)

    # Avoid duplicated handlers if setup_logging() is called multiple times.
    if root_logger.handlers:
        root_logger.handlers.clear()

    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(level)
    handler.setFormatter(
        logging.Formatter(
            fmt=DEFAULT_LOG_FORMAT,
            datefmt=DEFAULT_DATE_FORMAT,
        )
    )

    root_logger.addHandler(handler)


def get_logger(name: str) -> logging.Logger:
    """Small helper for consistent logger creation."""
    return logging.getLogger(name)