"""
utils.py - Common utilities (logging, UUID, memory tracking, date format).
"""
import logging
from logging.handlers import QueueHandler, QueueListener
import uuid
import psutil
import os
from datetime import datetime
from pathlib import Path
from typing import Optional

def setup_logging(log_file: Optional[str], verbose: bool = False) -> logging.Logger:
    """Configure root logger to write to pipeline.log and console."""
    logger = logging.getLogger("ligprep")
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")

    # File handler
    if log_file:
        fh = logging.FileHandler(log_file, mode="w")
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    #console handler (Info and above)
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    return logger

def start_queue_listener(logger: logging.Logger, log_queue) -> QueueListener:
    """
    Start a QueueListener that forwards log records placed on log_queue to the
    main process's real handlers (file + console).

    multiprocessing.Pool forks worker processes that would otherwise inherit
    the already-open FileHandler/StreamHandler, causing multiple processes to
    write to (and close) the same underlying file descriptor concurrently.
    Routing worker log records through a queue keeps all real I/O in the main
    process only.
    """
    listener = QueueListener(log_queue, *logger.handlers, respect_handler_level=True)
    listener.start()
    return listener

def worker_logging_init(log_queue, verbose: bool = False) -> None:
    """
    Pool(initializer=...) target: reconfigure the 'ligprep' logger inside each
    worker process to send records through log_queue instead of the handlers
    it inherited via fork.
    """
    worker_logger = logging.getLogger("ligprep")
    worker_logger.handlers.clear()
    worker_logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    worker_logger.addHandler(QueueHandler(log_queue))
    worker_logger.propagate = False

def unique_id() -> str:
    """return a UUID4 string"""
    return str(uuid.uuid4())

def current_date_tag() -> str:
    """Return date in DDMMYY format, eg. '01AUG26'. """
    return datetime.now().strftime("%d%b%y").upper()

def memory_usage_mb() -> float:
    """Return current process memory usage in MB."""
    process = psutil.Process(os.getpid())
    return process.memory_info().rss / (1024 * 1024)

def log_memory(logger: logging.Logger, tag: str = ""):
    """Log current memory usage with an optional tag."""
    mem = memory_usage_mb()
    logger.info(f"Memory usage {tag}: {mem:.2f} MB")