## ardal/utils/logger.py
import logging

LOGGER_NAME = "ardal"

def get_logger( submodule=None
                )-> logging.Logger:
    name = "ardal" if submodule is None else f"ardal.{submodule}"
    logger = logging.getLogger(name)

    if not logger.handlers:  ## prevent duplicate handlers on reload
        handler = logging.StreamHandler()
        formatter = logging.Formatter(
            fmt="[{asctime}] {levelname}: {message}", style="{"
        )
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        logger.setLevel(logging.WARNING)  ## set default level

    return logger
