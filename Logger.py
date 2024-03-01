from babel.numbers import format_decimal
from rich.console import Console
from rich.highlighter import JSONHighlighter
from rich.logging import RichHandler
from rich.theme import Theme
import locale

import json
import logging


class Logger:
    SUBPROC = 25  # Between INFO (20) and WARNING (30)
    HELP = 15  # Between DEBUG (10) and INFO (20)

    def __init__(self):
        self.user_locale = locale.getlocale()[0]  # Get the user's locale
        console = Console(
            stderr=True,
            theme=Theme(
                {
                    "logging.level.subproc": "bold blue",
                    "logging.level.help": "bold green",
                }
            ),
        )

        logging.basicConfig(
            level=logging.NOTSET,
            format="%(message)s",
            datefmt="[%X]",
            handlers=[RichHandler(console=console)],
        )
        self.logger = logging.getLogger(__name__)

        logging.addLevelName(self.SUBPROC, "SUBPROC")
        logging.addLevelName(self.HELP, "HELP")

    def format_numbers(self, message):

        if isinstance(message, str):
            lines = message.splitlines()
            for i, line in enumerate(lines):
                words = line.split()
                for j, word in enumerate(words):
                    try:
                        # Try to convert the word to a float
                        num = float(word)
                        # If successful, format the number and replace the word
                        words[j] = format_decimal(num, locale=self.user_locale)
                    except ValueError:
                        # If the word can't be converted to a float, ignore it
                        pass
                # Join the words back into a single string
                lines[i] = " ".join(words)
            # Join the lines back into a single string
            message = "\n".join(lines)
        elif isinstance(message, int):
            message = format_decimal(message, locale=self.user_locale)
        return message

    def info(self, message):
        message = self.format_numbers(message)
        self.logger.info(message)

    def debug(self, message):
        message = self.format_numbers(message)
        self.logger.debug(message)

    def warn(self, message):
        message = self.format_numbers(message)
        self.logger.warning(message)

    def error(self, message):
        message = self.format_numbers(message)
        self.logger.error(message)

    def subproc(self, message, *args, **kwargs):
        message = self.format_numbers(message)
        if not message:  # Check if the message is empty
            message = "No errors reported"
        if self.logger.isEnabledFor(self.SUBPROC):
            self.logger._log(self.SUBPROC, message, args, **kwargs)

    def help(self, message, *args, **kwargs):
        message = self.format_numbers(message)
        if not message:
            message = "No help available"
        if self.logger.isEnabledFor(self.HELP):
            self.logger._log(self.HELP, message, args, **kwargs)

    def json(self, data):
        json_str = json.dumps(data, indent=4)
        self.logger.info(json_str, extra={"highlighter": JSONHighlighter()})
