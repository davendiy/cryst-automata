
from dataclasses import dataclass
from typing import Generic, TypeVar


T = TypeVar('T')


class _noop:

    def __getattribute__(self, name):
        return self

    def __call__(self, *args, **kwargs):
        return self


@dataclass
class Option(Generic[T]):
    failed: bool
    error_msg: str
    result: T

    def wrapped(self):
        if self.failed:
            return _noop()
        else:
            return self.result

    @classmethod
    def error(cls, error_msg):
        return cls(True, error_msg, None)

    @classmethod
    def success(cls, result: T):
        return cls(False, '', result)
