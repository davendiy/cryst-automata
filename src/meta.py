
from dataclasses import dataclass
from typing import Generic, TypeVar


T = TypeVar('T')


@dataclass
class Option(Generic[T]):
    failed: bool
    error_msg: str
    result: T

    @classmethod
    def error(cls, error_msg):
        return cls(True, error_msg, None)

    @classmethod
    def success(cls, result: T):
        return cls(False, '', result)
