
from enum import Enum

from dataclasses import dataclass
from typing import Generic, TypeVar


T = TypeVar('T')


class Result(Enum):
    Success = 0
    Error = 1


@dataclass
class Option(Generic[T]):
    status: Result
    error_msg: str
    result: T

    @classmethod
    def error(cls, error_msg):
        return cls(Result.Error, error_msg, None)

    @classmethod
    def success(cls, result: T):
        return cls(Result.Success, '', result)
