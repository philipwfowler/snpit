from typing import Optional, List


class InvalidGenotypeString(Exception):
    pass


class UnexpectedGenotypeError(Exception):
    pass


class Genotype:
    def __init__(self, call1: Optional[int], call2: Optional[int] = None):
        self.call1 = call1 if call1 is not None else -1
        self.call2 = self.call1 if call2 is None else call2

    def __eq__(self, other):
        return sorted(self.call()) == sorted(other.call())

    def call(self) -> List[int]:
        return [self.call1, self.call2]

    def is_null(self) -> bool:
        return set(self.call()) == {-1}

    def is_heterozygous(self) -> bool:
        return self.call1 != self.call2

    def is_reference(self) -> bool:
        return self.call1 == 0 and (self.call1 == self.call2 or self.call2 == -1)

    def is_alt(self) -> bool:
        return (
            not self.is_heterozygous()
            and not self.is_reference()
            and not self.is_null()
        )

    @staticmethod
    def from_string(s: str) -> "Genotype":
        """Parse a VCF string into a genotype. Expected format is '0/0'"""
        delimiter = "/"
        if s.count(delimiter) != 1:
            raise InvalidGenotypeString(
                """Invalid genotype string received: {s}\nGenotype string should be 
                of the form 0/0""".format(
                    s=s
                )
            )

        calls = s.strip().replace(".", "-1").split("/")
        call1 = -1 if calls[0] == "." else int(calls[0])
        call2 = None if calls[1] == "" else int(calls[1])
        return Genotype(call1, call2)
