class InvalidGenotypeString(Exception):
    pass


class Genotype:
    def __init__(self, call1, call2=None):
        self.call1 = str(call1) if call1 is not None else None
        self.call2 = str(call2) if call2 is not None else None

    def __eq__(self, other):
        x = set([call for call in self.call() if call is not None])
        y = set([call for call in other.call() if call is not None])
        return x == y

    def call(self):
        return self.call1, self.call2

    def is_null(self):
        return self.call1 == "." or self.call1 is None

    def is_heterozygous(self):
        if self.is_null() or self.call1 == self.call2 or self.call2 is None:
            return False
        else:
            return True

    def is_reference(self):
        return str(self.call1) == "0" and (
            self.call1 == self.call2 or self.call2 is None
        )

    def is_alt(self):
        return (
            not self.is_heterozygous()
            and not self.is_reference()
            and not self.is_null()
        )

    @staticmethod
    def from_string(s):
        """Parse a VCF string into a genotype. Expected format is '0/0'"""
        delimiter = "/"
        if s.count(delimiter) != 1:
            raise InvalidGenotypeString(
                "Invalid genotype string received: {}\nGenotype string should be of the form 0/0".format(
                    s
                )
            )
        calls = s.strip().split("/")
        call1 = calls[0]
        call2 = calls[1] or None
        return Genotype(call1, call2)
