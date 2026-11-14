"""Exact binary linear algebra using Python integers as row bit vectors."""


def indices(row: int) -> list[int]:
    """Return the set-bit positions in ascending order."""
    row = int(row)
    if row < 0:
        raise ValueError("A binary row must be nonnegative")
    result = []
    while row:
        bit = row & -row
        result.append(bit.bit_length() - 1)
        row ^= bit
    return result


def echelon(rows) -> dict[int, int]:
    """Return a basis indexed by its highest set-bit pivot."""
    pivots = {}
    for row in rows:
        row = int(row)
        while row:
            pivot = row.bit_length() - 1
            if pivot in pivots:
                row ^= pivots[pivot]
            else:
                pivots[pivot] = row
                break
    return pivots


def reduce_row(row: int, pivots: dict[int, int]) -> int:
    """Reduce a row modulo the supplied echelon basis."""
    row = int(row)
    for pivot in sorted(pivots, reverse=True):
        if (row >> pivot) & 1:
            row ^= pivots[pivot]
    return row


def legacy_logicals(constraints, stabilizers, n: int) -> list[int]:
    """Construct a deterministic basis of ker(constraints) / stabilizers.

    Ascending-column elimination fixes the observable ordering used in the
    reference circuits. The returned logical operators need not be a paired
    symplectic X/Z basis: each memory experiment uses only one basis at a time.
    """
    low_pivots = {}
    for row in constraints:
        row = int(row)
        while row:
            pivot = (row & -row).bit_length() - 1
            if pivot in low_pivots:
                row ^= low_pivots[pivot]
            else:
                low_pivots[pivot] = row
                break
    quotient_pivots = echelon(stabilizers)
    result = []
    for column in range(n):
        if column in low_pivots:
            continue
        vector = 1 << column
        for pivot in sorted(low_pivots, reverse=True):
            if (low_pivots[pivot] & vector).bit_count() % 2:
                vector |= 1 << pivot
        reduced = reduce_row(vector, quotient_pivots)
        if reduced:
            result.append(vector)
            quotient_pivots[reduced.bit_length() - 1] = reduced
    return result
