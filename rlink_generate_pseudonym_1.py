#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Copyright (c) 2016-2018 CEA
#
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software. You can use,
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info".
#
# As a counterpart to the access to the source code and rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty and the software's author, the holder of the
# economic rights, and the successive licensors have only limited
# liability.
#
# In this respect, the user's attention is drawn to the risks associated
# with loading, using, modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean that it is complicated to manipulate, and that also
# therefore means that it is reserved for developers and experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and, more generally, to use and operate it in the
# same conditions as regards security.
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.

"""Generate pseudonyms to de-identify R-LiNK data.

These pseudonyms will be used to de-identify neuroimaging datasets and
biological samples at the source, in acquisition centres, before they
are sent to the relevant databank of biobank.

The algorithm is brain-dead. This is not reusable code, but it seems good
enough for a one-shot run.

Notes
-----
.. Damerau–Levenshtein distance
   https://en.wikipedia.org/wiki/Damerau%E2%80%93Levenshtein_distance

.. Lexicographic code
   https://en.wikipedia.org/wiki/Lexicographic_code

"""

from re import findall
from random import shuffle
from jellyfish import damerau_levenshtein_distance
from damm import encode

DIGITS = 5
MIN_DISTANCE = 3
PREFIXES = {
    '1000': 20,
    '1100': 20,
    '1200': 20,
    '1300': 20,
    '1400': 20,
    '1500': 20,
    '1600': 20,
    '1700': 20,
    '1800': 20,
    '1900': 20,
    '2000': 20,
    '2100': 20,
    '2200': 20,
    '2300': 20,
    '2400': 20,
}


def import_existing_pseudonyms(path):
    """Return a set of existign pseudonyms to avoid and remain distant enough from.

    Parameters
    ----------
    path : str
        Path of the file containing existing pseudonyms.

    Returns
    -------
    set
        Set of pseudonyms.

    """
    return set()


def pseudonym_generator(prefixes, digits, min_distance, existing=set()):
    """Generate "different enough" numeric pseudonyms (in the sense of Damerau-Levenshtein distance).

    The numeric pseudonyms are made of:
    - a prefix made from a `prefixes` key
    - a main code made of `digits` digits
    - a Damm check digit as a suffix

    Parameters
    ----------
    prefix : dict
        Keys are the expected prefixes - strings made of digits.
        Values represent the expected number of pseudonyms to generate
        for each prefix.

    digits : int
        Number of digits the main code is made of.

    min_distance : int
        Minimal Damerau-Levenshtein distance between generated strings.

    existing : set
        Set of existing, previously assigned pseudonyms. Newly generated
        pseudonyms will be added to `existing`.

    Returns
    -------
    dict
        Keys are the expected prefixes - strings made of digits.
        Values represent the generated pseudonyms for each prefix.

    """
    # avoid numbers starting with 0
    # for example for 5 digits, choose numbers between 10000 and 99999
    min_value = 10 ** (digits - 1)
    max_value = (10 ** digits) - 1

    # avoid more than 2 repeated consecutive characters
    candidates = list(x for x in range(min_value, max_value + 1)
                      if not findall(r'((\w)\2{2,})', str(x)))

    # randomize
    shuffle(candidates)

    generated = {x: [] for x in prefixes.keys()}
    for i in candidates:
        i = str(i)
        for prefix, n in prefixes.items():
            if len(generated[prefix]) < n:
                # append Damm decimal check digit
                code = i + str(encode(prefix + i))
                if code not in existing:
                    # calculate minimal Damerau-Levenshtein to other pseudonyms
                    distance = min((damerau_levenshtein_distance(code, l)
                                    for l in existing),
                                   default=min_distance)  # empty set
                    # keep pseudnyms "distant enough" from existing pseudonyms
                    if distance >= min_distance:
                        existing.add(code)
                        generated[prefix].append(code)
                        yield prefix + code
                        break  # from current inner loop over prefix


def main():
    existing = import_existing_pseudonyms(None)
    for pseudonym in pseudonym_generator(PREFIXES, DIGITS, MIN_DISTANCE, existing):
        print(pseudonym)


if __name__ == '__main__':
    main()
