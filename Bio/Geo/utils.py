# Copyright 2001 by Katharine Lindner.  All rights reserved.
# Copyright 2006 by PeterC.  All rights reserved.
# Copyright 2007 by Michiel de Hoon.  All rights reserved.
# Copyright 2011 by Phillip Garland <pgarland@gmail.com>.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# With the exception of the 'table_begin' and 'table_end' attributes, and
# Annotation entity indicators, all entity indicator, entity attribute, and data
# table row lines should have a label and a value. The label and value are
# separated by an '=' character.
def _read_key_value(line):
    words = line[1:].split("=", 1)
    try:
        key, value = words
        value = value.strip()
    except ValueError:
        key = words[0]
        value = ""
    key = key.strip()
    return key, value

def stringIsType(strObj, coercion):
    """Test if strObj can be converted by coercion, e.g. int or float"""
    try:
        coercion(strObj)
    except ValueError:
        return False 

    return True

def maybeConvertToNumber(strObj):
    """Convert str to a numeric type, if possible."""
    if stringIsType(strObj, int):
        return int(strObj)
    elif stringIsType(strObj, float):
        return float(strObj)
    else:
        return strObj
