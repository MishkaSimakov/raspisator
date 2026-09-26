"""LLDB data formatters for raspisator's exact-arithmetic types.

These decode the raw fields of BigInteger/Rational (see src/field/BigInteger.h)
instead of calling toString() in the inferior, so they also work in release
builds and stay fast for large containers of Rational.

Wiring: `command script import <path-to-this-file>` from ~/.lldbinit or the
LLDB console. CLion needs "Allow loading of .gdbinit/.lldbinit files" enabled
in Settings | Build, Execution, Deployment | Debugger.
"""

# BigInteger stores chunks little-endian in base 10^9 (kChunkWidth == 9).
kChunkWidth = 9


def _bigint_str(val):
    chunks = val.GetChildMemberWithName('chunks_')
    count = chunks.GetNumChildren()
    if count == 0:
        return '0'

    parts = [str(chunks.GetChildAtIndex(count - 1).GetValueAsSigned())]
    parts += ['%0*d' % (kChunkWidth, chunks.GetChildAtIndex(i).GetValueAsSigned())
              for i in range(count - 2, -1, -1)]
    result = ''.join(parts)

    is_positive = val.GetChildMemberWithName('is_positive_')
    if result.strip('0') != '' and not is_positive.GetValueAsUnsigned():
        result = '-' + result

    return result


def bigint_summary(val, internal_dict):
    return _bigint_str(val)


def rational_summary(val, internal_dict):
    numerator = _bigint_str(val.GetChildMemberWithName('numerator_'))
    denominator = _bigint_str(val.GetChildMemberWithName('denominator_'))

    if denominator == '1':
        return numerator

    try:
        approximation = ' (~%.6g)' % (int(numerator) / int(denominator))
    except (OverflowError, ZeroDivisionError):
        # value outside double's range, or a non-normalized/garbage Rational
        approximation = ''

    return '%s/%s%s' % (numerator, denominator, approximation)


def __lldb_init_module(debugger, internal_dict):
    for type_name, function in (('BigInteger', 'bigint_summary'),
                                ('Rational', 'rational_summary')):
        debugger.HandleCommand('type summary add -F %s.%s %s'
                               % (__name__, function, type_name))
