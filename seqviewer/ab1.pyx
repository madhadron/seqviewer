"""
ab1.py - Read Applied BioSystems AB1 files
Copyright 2011, Frederick Ross

"""
import collections
import contextlib
import threading
import cStringIO
import datetime
import struct
import pytest
import math
import sys

major_version = 1

def packed_handle(fmt, *args):
    """Wrap a call to struct.pack in a temporary handle.

    The struct module works with strings, while all the functions in
    this module work on handles, so this function resolves the
    mismatch.
    """
    return cStringIO.StringIO(struct.pack(fmt, *args))


def parse_header(h):
    """Parse the header of an AB1 file.

    Returns the version number, the offset to the directory of entries
    in the AB1 file, and the number of entries in that directory as a
    dictionary (with keys 'version', 'numelements', and 'dataoffset').
    """
    # AB1 files begin with four bytes encoding the ASCII characters 'ABIF'
    abif = h.read(4)
    if abif != "ABIF":
        raise ValueError("AB1 file does not begin with 'ABIF'")
    # Next come two bytes encoding the file format version (multiplied
    # by 100 so that 1.01 is stored as the integer 101).
    version, = struct.unpack('>h', h.read(2))
    if math.floor(version/100) != 1:
        raise ValueError("AB1 file version is not supported by this "
                         "library (found: %f, must be 1.xx)" % (version/100.0))
    # Then comes a directory entry pointing to the offset containing
    # all the entries for this file.
    tdir = parse_direntry(h)
    return {'version': version,
            'numelements': tdir['numelements'],
            'dataoffset': tdir['dataoffset']}

def test_header():
    valid_header = packed_handle('>cccc h ccccihhiiii',
                                 'A','B','I','F', 101,
                                 't', 'd', 'i', 'r',
                                 1, 1023, 28, 10, 280, 10233, 30)
    assert parse_header(valid_header) == {'version': 101,
                                          'numelements': 10,
                                          'dataoffset': 10233}
    with pytest.raises(ValueError):
        bad_abif = packed_handle('>cccch', 'A', 'B', 'Q', 'R', 101)
        parse_header(bad_abif)
        bad_version = packed_handle('>cccch', 'A','B','I','F', 205)


# Next we define functions to decode all the supported data types in
# the file.  All such decoder functions take two arguments, a handle
# to read from and an integer giving the number of bytes to read, and
# read a tuple of two values, the first being the value read and the
# second the number of bytes consumed in doing so.

# For types which decode to atomic objects like a single int, we use a
# factory function to produce the decoders.
def reader(fmt):
    """Return a function to decode according to *fmt*.

    The *fmt* argument is a string formatted according the usage of
    the ``struct`` module.  For example, to create a function to read
    a 32 bit, big Endian integer, you would call reader('>i').  The
    resulting function takes two arguments: a handle to read and the
    number of bytes which should be consumed.  It returns a tuple of
    the value read and the number of bytes actually consumed.  For
    instance,::

        reader('>i')(packed_handle('>i', 1066))

    evaluates to ``(1066, 4)``.
    """
    n = struct.calcsize(fmt)
    def reader(h, size):
        val = struct.unpack(fmt, h.read(n))[0]
        return (val, n)
    return reader

def test_reader():
    assert reader('>i')(packed_handle('>i', 1066), 4) == (1066, 4)


def read_date(h, size):
    """Read a date into a ``datetime.date`` object."""
    assert size == 4
    (year,month,day) = struct.unpack('>hBB', h.read(4))
    return (datetime.date(year,month,day), 4)

def test_read_date():
    d = datetime.date(1066,6,20)
    h = packed_handle('>hBB', d.year, d.month, d.day)
    assert read_date(h, 4) == (d, 4)


def read_time(h, size):
    """Read a time into a ``datetime.time`` object."""
    assert size == 4
    (hour,minute,second,centisecond) = struct.unpack('>BBBB', h.read(4))
    # time's last field is microseconds, not centiseconds
    return (datetime.time(hour,minute,second,centisecond*10), 4)

def test_read_time():
    t = datetime.time(13,22,56,120)
    h = packed_handle('>BBBB', t.hour, t.minute, t.second, t.microsecond/10)
    assert read_time(h, 4) == (t, 4)


def read_pstring(h, size):
    """Read a Pascal style string from handle *h*.

    The first byte read is taken to be the length of the string.  The
    *size* argument should be 1 for valid files, but is ignored for
    this function.  The number of bytes returned is the number of
    characters plus one for the byte giving the length.
    """
    assert size == 1
    first_byte = h.read(1)
    if first_byte == '':
        raise ValueError('Failed to read length for pstring')
    n_chars = struct.unpack('>B', first_byte)[0]
    return (h.read(n_chars), n_chars+1)

def test_read_pstring():
    s = 'abcdef'
    assert len(s) < 256
    h = cStringIO.StringIO(struct.pack('>B', len(s)) + s)
    assert read_pstring(h, 1) == (s, 7)

def read_cstring(h, size):
    """Read a C style string from handle *h*.

    The string continues until a terminating \x00 character is read.
    The length returned includes the terminating character.
    """
    accumulated_string = ''
    while True:
        next_char = h.read(1)
        if next_char == '\0' or next_char == '': # End of the string
            break
        accumulated_string += next_char
    return (accumulated_string, len(accumulated_string)+1)

def test_read_cstring():
    s = 'abcdef'
    h = cStringIO.StringIO(s + '\0')
    assert read_cstring(h, 1) == (s, 7)


# Legacy supported data types
def read_thumb(h, size):
    assert size == 10
    return (struct.unpack('>iiBB', h.read(10)), 10)

def read_bool(h, size):
    assert size == 1
    return (h.read(1) != 0, 1)

def read_user(h, size):
    """Returns raw strings of bytes."""
    return (h.read(size), size)

# We use a lookup table on the elementtype field to decide how to
# decode a given element.
_element_decoders = {1: reader('>B'), # unsigned byte
                     2: reader('>b'), # signed byte
                     3: reader('>H'), # unsigned short
                     4: reader('>h'), # signed short
                     5: reader('>i'), # signed 4 byte int
                     7: reader('>f'), # 4 byte float
                     8: reader('>d'), # 8 byte double precision float
                     10: read_date,  # Date (year,month,day)
                     11: read_time,   # Time (hour,minute,second,centisecond)
                     18: read_pstring, # Pascal-style strings
                     19: read_cstring, # C-style strings
                     12: read_thumb, # Thumbprints (obsolete)
                     13: read_bool, # Booleans (obsolete)
                     1024: read_user # User defined data type
                     }


def parse_direntry(h):
    """Parse a DirEntry from the AB1 file."""
    # The format is a C struct of:
    # struct DirEntry {
    #     char[4] name;
    #     signed int32 number;
    #     signed int16 elementtype;
    #     signed int16 elementsize;
    #     signed int32 numelements;
    #     signed int32 datasize;
    #     signed int32 dataoffset;
    #     signed int32 datahandle;
    # }
    fmt = ">ccccihhiiii"
    fields = struct.unpack(fmt, h.read(28))
    direntry = {'name': fields[0] + fields[1] + fields[2] + fields[3],
                'number': fields[4],
                'elementtype': fields[5],
                'elementsize': fields[6],
                'numelements': fields[7],
                'datasize': fields[8],
                'dataoffset': fields[9]}
    # If the data in this entry is no larger than 4 bytes, it is
    # stored directly in the dataoffset field instead of putting it
    # elsewhere and storing a pointer to it.
    if direntry['datasize'] <= 4:
        # We rencode the dataoffset field to be decoded by the
        # appropriate function.
        h2 = packed_handle('>i', direntry['dataoffset'])
        direntry['data'] = read_sequence(direntry, h2)
    assert direntry['datasize'] >= direntry['numelements'] * direntry['elementsize']
    return direntry

def test_parse_direntry():
    tdir = packed_handle('>ccccihhiiii',
                         't', 'd', 'i', 'r',
                         1, 1023, 28, 10, 280, 10233, 30)
    assert parse_direntry(tdir) == {'name': 'tdir', 'number': 1,
                                    'elementtype': 1023, 'elementsize': 28,
                                    'numelements': 10, 'datasize': 280,
                                    'dataoffset': 10233}


def read_sequence(direntry, handle):
    """Read the data encoded by *direntry* from *handle.

    The *direntry* should be a dictionary of the kind returned by
    ``parse_direntry``.  It properly loops through the provided
    handle, reading elements of type specified in direntry, then
    returns those elements as a list.
    """
    decoder = _element_decoders[direntry['elementtype']]
    values = []
    bytes_read = 0
    while bytes_read < direntry['datasize']:
        (val, n_read) = decoder(handle, direntry['elementsize'])
        values.append(val)
        bytes_read += n_read
    return values


class SequenceWithConfidence(object):
    def __init__(self, sequence, confidences):
        assert len(sequence) == len(confidences)
        self.sequence = sequence
        self.conf = confidences
    def __len__(self):
        return len(self.sequence)
    def __getitem__(self, i):
        return SequenceWithConfidence(self.sequence[i], [self.conf[i]])
    def __getslice__(self, i, j):
        return SequenceWithConfidence(self.sequence[i:j], self.conf[i:j])
    def __str__(self):
        return ' '.join(['%s/%s' % (b,c) for b,c in zip(self.sequence, self.conf)])
    def __repr__(self):
        return ' '.join(['%s/%s' % (b,c) for b,c in zip(self.sequence, self.conf)])
    def __add__(self, other):
        if isinstance(other, SequenceWithConfidence):
            return SequenceWithConfidence(self.sequence + other.sequence,
                                          self.conf + other.conf)
        elif isinstance(other, basestring):
            return self.sequence + other
        else:
            raise ValueError("Invalid thing to add to SequenceWithConfidence")
    def reverse(self):
        return SequenceWithConfidence(self.sequence[::-1], self.conf[::-1])
    def complement(self):
        new_seq = self.sequence. \
            replace('A','t'). \
            replace('T','a'). \
            replace('C','g'). \
            replace('G','c').upper()
        return SequenceWithConfidence(new_seq, self.conf[::-1])
    def revcomp(self):
        return self.reverse().complement()
    def high_quality_segment(self, threshold=40):
        high_quality_positions = [i for i,x in enumerate(self.conf)
                                  if x > threshold]
        if len(high_quality_positions) < 10:
            return SequenceWithConfidence('',[])
        else:
            high_quality_positions.reverse()
            while True:
                if len(high_quality_positions) < 10:
                    return SequenceWithConfidence('',[])
                start = high_quality_positions.pop()
                if all([x > threshold for x in self.conf[start:(start+5)]]):
                    break
            high_quality_positions.reverse()
            while True:
                if len(high_quality_positions) < 10:
                    return SequenceWithConfidence('',[])
                end = high_quality_positions.pop()
                if all([x > threshold for x in self.conf[(end-5):end]]):
                    break
            return self[start:end]

def reapply_confidences(seq, oldseq):
    old_conf = list(reversed(oldseq.conf)) # pop goes from the end in Python
    new_conf = []
    for c in seq:
        if c == '-':
            new_conf.append(None)
        else:
            new_conf.append(old_conf.pop())
    return SequenceWithConfidence(seq, new_conf)



class Ab1File(object):
    """Open an AB1 file for reading.

    The ``Ab1File`` object imitates a read only dictionary by
    providing ``keys()`` and lookup with ``[]``.  The result is always
    a list of lists.  The inner lists contain the actual values, the
    outer lists are each instance of a key in the file, in order.

    This class is thread safe *if* no other code uses the handle it
    uses.  Internally, all operations properly lock the handle before
    seeking to fetch their data.

    The class does not cache data read from disk.
    """
    def __init__(self, handle):
        self.handle = handle
        self.lock = threading.Lock()
        ab1_header = parse_header(self.handle)
        self.version = ab1_header['version']
        # Read enough information on each entry from disk to be able
        # to fetch it on demand.  The name of the entries is not
        # unique, so we store them in an association list.
        self.entries = collections.defaultdict(list)
        self.handle.seek(ab1_header['dataoffset'])
        for i in range(ab1_header['numelements']):
            direntry = parse_direntry(self.handle)
            key = direntry['name']
            self.entries[key].append(direntry)
        self.base_order = dict([(chr(x).upper(),i)
                                for i,x in enumerate(self['FWO_'][0])])

    def keys(self):
        return self.entries.keys()

    def __len__(self):
        """Return the number of entries in the file.

        Note that this is *not* the same as the number of distinct
        keys.  It is the number of instances of each key.
        """
        return sum(len(self.entries[k]) for k in self.keys())

    def __getitem__(self, key):
        """Fetch the values listed under name *key* as a list.

        The values are fetched from disk each time they are requested.
        No caching is attempted.
        """
        if not(key in self.entries):
            raise KeyError("No such key %s in AB1 file" % key)
        data = []
        for direntry in self.entries[key]:
            if 'data' in direntry: # The data is <= 4 bytes long, and stored locally
                data.append(direntry['data'])
            else: # The data must be fetched from disk
                try:
                    self.lock.acquire()
                    self.handle.seek(direntry['dataoffset'])
                    data.append(read_sequence(direntry, self.handle))
                finally:
                    self.lock.release()
        return data
    def bases(self):
        s = ''
        for x in self['PBAS'][0]:
            s += chr(x)
        return s
    def base_centers(self):
        return self['PLOC'][0]
    def base_confidences(self):
        return self['PCON'][0]
    def bases_with_confidence(self):
        return SequenceWithConfidence(self.bases(), self.base_confidences())
    def trace(self, base):
        assert base in 'ACTG'
        return self['DATA'][self.base_order[base] + 8]
    def _raw_trace(self, base):
        assert base in 'ACTG'
        return self['DATA'][self.base_order[base]]
    def signal_strength(self, base):
        assert base in 'ACTG'
        return self['S/N%'][0][self.base_order[base]]
    def run_date(self):
        return self['RUND'][0][0]
    def run_time(self):
        return self['RUNT'][0][0]

def test_ab1file():
    try:
        h = open('scrapings/tmpClCBJg-1.ab1', 'rb')
        a = Ab1File(h)
        assert len(a) == 122
        assert a['Scan'] == [[13960]]
        assert a['CTTL'] == [['Comment:']]
        assert isinstance(a.bases(), str)
        assert isinstance(a.base_centers(), list)
        assert isinstance(a.base_confidences(), list)
        assert all([x > 0 and x < 100 for x in a.base_confidences()])
        assert isinstance(a.trace('A'), list)
        assert isinstance(a._raw_trace('A'), list)
    finally:
        h.close()

