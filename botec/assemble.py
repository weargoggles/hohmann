#!/usr/bin/env python
#
# $Id: //projects/botec/assemble.py#6 $ $Date: 2006/06/15 $

"""
Data file parser for BOTEC.
"""

__package__ = 'botec'


import md5
import os
import sys
import time

import botec


SOL_DIRECTORY = 'data/sol'

DEFAULT_PREFIX = 'main'
DATA_EXTENSION = '.data'
MAIN_FILENAME = DEFAULT_PREFIX + DATA_EXTENSION

NAME_KEY = 'name'


class Parser:
    def __init__(self, directory='.', factory=botec.System):
        self.filename = '<unknown>'
        self.directory = directory
        self.system = factory()
        self.hasher = md5.new()

    def has(self, key): return self.system.has_key(key)
    def get(self, key): return self.system[key]
    def set(self, key, value): self.system.add(key, value)

    def add(self, data):
        name = data[NAME_KEY]
        world = botec.World(self.system, data)
        self.set(name, world)
        return world
        
    def parse(self, filename=MAIN_FILENAME, remapping=None):
        pathname = os.path.join(self.directory, filename)
        inputFile = open(pathname)
        if remapping is None:
            remapping = {}
        else:
            remapping = remapping.copy()
        defaults = {}
        current = None
        for line in inputFile.xreadlines():
            self.hasher.update(line)
            if line.find('#') >= 0:
                line = line.split('#', 1)[0]
            if line.find(';') >= 0:
                lines = line.split(';')
            else:
                lines = [line]
            for line in lines:
                line = line.strip()
                if not line:
                    break
                if line.startswith('@'):
                    assert current is None
                    include = line[1:].strip()
                    if include.find('.') < 0:
                        include += DATA_EXTENSION
                    self.parse(include, remapping)
                elif line.find('=') >= 0:
                    key, value = line.split('=', 1)
                    key = key.strip()
                    value = value.strip()
                    remapping[key] = value
                elif line.find('{') >= 0: # disallow override
                    current = {}
                    current.update(defaults)
                    name = line.split('{', 1)[0].strip()
                    if remapping.has_key(name):
                        name = remapping[name]
                    if self.has(name): print name # assert not self.has_key(name), name
                    current[NAME_KEY] = name
                elif line.find('[') >= 0: # override
                    name = line.split('[', 1)[0].strip()
                    if remapping.has_key(name):
                        name = remapping[name]
                    current = self.get(name).rawData
                    current.update(defaults)
                elif line.find('}') >= 0:
                    assert current is not None
                    name = current[NAME_KEY]
                    self.add(current) ###
                    current = None
                elif line.find(']') >= 0:
                    ### This works but is really hackish, there should be a
                    ### better way of updating previously created objects.
                    assert current is not None
                    name = current[NAME_KEY]
                    world = self.get(name)
                    world.data = current.copy()
                    world.ensure(self.system)
                elif line.find(':') >= 0:
                    key, value = map(lambda x: x.strip(), line.split(':', 1))
                    assert key in botec.LEGAL_KEYS, key
                    if current is None:
                        defaults[key] = value
                    else:
                        current[key] = value
                else:
                    key = line
                    assert key in botec.LEGAL_KEYS, key
                    if current is None:
                        defaults[key] = True
                    else:
                        current[key] = True
        inputFile.close()

    def tag(self):
        now = time.time()
        date = time.strftime("%Y%m%d", time.localtime(now))
        stamp = int(now)
        count = len(self.system)
        hash = self.hasher.hexdigest()
        pid = os.getpid()
        data = [self.filename, date, count, hash, botec.__version__]
        result = ':'.join(map(str, data))
        return result

    def save(self, filename=botec.DEFAULT_FILENAME, \
             protocol=1, compressed=False):
        self.filename = filename
        botec.save(self.system, self.tag(), filename, protocol, compressed)

def main():
    parser = Parser(SOL_DIRECTORY)
    parser.parse()
    parser.save()
    if botec.DEBUG:
        print >> sys.stderr, "System: %r" % parser.system
        print >> sys.stderr, "Tag: %s" % parser.tag()

if __name__ == '__main__': main()
