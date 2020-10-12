#!/usr/bin/env python
#
# $Id: //projects/botec/botec.py#54 $ $Date: 2009/06/25 $

"""
BOTEC is a simple astrophysical and orbital mechanics calculator,
including a database of all named Solar System objects.
"""

__program__ = 'botec'
__version__ = '0.3.5'
__url__ = 'http://www.alcyone.com/software/botec/'
__author__ = 'Erik Max Francis <max@alcyone.com>'
__copyright__ = 'Copyright (C) 2003-2009 Erik Max Francis'
__license__ = 'GPL'


import gzip
import math
import os
import re
import sys
import time
import types
import logging

try:
    import cPickle as pickle
except ImportError:
    import pickle

DEBUG = True

DEFAULT_FILENAME = 'default.botec'
COMPRESSED_EXTENSION = '.gz'

# Dimensionless constants
ONE_THIRD = 1.0/3.0
TWO_FIFTHS = 2.0/5.0
THREE_FIFTHS = 3.0/5.0
PI = math.pi
PI_OVER_TWO = PI/2.0
PI_OVER_THREE = PI/3.0
PI_OVER_TWELVE = PI/12.0
TWO_PI = 2.0*PI
FOUR_PI = 4.0*PI
FOUR_PI_SQUARED = 4.0*PI**2
FOUR_THIRDS_PI = 4.0*PI/3.0
PI_SQUARED_OVER_SIX = PI**2/6.0
FIFTH_ROOT_OF_100 = 100.0**(1/5.0)

# Dimensioned constants
C = 299792458.0 # m/s [exact]
G = 6.67e-11 # N m^2/kg^2
TWO_G = 2*G # N m^2/kg^2
SIGMA = 5.671e-8 # W/(m^2 K)
JOULES_TO_TONNES_TNT = 1/4.184e9 # J/t TNT
RADIANS_TO_DEGREES = 180.0/PI # deg/rad
DEGREES_TO_RADIANS = 1.0/RADIANS_TO_DEGREES # rad/deg
INTENSITY_OF_MAGNITUDE_ZERO = 2.66765e-8 # W/m^2

try:
    INFINITY = 1e1000 # this should overflow in IEEE
except:
    INFINITY = 1e100 # to be safe; just a Very Large Number (VLN)


# Conversions.
def degreesToRadians(degrees):
    """Convert an angle in degrees to radians."""
    return degrees*DEGREES_TO_RADIANS

def radiansToDegrees(radians):
    """Convert an angle in radians to degrees."""
    return radians*RADIANS_TO_DEGREES

def intensityToApparentMagnitude(intensity):
    """Convert an intensity [W/m^2] to an apparent magnitude."""
    return -math.log(intensity/INTENSITY_OF_MAGNITUDE_ZERO,
                     FIFTH_ROOT_OF_100)

LEGAL_KEYS = (
    'albedo', # visual geometric albedo
    'alternate', # list of alternate names
    'arcs', # ring consists of arcs
    'atmospheric composition',
    'chaotic', # rotation period is chaotic
    'class',
    'designation', # for comets
    'diameter',
    'distance', # distance from primary
    'eccentricity',
    'ellipticity', # flattening
    'family', # asteroid family association
    'inclination', #
    'inclined', # orbital plane is highly inclined
    'luminosity',
    'mass',
    'nonspherical', # object is highly nonspherical
    'normalized moment of inertia', # I/(M R^2)
    'number', # for asteroids
    'oblate', # object is highly oblate
    'oblique', # angle of rotation axis is highly oblique
    'obliquity', # angle of rotation axis relative to orbital plane normal
    'optical depth',
    'period', # rotation period
    'pressure', # pressure at "surface"
    'primary',
    'radius', # volumetric mean radius
    'retrograde', # object orbits retrograde
    'scale height',
    'spectral type', # spectral type for stars
    'pressure',
    'surfaceless', # object has no concrete surface; figures for cloud tops
    'synchronous', # object rotates in synchrony with its revolution
    'temperature', # temperature at "surface"
    'type',
    'zone', # asteroid zone association
)

NUMERIC_KEYS = ('distance', 'eccentricity', 'inclination', 'mass', 
                'radius', 'period', 'temperature', 'luminosity', 
                'pressure', 'scale height', 'optical depth', 'albedo', 
                'obliquity', 'ellipticity', 
                'normalized moment of inertia')

def defaultNumericKey(key):
    if key == 'albedo':
        return 1.0
    else:
        return 0.0

STRING_KEYS = ('name', 'type', 'class', 'zone', 'family', 'designation',
               'atmospheric composition')

#
# ReprMixin
#

class ReprMixin(object):

    """A common mixin class for handling reprs."""
    
    def __repr__(self):
        return "<%s @ 0x%x (%s)>" % \
               (self.__class__.__name__, id(self), str(self))

#
# SI
#

class SI(ReprMixin):

    """A class representing a value and an associated unit, capable of
    doing conversions to base SI units."""

    BASES = ['m/s', 'm', 'kg', 's', 'A', 'K', 'rad', 'Pa', 'W', 'sr']
    PREFIXES = {'y': 1e-24,
                'z': 1e-21,
                'a': 1e-18,
                'f': 1e-15,
                'p': 1e-12, 
                'n': 1e-9,
                'u': 1e-6,
                'm': 1e-3, 
                'c': 1e-2,
                'd': 1e-1,
                '': 1,
                'da': 1e1,
                'h': 1e2, 
                'k': 1e3,
                'M': 1e6,
                'G': 1e9, 
                'T': 1e12,
                'P': 1e15,
                'E': 1e18, 
                'Z': 1e21,
                'Y': 1e24}
    OTHERS = {'au': (149.60e9, 'm'),
              'h': (60*60, 's'),
              'd': (60*60*24, 's'),
              'mo': (60*60*24*30, 's'),
              'y': (60*60*24*365.25, 's'), 
              'deg': (DEGREES_TO_RADIANS, 'rad'), 
              'msol': (1.989e30, 'kg'),
              'bar': (1e5, 'Pa'),
              'gee': (9.80665, 'm/s^2')}
    ALTERNATES = {'m': ['km', 'au'],
                  'rad': ['deg'],
                  's': ['h', 'd', 'mo', 'y'],
                  'm/s': ['km/s'],
                  'm/s^2': ['gee']}

    significantDigits = 3

    #@staticmethod
    def Convert(raw):
        """Take a raw unit, do any necessary conversions, extract
        prefixes, and convert it to a 2-tuple of the converted value
        and the base unit."""
        # First check for non-SI units.
        if raw in SI.OTHERS:
            return SI.OTHERS[raw]
        # Next look for SI units.
        for unit in SI.BASES:
            if raw.endswith(unit):
                prefix, base = raw[:-len(unit)], raw[-len(unit):]
                break
        else:
            prefix, base = '', raw
        return SI.PREFIXES[prefix], base
    Convert = staticmethod(Convert)

    def __init__(self, value, unit=None):
        if type(value) is str:
            assert unit is None
            if value.find(' ') >= 0:
                valueS, rawUnit = value.split(' ', 1)
            else:
                valueS, rawUnit = value, None
            if valueS.find('..') >= 0:
                lower, higher = map(float, valueS.split('..', 1))
                value = (lower + higher)/2.0
            elif valueS.find('+-') >= 0:
                center, width = map(float, valueS.split('+-', 1))
                value = center
            elif valueS.startswith('~'):
                value = float(valueS[1:])
            else:
                value = float(valueS)
            if rawUnit is not None:
                factor, base = self.Convert(rawUnit)
            else:
                factor, base = 1, None
            self.value = value*factor
            self.unit = base
        elif isinstance(type(value), SI):
            self.value = value.value
            self.unit = value.unit
        else:
            self.value = float(value)
            self.unit = unit

    def __float__(self):
        return self.value

    def __str__(self):
        choices = []
        bestChoice = None
        if self.unit in self.ALTERNATES:
            for alternate in self.ALTERNATES[self.unit]:
                factor, oldBase = self.Convert(alternate)
                newValue = self.value/factor
                choices.append((newValue, alternate))
            for choice in choices:
                if abs(choice[0]) >= 1.0 and abs(choice[0]) < 1000.0:
                    bestChoice = choice
        if bestChoice is not None:
            # Figure out how many digits to show.
            bestValue, bestUnit = bestChoice
            leftSide = int(abs(bestValue))
            if leftSide == 0:
                leftDigits = 0
            else:
                leftDigits = len(str(leftSide))
            rightDigits = self.significantDigits - leftDigits
            actualRightDigits = max(rightDigits, 0)
            equiv = (" = %%.%df %%s" % actualRightDigits) % (round(bestValue, rightDigits), bestUnit)
        else:
            equiv = ""
        if self.unit is None:
            return (("%%.%de" % (self.significantDigits - 1,)) % (self.value) + equiv)
        else:
            return (("%%.%de %%s" % (self.significantDigits - 1,)) % (self.value, self.unit) + equiv)

#
# Location
#

class Location(ReprMixin):

    """A location represents an orbital position around a primary
    object at a given distance."""
    
    def __init__(self, primary, distance):
        assert distance >= 0.0
        assert primary is None or isinstance(primary, World)
        self.__primary = primary
        self.__distance = distance

    def primary(self): return self.__primary
    def distance(self): return self.__distance

    def phaseAngle(self):
        """The base phase angle. [rad]"""
        return 0.0

    def altitude(self):
        """The altitude is the height above the surface. [m]"""
        return self.__distance - self.__primary.radius()

    def primaries(self, seq=None):
        """Return a list of primaries, from most distant to least distant."""
        if seq is None:
            top = 1
            seq = []
        else:
            top = 0
        if self.__primary is not None:
            seq.append(self.__primary)
            self.__primary.primaries(seq)
        if top:
            return seq

    def hierarchy(self):
        """A list of the objects that are enumerated as one moves up the
        hierarchy, starting with this object."""
        return [self] + self.primaries()

    def isDescendentOf(self, other):
        """Is this object an descendent (child, grandchild, etc.) of the
        other object?"""
        location = self
        while location is not None:
            if location is other:
                return True
            location = location.primary()
        return False

    def isAncestorOf(self, other):
        """Is this object an ancestor (parent, grandparent, etc.) of the
        other object?"""
        return other.isDescendentOf(self)

    def findCommonAncestor(self, other):
        """Find the lowest world in the hierarchy that is an ancestor
        to both objects."""
        common = None
        selfPrimaries = self.primaries()
        selfPrimaries.reverse()
        otherPrimaries = other.primaries()
        otherPrimaries.reverse()
        for p, q in zip(selfPrimaries, otherPrimaries):
            if p is not q:
                break
            common = p # either one, they're identical
        assert common is not None
        return common

    def distanceFrom(self, other):
        """Compute the (average) distance between two bodies, either both of
        which have the same parent or one is a parent of the other. [m]"""
        if self.primary() is other:
            return self.distance()
        elif other.primary() is self:
            return other.distance()
        elif self.primary() is other.primary():
            return abs(self.distance() - other.distance())
        else:
            assert False

    def orbitalCircumference(self):
        """The circumference of a circular orbit at this distance. [m]"""
        return TWO_PI*self.distance()

    def orbitalAngularSpeedAroundPrimary(self):
        """The orbital angular speed around the primary. [rad/s]"""
        return self.__primary.orbitalAngularSpeedAtDistance(self.__distance)

    def orbitalSpeedAroundPrimary(self):
        """The orbital speed around the primary. [m/s]"""
        return self.__primary.orbitalSpeedAtDistance(self.__distance)

    def orbitalPeriodAroundPrimary(self):
        """The orbital period around the primary. [s]"""
        return self.__primary.orbitalPeriodAtDistance(self.__distance)

    def escapeSpeedFromPrimary(self):
        """The speed required to escape the primary. [m/s]"""
        return self.__primary.escapeSpeedFromDistance(self.__distance)

    def escapeSpeedExcess(self):
        """The difference between the escape speed and the orbital speed.
        [m/s]"""
        return (self.escapeSpeedFromPrimary() - 
                self.orbitalSpeedAroundPrimary())

    def nominalSpeed(self):
        """This is the nominal speed of the object at this location.
        For orbits it is orbital speed; for non-orbits it is the
        instantaneous speed.  [m/s]"""
        return self.orbitalSpeedAroundPrimary()

    def closingSpeeds(self):
        """Return a 2-tuple providing the range of possible impactor closing
        speeds given that the impactor is travelling at escape speed when
        it hits. [m/s, m/s]"""
        escape = self.escapeSpeedFromPrimary()
        orbital = self.orbitalSpeedAroundPrimary()
        assert escape > orbital
        return escape - orbital, escape + orbital

    def massicAngularMomentum(self):
        """Orbital angular momentum per unit mass. [m^2/s]"""
        # l/m = r v
        return self.__distance*self.orbitalSpeedAroundPrimary()

    def massicKineticEnergy(self):
        """Orbital kinetic energy per unit mass. [J/kg]"""
        # K/m = (1/2) v^2
        return 0.5*self.orbitalSpeedAroundPrimary()**2

    def massicPotentialEnergy(self):
        """Orbital potential energy per unit mass. [J/kg]"""
        # U/m = -G M/r
        return -G*self.__primary.mass()/self.distance()

    def massicTotalEnergy(self):
        """Orbital total energy per unit mass. [J/kg]"""
        # E/m = K/m + U/m
        return self.massicKineticEnergy() + self.massicPotentialEnergy()

    def effectiveDeltaveeDifference(self, other):
        """The effective (ideal) deltavee between the two locations
        around the same primary based purely on energy arguments.
        [m/s]"""
        # Estimating deltavee difference by the difference in massic kinetic
        # energy.
        assert self.primary() is other.primary()
        massicEnergyDifference = abs(self.massicTotalEnergy() - 
                                     other.massicTotalEnergy())
        return math.sqrt(2.0*massicEnergyDifference)

    def relativeOrbitalAngularSpeedWith(self, other):
        """The relative orbital angular speed of the two objects [rad/s]."""
        # deltaomega = omega_2 - omega_1
        assert self.primary() is other.primary()
        return (other.orbitalAngularSpeedAroundPrimary() - 
                self.orbitalAngularSpeedAroundPrimary())

    def opportunityDisplacementAngleWith(self, other):
        """The opportunity angle that this body makes with the other.
        Note that the angle can be greater than 2 pi or less than 2 pi
        to indicate multiple opportunity orbits are required. [rad]"""
        # PI - omega t_Hohmann
        assert self.primary() is other.primary()
        transfer = HohmannTransfer(self, other)
        theta = other.orbitalAngularSpeedAroundPrimary()*transfer.duration()
        return PI - theta

    def opportunityAngleWith(self, other):
        """The opportunity angle that this body makes with the other.
        The angle is normalized to be in the range (-2 pi, 2
        pi). [rad]"""
        theta = self.opportunityDisplacementAngleWith(other)
        if theta > TWO_PI or theta < -TWO_PI:
            theta %= TWO_PI
        return PI - theta

    def opportunitiesPeriodWith(self, other):
        """The period at which opportunities recur between the two bodies.
        [s]"""
        # 2 PI/|deltaomega|
        return TWO_PI/abs(self.relativeOrbitalAngularSpeedWith(other))

    def roundTripWaitTimeWith(self, other):
        """For a round trip of Hohmann transfers from self to other, a
        wait for opportunity, and a return from other to self, the
        time required to wait for the return opportunity. [s]"""
        transfer = HohmannTransfer(self, other)
        transferDuration = transfer.duration()
        opportunityAngle = other.opportunityAngleWith(self)
        relativeAngularSpeed = self.relativeOrbitalAngularSpeedWith(other)
        relativeAngle = PI + relativeAngularSpeed*transferDuration
        # If we missed this opportunity, we have to wait for another cycle.
        if relativeAngularSpeed > 0.0:
            if relativeAngle > opportunityAngle:
                relativeAngle -= TWO_PI
        else:
            if relativeAngle < opportunityAngle:
                relativeAngle += TWO_PI
        waitingAngle = abs(opportunityAngle - relativeAngle)
        waitingDuration = waitingAngle/abs(relativeAngularSpeed)
        return waitingDuration

    def roundTripTimeWith(self, other):
        """The round trip time for a departure from self, arrival at
        other, wait for opportunity to return, and transfer from other
        back to self. [s]"""
        # Hohmann transfer times are symmetric, so this is the way delay
        # plus twice the Hohmann transfer duration.
        waitingDuration = self.roundTripWaitTimeWith(other)
        transferDuration = HohmannTransfer(self, other).duration()
        return waitingDuration + 2*transferDuration

    def instability(self):
        """The inherent instability of this orbital location, measured
        as the e-folding time of the positional error.  None indicates
        no instability (an infinite e-folding time). [s, or None]"""
        return None

    # TODO: General orbital instability might be reasonably quickly
    # estimated by doing a numerical integration of the gravitational
    # effects of nearby/large masses (keeping in mind the relative
    # angular speeds).
    
    def stationKeeping(self):
        """The amount of station keeping that must be exerted to keep a
        body from drifting here, measured as an acceleration. [m/s^2]"""
        # This is estimated from a simple exponential position law
        # tending to draw particles away from the point, and analyzing the
        # average acceleration of the particle which needs to be cancelled.
        if self.instability() is None:
            return 0.0
        # a = r_0/tau^2
        R_0 = 100e3
        return R_0/self.instability()**2

    def pressureAt(self):
        """The pressure from the primary at this altitude. [Pa]"""
        return self.primary().pressureAtAltitude(self.altitude())

    def __str__(self):
        return "%s around `%s'" % \
               (SI(self.__distance, 'm'), str(self.__primary))

class StationaryLocation(Location):

    """A stationary location is one that is in an orbit that has the
    same period as the primary's rotation period, and is circular."""

    def __init__(self, primary):
        Location.__init__(self, primary, primary.apostationaryDistance())


class AltitudeLocation(Location):

    """An altitude location is simply a location expressed in terms of
    altitude (distance from surface) rather than radius (distance from
    center)."""
    
    def __init__(self, primary, height):
        assert primary.radius() > 0.0
        Location.__init__(self, primary, primary.radius() + height)


class GrazingLocation(AltitudeLocation):

    """A grazing location is a location that grazes the top of a
    world's surface or atmosphere."""

    def __init__(self, primary):
        AltitudeLocation.__init__(self, primary, 0.0)


class World(Location):

    """A world represents a physical object in the system."""
    
    def __init__(self, system, data=None, **keywords):
        if data is None:
            data = keywords
        self.data = data
        self.rawData = data.copy()
        primary, distance = self.ensure(system)
        Location.__init__(self, primary, distance)
        self.__secondaries = []
        if primary is not None:
            primary.__secondaries.append(self)

    def ensure(self, system):
        """Normalize all the figures provided; return the primary (object)
        and distance (in metres) as a 2-tuple."""
        # Default Nones.
        for key in ('alternate',):
            if self.data.get(key, ''):
                self.data[key] = [x.strip()
                                  for x in self.data.get(key, '').split(',')]
            else:
                self.data[key] = []
        # Make sure the data representing figures are properly converted.
        for key in NUMERIC_KEYS:
            figure = self.data.get(key, None)
            if figure is None:
                figure = defaultNumericKey(key)
            else:
                figure = float(SI(figure))
            self.data[key] = figure
        # Default integers.
        for key in ('number',):
            self.data.setdefault(key, '0')
        # Default string names.
        for key in STRING_KEYS:
            self.data.setdefault(key, 'unspecified')
        # Sort out the primary.
        if 'primary' in self.data:
            primary = system[self.data['primary']]
            distance = self.data['distance']
        else:
            primary = None
            distance = 0.0
        # If figures are missing, replace them.
        if not self.data['radius']:
            if 'diameter' in self.data:
                diameter = float(SI(self.data['diameter']))
                self.data['radius'] = diameter/2.0
        if not self.data['mass']:
            if 'density' in self.data:
                density = float(SI(self.data['density']))
                self.data['mass'] = FOUR_THIRDS_PI*density*self.radius()**3
        if not self.data['period']:
            if 'synchronous' in self.data:
                self.data['period'] = primary.orbitalPeriodAtDistance(distance)
        if not self.data['normalized moment of inertia']:
            self.data['normalized moment of inertia'] = TWO_FIFTHS
        # Done.
        return primary, distance

    def __getitem__(self, key): return self.data[key]
    def get(self, key, default=None): return self.data.get(key, default)
    def has(self, key): return key in self.data
    def had(self, key): return key in self.rawData

    def isRoot(self): return self.primary() is None
    def secondaries(self): return self.__secondaries[:]

    def select(self, filterFunc):
        return filter(filterFunc, self.__secondaries)

    def planets(self):
        """Get the list of this body's planets."""
        return self.select(lambda x: x.type() == 'planet')

    def asteroids(self):
        """Get the list of this body's asteroids."""
        return self.select(lambda x: x.type() == 'asteroid')

    def comets(self):
        """Get the list of this body's comets."""
        return self.select(lambda x: x.type() == 'comet')

    def satellites(self):
        """Get the list of this body's satellites."""
        return self.select(lambda x: x.type() == 'satellite')

    def rings(self):
        """Get the list of this body's rings."""
        return self.select(lambda x: x.type() == 'ring')

    def name(self): return self.data['name']
    def alternateNames(self): return self.data['alternate']

    def allNames(self):
        names = [self.name()] + self.alternateNames()
        if self.type() == 'asteroid' and self.number():
            names.append(str(self.number()))
        if self.type() == 'comet' and self.number():
            names.append(self.designation())
        return names
    
    def number(self): return int(self.data['number'])
    def type(self): return self.data['type']
    def class_(self): return self.data['class']
    def zone(self): return self.data['zone']
    def designation(self): return self.data['designation']

    def isSynchronous(self):
        """Is the body synchronous with respect to its primary?"""
        return self.has('synchronous')
    
    def hasSurface(self):
        """Does the body have a physical surface?"""
        return not self.has('surfaceless')
    
    def isLuminous(self):
        """Is the body self-luminous?"""
        return self.luminosity() > 0.0
    
    def hasRetrogradeRevolution(self):
        """Does the body orbit retrograde?"""
        return self.has('retrograde')

    def hasRetrogradeRotation(self):
        """Does the body rotate retrograde?"""
        return self.period() < 0.0 or self.obliquity() > PI_OVER_TWO

    def isInclined(self):
        """Is the body highly inclined?"""
        return self.has('inclined') or self.inclination() > PI_OVER_TWELVE

    def isOblique(self):
        """Is the body obliquely tilted?"""
        return self.has('oblique') or self.obliquity() > PI_OVER_TWELVE

    def isOblate(self):
        """Is the body significantly oblate?"""
        return self.has('oblate') or self.ellipticity() > 0.05

    def isNonSpherical(self):
        """Is the body decidedly non-spherical?"""
        return self.has('nonspherical') or self.ellipticity() > 0.05

    def mass(self):
        """The mass of the body. [kg]"""
        return self.data['mass']
    
    def radius(self):
        """The (volumetric mean) radius of the body. [m]"""
        return self.data['radius']
    
    def eccentricity(self):
        """The eccentricity of the body. [dimensionless]"""
        return self.data['eccentricity']

    def inclination(self):
        """The inclination of this body's orbit with the ecliptic. [rad]"""
        return self.data['inclination']
    
    def obliquity(self):
        """The obliquity of this body's axial tilt with respect to its
        orbital plane. [rad]"""
        return self.data['obliquity']
    
    def period(self):
        """The body's rotation period. [s]"""
        return self.data['period']
    
    def temperature(self):
        """The body's surface temperature. [K]"""
        return self.data['temperature']
    
    def luminosity(self):
        """The body's luminosity. [W]"""
        return self.data['luminosity']
    
    def pressure(self):
        """The body's pressure at its surface. [Pa]"""
        return self.data['pressure']
    
    def scaleHeight(self):
        """This body's atmosphere's scale height. [m]"""
        return self.data['scale height']
    
    def opticalDepth(self):
        """The ring system's optical depth. [dimensionless]"""
        return self.data['optical depth']
    
    def albedo(self):
        """The body's visual geometric albedo. [dimensionless]"""
        return self.data['albedo']

    def ellipticity(self):
        """The ellipticity of this body; the ratio of the difference
        in equatorial and polar radii to the equatorial
        radius. [dimensionless]"""
        return self.data['ellipticity']

    def normalizedMomentOfInertia(self):
        """The normalized moment of inertia, that is the coefficient k
        when the moment of inertia is written I = k M R^2.  [dimensionless]"""
        return self.data['normalized moment of inertia']

    def equatorialRadius(self):
        """The equatorial radius of the body. [m]"""
        if self.ellipticity():
            return self.radius()/(1 - self.ellipticity())**ONE_THIRD
        else:
            return self.radius()

    def polarRadius(self):
        """The polar radius of the body. [m]"""
        if self.ellipticity():
            return self.equatorialRadius()*(1 - self.ellipticity())
        else:
            return self.radius()

    def emissivity(self):
        """The emissivity of the body. [dimensionless]"""
        return 1.0 - self.albedo()

    def diameter(self):
        """The body's diameter. [m]"""
        return 2*self.radius()
    
    def circumference(self):
        """The body's circumference. [m]"""
        return TWO_PI*self.radius()
    
    def volume(self):
        """The body's volume. [m^3]"""
        return FOUR_THIRDS_PI*self.radius()**3
    
    def crossSectionalArea(self):
        """The body's cross-sectional area. [m^2]"""
        return PI*self.radius()**2
    
    def surfaceArea(self):
        """The body's total surface area. [m^2]"""
        return FOUR_PI*self.radius()**2
    
    def density(self):
        """The density of this object. [kg/m^3]"""
        return self.mass()/self.volume()

    def periapsis(self):
        """The periapsis of this orbit. [m]"""
        return self.distance()*(1.0 - self.eccentricity())

    def apoapsis(self):
        """The periapsis of this orbit. [m]"""
        return self.distance()*(1.0 + self.eccentricity())

    def orbitalSpeedAroundPrimaryAtDistance(self, distance):
        """Return the instantaneous orbital speed of this body in orbit
        when its speed is at the specified distance. [m/s]"""
        return math.sqrt(TWO_G*self.primary().mass()*(1.0/distance -
                                                      1.0/(2.0*self.distance())))

    def orbitalKineticEnergy(self):
        """The orbital kinetic energy. [J]"""
        return self.mass()*self.massicKineticEnergy()

    def orbitalPotentialEnergy(self):
        """The orbital potential energy. [J]"""
        return self.mass()*self.massicPotentialEnergy()

    def orbitalTotalEnergy(self):
        """The sum of the kinetic and potential energies. [J]"""
        return self.mass()*self.massicTotalEnergy()

    def massEnergy(self):
        """The mass-energy of the object according to Einstein. [J]"""
        return self.mass()*C**2

    def gravitationalBindingEnergy(self):
        """The gravitational binding energy of the object. [J]"""
        return THREE_FIFTHS*G*self.mass()**2/self.radius()

    def orbitalMomentOfInertia(self):
        """The moment of inertia for revolution. [kg m^2]"""
        return self.mass()*self.radius()**2

    def orbitalAngularMomentum(self):
        """Orbital angular momentum around the primary. [kg m^2/s]"""
        return self.mass()*self.massicAngularMomentum()

    def totalAngularMomentum(self):
        """The sum of the rotational and orbital angular momenta. [kg m^2/s]"""
        return self.rotationalAngularMomentum() + self.orbitalAngularMomentum()

    def massRatio(self, barycentric=False):
        """The mass of this body expressed in units of its primary, optionally
        barycentric (expressed in units of the sum of the two).
        [dimensionless]"""
        if barycentric:
            return self.mass()/(self.mass() + self.primary().mass())
        else:
            return self.mass()/self.primary().mass()

    def rotationalMomentOfInertia(self):
        """The rotational moment of inertia. [kg m^2]"""
        return self.normalizedMomentOfInertia()*self.mass()*self.radius()**2

    def rotationSpeed(self):
        """The tangential rotation speed at the equator. [m/s]"""
        return self.rotationalAngularSpeed()*self.radius()

    def rotationalAngularSpeed(self):
        """The rotational angular speed. [rad/s]"""
        return TWO_PI/self.period()

    def rotationalAngularMomentum(self):
        """The rotational angular momentum. [kg m^2/s]"""
        return self.rotationalMomentOfInertia()*self.rotationalAngularSpeed()

    def rotationalKineticEnergy(self):
        """The rotational kinetic energy of the object. [J]"""
        return 0.5*self.rotationalMomentOfInertia()*self.rotationalAngularSpeed()**2

    def orbitalArea(self):
        """The plane area of the body's orbit. [m^2]"""
        return PI*self.distance()**2

    def orbitalAngularSpeedAtDistance(self, distance):
        """The angular speed an object would have at the given distance
        while in orbit around this body. [rad/s]"""
        return math.sqrt(G*self.mass()/distance**3)

    def orbitalSpeedAtDistance(self, distance):
        """The orbital speed an object would have at the given distance
        while in orbit around this body. [m/s]"""
        return self.orbitalAngularSpeedAtDistance(distance)*distance

    def orbitalSpeedAtSurface(self):
        """The orbital speed an object would have to have if just at
        the surface of this body, ignoring atmosphere. [m/s]"""
        return self.orbitalSpeedAtDistance(self.radius())

    def orbitalPeriodAtDistance(self, distance):
        """The orbital period an object would have at the given distance
        while in orbit around this body. [s]"""
        return TWO_PI/self.orbitalAngularSpeedAtDistance(distance)

    def orbitalPeriodAtSurface(self):
        """The orbital period an object would have to have if just at
        the surface of this body, ignoring atmosphere. [s]"""
        return self.orbitalPeriodAtDistance(self.radius())

    def gravityAtDistance(self, distance):
        """The gravitational field strength created by this object
        at the given distance. [m/s^2]"""
        return G*self.mass()/distance**2

    def gravityAtSurface(self):
        """The gravitational field of this body at its surface. [m/s^2]"""
        return self.gravityAtDistance(self.radius())

    def centripetalAccelerationAtSurface(self):
        """The centripetal acceleration a particle has to experience at
        this body's surface to maintain constant distance from its center.
        [m/s^2]"""
        return self.rotationalAngularSpeed()**2*self.radius()

    def netAccelerationAtSurface(self):
        """The net 'weight' (in acceleration) of a particle at rest on the
        surface, taking into account centrifugal pseudoforces. [m/s^2]"""
        return (self.gravityAtSurface() - 
                self.centripetalAccelerationAtSurface())

    def escapeSpeedFromDistance(self, distance):
        """The speed required to escape the gravitational field of this
        object at the specified distance. [m/s]"""
        return math.sqrt(2.0*G*self.mass()/distance)

    def escapeSpeedFromSurface(self):
        """The speed required to escape this body from its surface. [m/s]"""
        return self.escapeSpeedFromDistance(self.radius())

    def escapeSpeedFromPrimary(self):
        """The speed required for this body to escape its primary. [m/s]"""
        return self.primary().escapeSpeedFromDistance(self.distance())

    def angularRadiusFromDistance(self, distance):
        """The angular radius of this object when viewed from the given
        distance. [rad]"""
        return math.atan(self.radius()/distance)

    def angularRadiusOfPrimary(self):
        """The apparent angular radius of the primary as seen from this
        body. [rad]"""
        return self.primary().angularRadiusFromDistance(self.distance())

    def angularRadiusFromPrimary(self):
        """The apparent angular radius of this body as seen from its
        primary. [rad]"""
        return self.angularRadiusFromDistance(self.distance())

    def angularRadiusOfSecondary(self, secondary):
        """The apparent angular radius of the secondary body as seen from
        this body. [rad]"""
        assert secondary in self.__secondaries
        return secondary.angularRadiusFromDistance(secondary.distance())

    def angularSeparationFromPrimaryAtDistance(self, viewingDistance):
        """The angular separate between this body and its primary as
        seen from a distance. [rad]"""
        return math.atan(self.distance()/viewingDistance)

    def solidAngleFromDistance(self, distance):
        """The solid angle subtended by this object when viewed from the
        given distance (for small angles). [sr]"""
        return self.crossSectionalArea()/distance**2

    def solidAngleOfPrimary(self):
        """The solid angle subtended by the primary as viewed from this
        body. [sr]"""
        return self.primary().solidAngleFromDistance(self.distance())

    def solidAngleFromPrimary(self):
        """The solid angle subtended by this body as viewed from its
        primary. [sr]"""
        return self.solidAngleFromDistance(self.distance())

    def solidAngleOfSecondary(self, secondary):
        """The solid angle of the specified secondary body as viewed
        from this body. [sr]"""
        assert secondary in self.__secondaries
        return secondary.solidAngleFromDistance(secondary.distance())

    def distanceToSecondary(self, secondary):
        """The distance to the specified secondary body. [m]"""
        assert secondary in self.__secondaries
        return secondary.distance()

    def barycenterDistanceWithSecondary(self, secondary):
        """Compute the barycenter distance between this body and its
        secondary (measured from the primary). [m]"""
        assert secondary in self.__secondaries
        return secondary.distance()*secondary.mass()/(self.mass() +
                                                      secondary.mass())

    def apostationaryDistance(self):
        """The distance an object orbiting this body would have to be in
        order for its orbital period to equal the rotation period of the
        body itself. [m]"""
        assert self.period() != 0.0
        period = abs(self.period())
        distance = (G*self.mass()*period**2/FOUR_PI_SQUARED)**ONE_THIRD
        assert self.primary() is None or distance < self.distance()
        return distance

    def apostationaryAltitude(self):
        """The apostationary distance expressed as an altitude. [m]"""
        return self.apostationaryDistance() - self.radius()

    def tidalStrengthAtDistance(self, distance):
        """The "tidal strength" of this object due to its gravity at the
        specified distance.  Note that this tidal strength does not
        represent a physical quantity but is instead intended to be used
        to scale relative to other objects. [kg/m^3]"""
        return self.mass()/distance**3

    def tidalStrengthFromSelfAtSurface(self):
        """The "tidal strength" experienced by an object on the surface
        of this body from the body's gravitational field itself.
        [kg/m^3]"""
        return self.tidalStrengthAtDistance(self.radius())

    def tidalStrengthFromPrimary(self):
        """The "tidal strength" on this body induced by its primary.
        [kg/m^3]"""
        return self.primary().tidalStrengthAtDistance(self.distance())

    def tidalStrengthFromSecondary(self, secondary):
        """The "tidal strength" on this body induced by the specified
        secondary body. [kg/m^3]"""
        assert secondary.primary() is self
        return secondary.tidalStrengthAtDistance(secondary.distance())

    def maximumTidalStrengthFromSecondaries(self):
        """The maximum "tidal strength" on this body inducde by all of
        its secondaries. [kg/m^3]"""
        return sum([self.tidalStrengthFromSecondary(x)
                    for x in self.secondaries()])
        
    def maximumTidalStrengthFromSibling(self, sibling):
        """The maximum possible "tidal strength" on this body from the
        specified sibling (in orbit around the same primary). [kg/m^3]"""
        assert self.primary() is sibling.primary()
        inner, outer = self.distance(), sibling.distance()
        if inner > outer:
            inner, outer = outer, inner
        return sibling.tidalStrengthAtDistance(outer - inner)

    def tidalStrengthRangeFromSibling(self, sibling):
        """The range in "tidal strengths" on this body from the
        specified sibling (in orbit around the same primary). [kg/m^3]"""
        assert self.primary() is sibling.primary()
        inner, outer = self.distance(), sibling.distance()
        if inner > outer:
            inner, outer = outer, inner
        return (sibling.tidalStrengthAtDistance(outer - inner),
                sibling.tidalStrengthAtDistance(outer + inner))

    def rocheLimitForDensity(self, density):
        """The Roche limit for this body for an object of the specified
        density. [m]"""
        return self.radius()*(2*self.density()/density)**ONE_THIRD

    def hillRadius(self):
        """The Hill radius, which is the distance at which the tidal strength
        of this body and its primary are equal. [m]"""
        return self.distance()*(ONE_THIRD*self.massRatio())**ONE_THIRD

    def sphereOfInfluence(self):
        """The sphere of influence for this body; for a satellite of this
        body to be relatively stable its orbit needs to be completely
        within the sphere of influence. [m]"""
        return self.distance()*self.massRatio()**0.4

    def tisserandParameter(self, major):
        """The Tisserand parameter for this object with respect to the
        given major planet. [dimensionless]"""
        assert self.primary() is major.primary()
        a_P = major.distance()
        a = self.distance()
        e = self.eccentricity()
        i = self.inclination()
        if i == 0.0:
            cosineI = 1.0
        else:
            cosineI = math.cos(i)
        return a_P/a + 2*cosineI*math.sqrt(a*(1 - e**2)/a_P)

    def areicAtmosphereMass(self):
        """The mass of the atmosphere per unit area. [kg/m^2]"""
        assert self.hasSurface()
        return self.pressure()/self.gravityAtSurface()

    def atmosphereMass(self):
        """The total mass of the atmosphere. [kg]"""
        return self.areicAtmosphereMass()*self.surfaceArea()

    def intensityAtDistance(self, distance):
        """The direct irradiance from this (luminous) object at the given
        distance.  [W/m^2]"""
        return self.luminosity()/(FOUR_PI*distance**2)

    def intensityFromPrimary(self):
        """The direct irradiance from this object's (luminous) parent at
        the distance of this body.  [W/m^2]"""
        return self.primary().intensityAtDistance(self.distance())

    def intensityFromRoot(self):
        """The brightness of this object due to its luminous root parent.
        [W/m^2]"""
        if self.primary().isLuminous():
            return self.intensityFromPrimary()
        else:
            return self.primary().intensityFromRoot()

    def interceptedPower(self):
        """The total rate of intercepted power on this body. [W]"""
        return self.intensityFromRoot()*self.crossSectionalArea()

    def interceptedPowerAtSurface(self):
        """The total rate of intercepted power on this body at its
        surface. [W]"""
        return self.emissivity()*self.intensityFromRoot()*self.crossSectionalArea()

    def reflectedPower(self):
        """The total rate of reflected power on this body.  [W]"""
        return self.albedo()*self.intensityFromRoot()*self.crossSectionalArea()

    def indirectIntensityAtDistance(self, distance, phase=1.0):
        """The indirect intensity from reflected power on this object (which
        is at the given phase) at the given distance.  [W/m^2]"""
        return phase*self.reflectedPower()/(TWO_PI*distance**2)

    def surfaceBrightness(self):
        """The total surface brightness of this object, either by
        self-luminosity or reflection.  [W/m^2]"""
        if self.isLuminous():
            return self.luminosity()/self.surfaceArea()
        else:
            return self.albedo()*self.intensityFromPrimary()/2.0

    def idealTemperature(self, emissivity=None):
        """The equilibrium thermodynamic temperature of this body if it
        were treated as an ideal emitter with the specified emissivity. [K]"""
        if emissivity is None:
            emissivity = self.emissivity()
        denominator = emissivity*SIGMA*self.surfaceArea()
        return (self.interceptedPower()/denominator)**0.25

    def distanceForIntensity(self, intensity):
        """The distance at which one would receive this specified
        intensity from this (self-luminous) body. [m]"""
        assert self.isLuminous()
        return math.sqrt(self.luminosity()/(FOUR_PI*intensity))

    def biozone(self):
        """The range of distances that roughly correspond to the habitable
        zone for this self-luminous object. [m, m]"""
        # Consider the biozone to be the range around a star that ranges plus
        # and minus 20% of Earth's base intensity.
        BASE_INTENSITY = 1400.0
        DISTANCE_FACTOR = 0.20
        minIntensity, maxIntensity = (1/(1 + DISTANCE_FACTOR)**2, 
                                      1/(1 - DISTANCE_FACTOR)**2)
        return (self.distanceForIntensity(maxIntensity*BASE_INTENSITY), 
                self.distanceForIntensity(minIntensity*BASE_INTENSITY))

    def pressureAtAltitude(self, altitude):
        """The pressure at the specified altitude. [Pa]"""
        assert self.scaleHeight() > 0.0
        return self.pressure()*math.exp(-altitude/self.scaleHeight())

    def altitudeAtPressure(self, pressure):
        """The altitude for the specified pressure. [m]"""
        assert self.scaleHeight() > 0.0
        return -self.scaleHeight()*math.log(pressure/self.pressure())

    def distanceAndPhaseAtCentralAngle(self, other, centralAngle):
        """For another sibling body at some central angle relative to this
        body, compute the distance to the body and its phase and return
        them.  [m, 1]"""
        assert self.primary() is other.primary()
        # Use the law of cosines.
        a = self.distance()
        b = other.distance()
        C = centralAngle
        c = math.sqrt(a**2 + b**2 - 2*a*b*math.cos(C))
        intermediate = (-a**2 + b**2 + c**2)/(2*b*c)
        # Due to rounding errors, the cosine can be just slightly > 1 or
        # < -1.  Adjust appropriately.
        if intermediate > 1.0:
            intermediate = 1.0
        if intermediate < -1.0:
            intermediate = -1.0
        A = math.acos(intermediate)
        return c, (-1/PI)*A + 1

    def portionOfOrbitInShadow(self):
        """Calculate the portion (dimensionless) of this satellite's
        orbit during which it is in the shadow of its primary.  This does not
        take into account the motion of the primary around the Sun itself.
        [1]"""
        return (1/PI)*math.atan(self.primary().radius()/self.distance())

    def __str__(self): return self.data['name']

class SurfaceLocation(Location):

    """A location on the surface of a world."""
    
    def __init__(self, world):
        assert isinstance(world, World)
        assert world.hasSurface() ###
        Location.__init__(self, world, world.radius())
        self.__world = world
        self.__orbit = world

    def world(self): return self.__world
    def orbit(self): return self.__orbit

    def atmosphericFactor(self): return 1.0 ###

    def nominalSpeed(self):
        return self.__world.rotationSpeed()

    def escapeSpeedExcess(self):
        """The difference between the escape speed from the object and
        its rotation speed."""
        return (self.escapeSpeedFromPrimary() -
                self.__world.rotationSpeed())

    def __str__(self):
        return "surface of `%s'" % (str(self.primary()))

if False:
    class OrbitalLocation(Location):

        """An orbital location, ostensibly above a world, with a quantized
        set of altitudes."""

        ATMOSPHERIC, LOW, MEDIUM, HIGH = range(4)
        STATIONARY = -1

        ORBITAL_ALTITUDES = (10e3, 200e3, 1000e3, 10000e3) ###
        NAMES = {ATMOSPHERIC: 'atmospheric',
                 LOW: 'low',
                 MEDIUM: 'medium',
                 HIGH: 'high',
                 STATIONARY: 'stationary'}

        def __init__(self, primary, which):
            if which == OrbitalLocation.STATIONARY:
                distance = primary.apostationaryDistance()
            else:
                distance = primary.radius() + self.ORBITAL_ALTITUDES[which]
            Location.__init__(self, primary, distance)
            self.__which = which

        def instability(self):
            raise NotImplementedError ### atmospheric drag should be modelled

        def __str__(self):
            return "%s orbit around `%s'" % \
                   (self.NAMES[self.__which], str(self.primary()))

class PhaseLocation(Location):

    """A location at a phase angle relative to another location (i.e.,
    in the same orbit but some number of degrees ahead or behind."""
    
    def __init__(self, relative, phaseAngle):
        Location.__init__(self, relative.primary(), relative.distance())
        self.__relative = relative
        self.__phaseAngle = phaseAngle
        assert phaseAngle != relative.phaseAngle()

    def relative(self): return self.__relative
    def phaseAngle(self): return self.__phaseAngle

    def angularDisplacement(self):
        """The phase angle, normalized to be in the range (-PI/2, PI/2)."""
        angle = abs(self.__phaseAngle)
        if angle > PI:
            angle = TWO_PI - angle
        return angle

    def displacement(self):
        return self.angularDisplacement()*self.distance()

    def __str__(self):
        return "%s from `%s' around `%s'" % \
               (SI(self.__phaseAngle, 'rad'), 
                str(self.__relative), str(self.primary()))
    
class LibrationLocation(Location):

    """A location representing one of the L1 or L2 points."""
    
    L1, L2 = range(2)
    
    def __init__(self, relative, which):
        Location.__init__(self, relative.primary(), relative.distance())
        self.__relative = relative
        self.__which = which

    def relative(self): return self.__relative
    def which(self): return self.__which

    def instability(self):
        # The e-folding time for L1 and L2 is 1/(2.508 Omega) ~ 2/(5 Omega).
        return 1.0/(2.508*self.orbitalAngularSpeedAroundPrimary())

    def displacement(self):
        """The distance from the corresponding secondary body."""
        massRatio = self.__relative.massRatio()
        return math.pow(ONE_THIRD*massRatio, ONE_THIRD)*self.distance()

    def __str__(self):
        return "%s between `%s' and `%s'" % \
               (('L1', 'L2')[self.__which], 
                str(self.__relative.primary()), str(self.__relative))
    
class LagrangeLocation(PhaseLocation):

    """A location representing one of the remaining three Lagrange
    points: L3, L4, or L5."""
    
    L3, L4, L5 = range(3)
    PHASE_ANGLES = (PI, PI_OVER_THREE, -PI_OVER_THREE)
    
    def __init__(self, relative, which):
        phaseAngle = self.PHASE_ANGLES[which]
        PhaseLocation.__init__(self, relative, phaseAngle)
        self.__which = which

    def which(self): return self.__which

    def instability(self):
        massRatio = self.relative().massRatio(True)
        if self.__which == self.L3:
            ### L3 e-folding time is 1/[Omega (2.625 u)^(1/2)], according to
            # pervect.  This gives an e-folding time for the Sun-Earth L3 of
            # 56 y, whereas references put it at 130 y.
            angularSpeed = self.relative().orbitalAngularSpeedAroundPrimary()
            return 1/(angularSpeed*math.sqrt(2.652*massRatio))
        else: # for L4, L5
            if massRatio < 1.0/25.0:
                return 0.0 # stable
            else:
                raise NotImplementedError # unstable but not sure how

    def __str__(self):
        return "%s between `%s' and `%s'" % \
               (('L3', 'L4', 'L5')[self.__which], 
                str(self.relative().primary()), str(self.relative()))

#
# Orbit
#

class Orbit(ReprMixin):

    """A represention of an orbit, used as a supplementary class."""
    
    def __init__(self, primary, semimajorAxis, eccentricity):
        assert primary is not None
        assert semimajorAxis > 0.0
        assert 0.0 <= eccentricity < 1.0
        self.__primary = primary
        self.__a = semimajorAxis
        self.__e = eccentricity

    def semimajorAxis(self): return self.__a
    def eccentricity(self): return self.__e

    def semiFocalDistance(self): return self.__a*self.__e

    def focusCenterDistance(self): return self.__a*self.__e
    def focusDirectrixDistance(self): return self.__a*(1/self.__e - 1)
    def focalWidth(self): return 2*self.__a*(1 - self.__e**2)

    def semiminorAxis(self):
        return math.sqrt(self.__a**2 - self.focusCenterDistance()**2)

    def massicAngularMomentum(self):
        return self.periapsis()*self.orbitalSpeedAtPeriapsis()
    
    def periapsis(self): return self.__a*(1 - self.__e)
    def apoapsis(self): return self.__a*(1 + self.__e)

    def orbitalSpeedAtPeriapsis(self):
        primaryMass = self.__primary.mass()
        h = (1 + self.__e)/(2*self.__a*(1 - self.__e))
        return math.sqrt(TWO_G*primaryMass*h)
                         
    def orbitalSpeedAtApoapsis(self):
        assert self.__e < 1.0
        primaryMass = self.__primary.mass()
        h = (1 - self.__e)/(2*self.__a*(1 + self.__e))
        return math.sqrt(TWO_G*primaryMass*h)
                         
    def area(self): return PI*self.__a*self.semiminorAxis()

    def isCircular(self): return self.__e == 0.0

    def period(self):
        assert self.__e < 1.0
        primaryMass = self.__primary.mass()
        return math.sqrt(FOUR_PI_SQUARED*self.__a**3/(G*primaryMass))

    def __str__(self):
        return "a = %s, e = %.3f around `%s'" % \
               (SI(self.__a, 'm'), self.__e, str(self.__primary))

#
# Transfer
#

class Transfer(ReprMixin):

    """A representation of a portion of a course identified by a
    single conic section."""
    
    def __init__(self):
        if self.__class__ is Transfer:
            raise NotImplementedError

    def isImpulsive(self):
        """Is this transfer impulsive?"""
        return False
    
    def duration(self):
        """The duration of the transfer. [s]"""
        raise NotImplementedError
    
    def burns(self):
        """A list of burns required to commit the transfer. [m/s]"""
        raise NotImplementedError

    def deltavee(self):
        return sum(self.burns())

    def massicKineticEnergy(self):
        return 0.5*self.deltavee()**2

class ImpulsiveTransfer(Transfer):

    """An impulsive transfer approximation is one where the duration
    is taken to be negligible."""

    def __init__(self):
        if self.__class__ is ImpulsiveTransfer:
            raise NotImplementedError

    def isImpulsive(self): return True

    def duration(self): return 0.0

class VectorTransfer(ImpulsiveTransfer):

    """Simple vector change."""

    def __init__(self, start, finish):
        self.__start = start
        self.__finish = finish

    def start(self): return self.__start
    def finish(self): return self.__finish

    def burns(self):
        startX, startY, startZ = self.__start
        finishx, finishY, finishZ = self.__finish
        deltaX = finishX - startX
        deltaY = finishY - startY
        deltaZ = finishZ - startZ
        return [math.sqrt(deltaX**2 + deltaY**2 + deltaZ**2)]

class CosinesTransfer(ImpulsiveTransfer):

    """Given two speeds and an angle between their vectors, compute the
    deltavee change."""

    def __init__(self, startSpeed, finishSpeed, angle):
        self.__startSpeed = startSpeed
        self.__finishSpeed = finishSpeed
        self.__angle = angle

    def startSpeed(self): return self.__startSpeed
    def finishSpeed(self): return self.__finishSpeed
    def angle(self): return self.__angle
    
    def burns(self):
        a = self.__startSpeed
        b = self.__finishSpeed
        C = self.__angle
        c = math.sqrt(a**2 + b**2 - 2*a*b*math.cos(C))
        return [c]

class PlaneChangeTransfer(ImpulsiveTransfer):

    """Exceute a plane change transfer."""

    def __init__(self, orbitalSpeed, angle):
        self.__orbitalSpeed = orbitalSpeed
        self.__angle = angle

    def orbitalSpeed(self): return self.__orbitalSpeed
    def angle(self): return self.__angle

    def burns(self):
        return [2.0*self.__orbitalSpeed*math.sin(self.__angle/2.0)]
        
class ExtractionTransfer(ImpulsiveTransfer):

    """A transfer out of the theatre corresponding to the given
    location."""
    
    def __init__(self, location):
        ImpulsiveTransfer.__init__(self)
        self.__location = location

    def location(self): return self.__location

    def burns(self): return [self.__location.escapeSpeedExcess()]

    def __str__(self):
        return "out of `%s'" % str(self.__location)
    
class InsertionTransfer(ImpulsiveTransfer):

    """A transfer into the theatre corresponding to the given
    location."""
    
    def __init__(self, location):
        ImpulsiveTransfer.__init__(self)
        self.__location = location

    def location(self): return self.__location

    def burns(self): return [self.__location.escapeSpeedExcess()]

    def __str__(self):
        return "into `%s'" % str(self.__location)

class LaunchTransfer(ImpulsiveTransfer):

    """A transfer launching from the surface of a world."""
    
    def __init__(self, surface):
        ImpulsiveTransfer.__init__(self)
        self.__surface = surface

    def surface(self): return self.__surface
    def location(self): return self.__surface

    def burns(self): return [self.__surface.atmosphericFactor()*
                             self.__surface.escapeSpeedExcess()]

    def __str__(self):
        return "launch from `%s'" % str(self.__surface)

class LandTransfer(ImpulsiveTransfer):

    """A transfer landing on the surface of a world."""
    
    def __init__(self, surface):
        ImpulsiveTransfer.__init__(self)
        self.__surface = surface

    def surface(self): return self.__surface
    def location(self): return self.__surface

    def burns(self): return [self.__surface.atmosphericFactor()*
                             self.__surface.escapeSpeedExcess()]

    def __str__(self):
        return "land on `%s'" % str(self.__surface)

class SendTransfer(ImpulsiveTransfer):

    """A transfer sending an object from a location to an L1 or L2
    libration point."""
    
    def __init__(self, libration):
        Transfer.__init__(self)
        self.__libration = libration

    def libration(self): return self.__libration

    def duration(self): ###
        relative = self.__libration.relative()
        return 0.11*relative.orbitalPeriodAroundPrimary()

    def burns(self): ###
        relative = self.__libration.relative()
        return [0.11*relative.orbitalSpeedAroundPrimary()]

    def __str__(self):
        return "send to `%s'" % str(self.__libration)

class ReturnTransfer(ImpulsiveTransfer):

    """A transfer returning an object to a location from its L1 or L2
    libration point."""
    
    def __init__(self, libration):
        ImpulsiveTransfer.__init__(self)
        self.__libration = libration

    def duration(self): ###
        relative = self.__libration.relative()
        return 0.11*relative.orbitalPeriodAroundPrimary()

    def burns(self): ###
        relative = self.__libration.relative()
        return [0.11*relative.orbitalSpeedAroundPrimary()]

    def __str__(self):
        return "return from `%s'" % str(self.__libration)

class LateralTransfer(Transfer):

    """The lateral or primary transfer that takes place and consists
    of moving between two objects in the same theatre.  For instance,
    in plotting a course from the surface of the Moon to the surface
    of Europa, the lateral transfer would be the transfer from Earth
    to Jupiter in the theatre of the Sun."""
    
    def __init__(self, source, destination):
        if self.__class__ is LateralTransfer:
            raise NotImplementedError
        Transfer.__init__(self)
        assert source.primary() is destination.primary()
        self.__source = source
        self.__destination = destination

    def source(self): return self.__source
    def destination(self): return self.__destination

    def computeOrbitalElements(self, first, second):
        """Compute the semimajor axis and eccentricity, and return
        them as a 2-tuple, required to create an orbit that is
        tangential to the two specified distances."""
        p, q = first, second
        assert p != q
        if p > q:
            p, q = q, p
        a = (p + q)/2.0
        e = (q - p)/(q + p)
        return a, e

class PhaseTransfer(LateralTransfer):

    """A phase transfer as the primary transfer."""
    
    def __init__(self, source, destination):
        LateralTransfer.__init__(self, source, destination)
        assert source.distance() == destination.distance()

    def isImpulsive(self): return True

    def angularDisplacement(self):
        return self.destination().phaseAngle() - self.source().phaseAngle()

    def angularDistance(self):
        angle = abs(self.angularDisplacement())
        if angle > PI:
            angle = TWO_PI - angle
        return angle/TWO_PI

    def duration(self):
        orbitalPeriod = self.source().orbitalPeriodAroundPrimary()
        return self.angularDistance()*orbitalPeriod

    def burns(self):
        orbitalSpeed = self.source().orbitalSpeedAroundPrimary()
        return [self.angularDistance()*orbitalSpeed]

    def __str__(self):
        return "%s around `%s'" % \
               (SI(self.angularDisplacement(), 'rad'), 
                str(self.source().primary()))

class HohmannTransfer(LateralTransfer):

    """A low-energy Hohmann transfer."""
    
    def __init__(self, source, destination):
        LateralTransfer.__init__(self, source, destination)
        a, e = self.computeOrbitalElements(source.distance(), 
                                           destination.distance())
        self.__orbit = Orbit(source.primary(), a, e)

    def isImpulsive(self): return True

    def orbit(self): return self.__orbit
    def duration(self): return self.__orbit.period()/2.0

    def burns(self):
        isReversed = False
        inner, outer = self.source(), self.destination()
        if inner.distance() > outer.distance():
            isReversed = True
            inner, outer = outer, inner
        result = [self.__orbit.orbitalSpeedAtPeriapsis() - 
                  inner.orbitalSpeedAroundPrimary(), 
                  self.__orbit.orbitalSpeedAtApoapsis() - 
                  outer.orbitalSpeedAroundPrimary()]
        if isReversed:
            result.reverse()
        return list(map(abs, result))

    def __str__(self):
        return "Hohmann from `%s' to `%s'" % \
               (str(self.source()), str(self.destination()))

class BiparabolicTransfer(LateralTransfer):

    """A biparabolic transfer, which is only included here for
    completeness since it requires infinite time."""

    ECONOMICAL_DISTANCE_RATIO = 11.94
    
    def __init__(self, source, destination):
        LateralTransfer.__init__(self, source, destination)

    def isImpulsive(self): return True

    def burns(self):
        return [self.source().escapeSpeedExcess(),
                self.destination().escapeSpeedExcess()]

    def duration(self): return INFINITY

    def __str__(self):
        return "biparabolic from `%s' to `%s'" % \
               (str(self.source()), str(self.destination()))

class BiellipticTransfer(LateralTransfer):

    """A bielliptic transfer, which can in some cases take less
    deltavee than a Hohmann transfer but will require much more
    duration."""

    ECONOMICAL_DISTANCE_RATIO = 15.58
    
    def __init__(self, source, destination, reference):
        assert (reference > source.distance() and 
                reference > destination.distance())
        LateralTransfer.__init__(self, source, destination)
        self.__reference = reference
        a1, e1 = self.computeOrbitalElements(source.distance(), reference)
        a2, e2 = self.computeOrbitalElements(reference, destination.distance())
        primary = source.primary()
        self.__orbits = Orbit(primary, a1, e1), Orbit(primary, a2, e2)

    def isImpulsive(self): return True

    def reference(self): return self.__reference
    def orbits(self): return self.__orbits[:]
    def outboundOrbit(self): return self.__orbits[0]
    def inboundOrbit(self): return self.__orbits[1]

    def burns(self):
        start, end = self.source(), self.destination()
        outbound, inbound = self.__orbits
        return [outbound.orbitalSpeedAtPeriapsis() - 
                start.orbitalSpeedAroundPrimary(),
                abs(inbound.orbitalSpeedAtApoapsis() - 
                    outbound.orbitalSpeedAtApoapsis()),
                inbound.orbitalSpeedAtPeriapsis() - 
                end.orbitalSpeedAroundPrimary()]

    def duration(self):
        outbound, inbound = self.__orbits
        return outbound.period()/2.0 + inbound.period()/2.0
    
    def __str__(self):
        return "bielliptic from `%s' to `%s' using %s" % \
               (str(self.source()), str(self.destination()), 
                SI(self.__reference, 'm'))

#
# Maneuver
#

class Maneuver(ReprMixin):

    """A maneuver involves a list of burns (with some total deltavee)
    which ultimately connects several transfers."""

    def __init__(self):
        if self.__class__ is Maneuver:
            raise NotImplementedError

    def transfers(self):
        raise NotImplementedError

    def burns(self):
        raise NotImplementedError

    def deltavee(self):
        return sum(self.burns())

class TransferManeuver(Maneuver):

    """A maneuver merely consisting of those burns associated with a
    transfer."""

    def __init__(self, transfer):
        Maneuver.__init__(self)
        self.__transfer = transfer

    def transfers(self): return [self.__transfer]

    def burns(self): return self.__transfer.burns()

    def __str__(self):
        return "involving `%s'" % (', '.join(map(str, self.transfers())),)

class PartialTransferManeuver(Maneuver):

    """A maneuver merely consisting of one burn associated with a
    transfer."""

    def __init__(self, transfer, whichBurn):
        Maneuver.__init__(self)
        self.__transfer = transfer
        self.__whichBurn = whichBurn

    def whichBurn(self): return self.__whichBurn
    
    def transfers(self): return [self.__transfer]

    def burns(self): return [self.__transfer.burns()[self.__whichBurn]]

    def __str__(self):
        return "involving burn %d of `%s'" % \
               (self.__whichBurn, ', '.join(map(str, self.transfers())))

class OberthManeuver(Maneuver):

    """An Oberth maneuver (possibly chained) is an insertion or
    extraction maneuver which takes into account the benefit one gets
    from doing a burn low in a gravity well."""

    def __init__(self, excessSpeed, innerLocation, outerLocation,
                 transfers=None):
        Maneuver.__init__(self)
        self.__excessSpeed = excessSpeed
        self.__backwards = False
        if not innerLocation.isDescendentOf(outerLocation):
            self.__backwards = True
            innerLocation, outerLocation = outerLocation, innerLocation
            assert innerLocation.isDescendentOf(outerLocation)
        self.__innerLocation = innerLocation
        self.__outerLocation = outerLocation
        if transfers is None:
            transfers = []
        self.__transfers = transfers

    def excessSpeed(self): return self.__excessSpeed
    def innerLocation(self): return self.__innerLocation
    def outerLocation(self): return self.__outerLocation

    def transfers(self): return self.__transfers[:]

    def burns(self):
        hierarchy = []
        location = self.__innerLocation
        while location is not self.__outerLocation:
            assert location is not None
            hierarchy.append(location)
            location = location.primary()
        hierarchy.reverse()
        # Iterate through the different theatres, plugging the deltavee of
        # the last into the excess speed of the next.
        deltavee = self.__excessSpeed
        for location in hierarchy:
            totalSpeed = math.sqrt(deltavee**2 +
                                   location.escapeSpeedFromPrimary()**2)
            deltavee = totalSpeed - location.nominalSpeed()
            location = location.primary()
        return [deltavee]
    
    def __str__(self):
        if self.__backwards:
            return "Oberth from `%s' to `%s' with %s excess" % \
                   (str(self.__outerLocation), str(self.__innerLocation), 
                    SI(self.__excessSpeed, 'm/s'))
        else:
            return "Oberth from `%s' to `%s' with %s excess" % \
                   (str(self.__innerLocation), str(self.__outerLocation), 
                    SI(self.__excessSpeed, 'm/s'))
    

#
# Opportunity
#

class Opportunity(ReprMixin):

    """The representation of an opportunity between two objects in the
    same theatre."""
    
    def __init__(self, source, destination):
        assert source.primary() is destination.primary()
        self.__source = source
        self.__destination = destination

    def source(self): return self.__source
    def destination(self): return self.__destination

    def angle(self):
        """The opportunity angle."""
        return self.__source.opportunityAngleWith(self.__destination)

    def period(self):
        """The rate at which these opportunities will occur."""
        return self.__source.opportunitiesPeriodWith(self.__destination)

    def __str__(self):
        return "opportunity from `%s' to `%s'" % \
               (self.__source, self.__destination)

#
# Course
#

class Course(ReprMixin):

    """The class which plots a series of transfers from a source
    location to a destination location."""
    
    def __init__(self, source, destination, factory=None, additionalArgs=(),
                 useOberthEffect=True):
        self.__source = source
        self.__destination = destination
        self.__opportunities = []
        self.__transfers = []
        self.__maneuvers = []
        self.__lateralTransfer = None
        self.__factory = factory
        self.__additionalArgs = additionalArgs
        self.__useOberthEffect = useOberthEffect
        self.plot()

    def plot(self):
        source, destination = self.__source, self.__destination
        before, after = [], []
        # Handle launches.
        if isinstance(source, SurfaceLocation):
            before.append(LaunchTransfer(source))
            source = source.world()
        # Handle landings.
        if isinstance(destination, SurfaceLocation):
            after.append(LandTransfer(destination))
            destination = destination.world()
        # Handle returns.
        if isinstance(source, LibrationLocation):
            before.append(ReturnTransfer(source))
            source = source.relative()
        # Handle sends.
        if isinstance(destination, LibrationLocation):
            after.append(SendTransfer(destination))
            destination = destination.relative()
        # Arrange into a hierarchy.
        sourcePrimaries = source.primaries()
        sourcePrimaries.reverse()
        destinationPrimaries = destination.primaries()
        destinationPrimaries.reverse()
        theatre = source.findCommonAncestor(destination)
        assert theatre is not None
        sourceIndex = sourcePrimaries.index(theatre)
        destinationIndex = destinationPrimaries.index(theatre)
        # Extractions.
        try:
            transferSource = sourcePrimaries[sourceIndex + 1]
            this = source
            while this is not transferSource:
                before.append(ExtractionTransfer(this))
                this = this.primary()
        except IndexError:
            transferSource = source
        # Insertions.
        try:
            transferDestination = destinationPrimaries[destinationIndex + 1]
            that = destination
            while that is not transferDestination:
                after.append(InsertionTransfer(that))
                that = that.primary()
        except IndexError:
            transferDestination = destination
        after.reverse()
        # Now for the main (lateral) transfer.
        factory = self.__factory
        additionalArgs = self.__additionalArgs
        # If a lateral factory hasn't been provided, dynamically try to
        # choose one.
        if factory is None:
            if transferSource.distance() == transferDestination.distance():
                # If they're in the same orbit, do a phase transfer, provided
                # they have differing phase angles; otherwise do not do any
                # transfer (they're already in the same place).
                if transferSource.phaseAngle() != transferDestination.phaseAngle():
                    factory = PhaseTransfer
            else:
                # Otherwise, do a Hohmann transfer.
                factory = HohmannTransfer
        if factory is not None:
            self.__factory = factory
            args = (transferSource, transferDestination) + additionalArgs
            self.__lateralTransfer = factory(*args)
        # Arrange the final sequence of transfers.
        if self.__lateralTransfer is None:
            middle = []
        else:
            middle = [self.__lateralTransfer]
            # Register the opportunities.
            opportunity = Opportunity(transferSource, transferDestination)
            self.__opportunities.append(opportunity)
        self.__transfers.extend(before + middle + after)
        # Now back up.  If we're going to use the Oberth effect, we'll want
        # to take advantage of it by reformulating the extractions and
        # insertions separately as a single (possibly chained) Oberth
        # maneuver.  To do that, though, we need to have the information on
        # the lateral transfer's burns.
        if (self.__lateralTransfer is not None and
            self.__useOberthEffect and
            self.__lateralTransfer.isImpulsive() and
            len(self.__lateralTransfer.burns()) == 2):
            # The lateral transfer has two impulsive burns, so we can proceed
            # with applying the Oberth effect.
            initialBurn, finalBurn = self.__lateralTransfer.burns()
            if before:
                firstManeuver = OberthManeuver(initialBurn,
                                               before[0].location(),
                                               before[-1].location().primary(),
                                               before[:])
            else:
                firstManeuver = PartialTransferManeuver(self.__lateralTransfer,
                                                        0)
            if after:
                secondManeuver = OberthManeuver(finalBurn,
                                                after[0].location().primary(),
                                                after[-1].location(),
                                                after[:])
            else:
                secondManeuver = PartialTransferManeuver(self.__lateralTransfer,
                                                         1)
            self.__maneuvers = [firstManeuver, secondManeuver]
        else:
            self.__maneuvers = list(map(TransferManeuver, self.__transfers))

    def source(self): return self.__source
    def destination(self): return self.__destination
    def opportunities(self): return self.__opportunities[:]

    def transfers(self): return self.__transfers[:]
    def lateralTransfer(self): return self.__lateralTransfer

    def maneuvers(self): return self.__maneuvers[:]

    def duration(self):
        durations = [x.duration() for x in self.__transfers]
        return sum(durations)
    
    def deltavee(self):
        deltavees = [x.deltavee() for x in self.__maneuvers]
        return sum(deltavees)

    def __str__(self):
        return "from `%s' to `%s'" % \
               (str(self.__source), str(self.__destination))
    
#
# System
#

class System(dict, ReprMixin):

    """The fundamental system class which holds all data on the
    stellar system objects."""
    
    def __init__(self, identifier='<default>', payload=None):
        self.identifier = identifier
        if payload is not None:
            self.update(payload)

    def add(self, name, world):
        """Add a world to the system, along with all its alternate
        names."""
        assert isinstance(world, World)
        assert name == world.name()
        for key in world.allNames():
            assert not key in self, (world, key)
            self[key] = world

    def objects(self):
        """Return all objects in the system."""
        return self.values()

    def __str__(self):
        return "`%s', %d objects" % (self.identifier, len(self))

    def __repr__(self):
        return ReprMixin.__repr__(self)
    

def load(filename=None):
    if filename is None:
        filename = DEFAULT_FILENAME
    try:
        inputFile = gzip.GzipFile(filename + COMPRESSED_EXTENSION, 'rb')
    except IOError:
        inputFile = open(filename, 'rb')
    tag = pickle.load(inputFile)
    logging.debug("Tag: %s" % tag)
    system = pickle.load(inputFile)
    inputFile.close()
    return system

def save(system, tag, filename=DEFAULT_FILENAME, protocol=1, compressed=False):
    if compressed:
        outputFile = gzip.GzipFile(filename + COMPRESSED_EXTENSION, 'wb')
    else:
        outputFile = open(filename, 'wb')
    # There's some Windows bug where you get a "can't import copy_reg" error
    # if you write pickles that contain more than one pickle protocol.
    pickle.dump(tag, outputFile, protocol)
    pickle.dump(system, outputFile, protocol)
    outputFile.close()
