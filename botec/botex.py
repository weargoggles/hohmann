#!/usr/local/bin/python -i
#
# $Id: //projects/botec/botex.py#20 $ $Date: 2008/09/06 $

"""
Back-of-the-Envelope Calculator executive.
"""

__package__ = 'botec'

import sys

from botec import *


def init(filename):
    return load(filename)

def info(thing, verbose=False):
    if type(thing) is str:
        thing = SYSTEM[thing]
    print "Object: %s" % thing
    if isinstance(thing, Location):
        print " Location:"
        print "  Primary: %s" % thing.primary()
        print "  Distance: %s" % SI(thing.distance(), 'm')
        print "  Phase angle: %s" % SI(thing.phaseAngle(), 'rad')
        if verbose:
            print "  Altitude: %s" % SI(thing.altitude(), 'm')
            print "  Escape speed from primary: %s" % SI(thing.escapeSpeedFromPrimary(), 'm/s')
            print "  Intensity from primary: %s" % SI(thing.intensityFromPrimary(), 'W/m^2')
            print "  Orbital angular speed around primary: %s" % SI(thing.orbitalAngularSpeedAroundPrimary(), 'rad/s')
            print "  Orbital period around primary: %s" % SI(thing.orbitalPeriodAroundPrimary(), 's')
            print "  Orbital speed around primary: %s" % SI(thing.orbitalSpeedAroundPrimary(), 'm/s')
    if isinstance(thing, World):
        print " World:"
        print "  Name: %s" % thing.name()
        print "  Alternate names: %s" % ', '.join(thing.alternateNames())
        print "  Number: %s" % thing.number()
        print "  Type: %s" % thing.type()
        print "  Class: %s" % thing.class_()
        print "  Zone: %s" % thing.zone()
        print "  Designation: %s" % thing.designation()
        print "  Eccentricity: %s" % SI(thing.eccentricity())
        print "  Inclination: %s" % SI(thing.inclination(), 'rad')
        print "  Obliquity: %s" % SI(thing.obliquity(), 'rad')
        print "  Mass: %s" % SI(thing.mass(), 'kg')
        print "  Radius: %s" % SI(thing.radius(), 'm')
        print "  Ellipticity: %s" % SI(thing.ellipticity())
        print "  Period: %s" % SI(thing.period(), 's')
        print "  Temperature: %s" % SI(thing.temperature(), 'K')
        print "  Pressure: %s" % SI(thing.pressure(), 'Pa')
        print "  Scale height: %s" % SI(thing.scaleHeight(), 'm')
        print "  Optical depth: %s" % SI(thing.opticalDepth())
        print "  Albedo: %s" % SI(thing.albedo())
        print "  Luminosity: %s" % SI(thing.luminosity(), 'W')
        if verbose:
            try:
                print "  Apostationary distance: %s" % SI(thing.apostationaryDistance())
            except AssertionError:
                pass
            try:
                print "  Density: %s" % SI(thing.density(), 'kg/m^3')
                print "  Graviational binding energy: %s" % SI(thing.gravitationalBindingEnergy(), 'J')
                print "  Gravity at surface: %s" % SI(thing.gravityAtSurface(), 'm/s^2')
            except ZeroDivisionError:
                pass
            print "  Hill radius: %s" % SI(thing.hillRadius(), 'm')
            print "  Intercepted power: %s" % SI(thing.interceptedPower(), 'W')
            print "  Orbital angular momentum: %s" % SI(thing.orbitalAngularMomentum(), 'kg m^2/s')
            print "  Orbital kinetic energy: %s" % SI(thing.orbitalKineticEnergy(), 'J')
            print "  Orbital potential energy: %s" % SI(thing.orbitalPotentialEnergy(), 'J')
            print "  Orbital moment of inertia: %s" % SI(thing.orbitalMomentOfInertia(), 'kg m^2')
            print "  Reflected power: %s" % SI(thing.reflectedPower(), 'W')
            try:
                print "  Rotational angular momentum: %s" % SI(thing.rotationalAngularMomentum(), 'kg m^2/s')
                print "  Rotational angular speed: %s" % SI(thing.rotationalAngularSpeed(), 'rad/s')
                print "  Rotational kinetic energy: %s" % SI(thing.rotationalKineticEnergy(), 'J')
            except ZeroDivisionError:
                pass
            print "  Rotational moment of inertia: %s" % SI(thing.rotationalMomentOfInertia(), 'kg m^2')
            print "  Angular radius of primary: %s" % SI(thing.angularRadiusOfPrimary(), 'rad')
            print "  Solid angle of primary: %s" % SI(thing.solidAngleOfPrimary(), 'sr')
            print "  Sphere of influence: %s" % SI(thing.sphereOfInfluence(), 'm')

def extended(thing):
    info(thing, True)

def plot(source, destination, factory=None, additionalArgs=(),
         useOberthEffect=True):
    course = Course(source, destination, factory, additionalArgs,
                    useOberthEffect)
    print "Course: %s" % course
    print " Source: %s" % source
    print " Destination: %s" % destination
    if course.opportunities():
        print "Opportunities:"
        for opportunity in course.opportunities():
            print " Opportunity: %s" % opportunity
            print "  Angle: %s" % SI(opportunity.angle(), 'rad')
            print "  Period: %s" % SI(opportunity.period(), 's')
    print " Transfers:"
    for transfer in course.transfers():
        print "  Transfer: %s" % transfer
        #print "   Duration: %s" % SI(transfer.duration(), 's')
        #print "   Nominal burns: %d, totalling %s" % \
        #      (len(transfer.burns()), SI(transfer.deltavee(), 'm/s'))
    print " Maneuvers:"
    for maneuver in course.maneuvers():
        print "  Maneuver: %s" % maneuver
        burns = [SI(x, 'm/s') for x in maneuver.burns()]
        print "   Burns: %s" % (', '.join(map(str, burns)),)
        print "   Involving:"
        for transfer in maneuver.transfers():
            print "    Transfer: %s" % transfer
    print " Duration: %s" % SI(course.duration(), 's')
    print " Deltavee: %s" % SI(course.deltavee(), 'm/s')

# This is the key to the enterprise.  The 'SYSTEM' object is the main entry
# point to all BOTEC functionality.
if len(sys.argv) > 1:
    filename = sys.argv[1]
else:
    filename = None
SYSTEM = init(filename)

if filename is None:
    sun = SYSTEM['Sun']
    mercury = SYSTEM['Mercury']
    venus = SYSTEM['Venus']
    earth = SYSTEM['Earth']
    icarus = SYSTEM['Icarus']
    moon = SYSTEM['Moon']
    mars = SYSTEM['Mars']
    phobos = SYSTEM['Phobos']
    ceres = SYSTEM['Ceres']
    vesta = SYSTEM['Vesta']
    jupiter = SYSTEM['Jupiter']
    io = SYSTEM['Io']
    europa = SYSTEM['Europa']
    saturn = SYSTEM['Saturn']
    iapetus = SYSTEM['Iapetus']
    chiron = SYSTEM['Chiron']
    uranus = SYSTEM['Uranus']
    neptune = SYSTEM['Neptune']
    triton = SYSTEM['Triton']
    pluto = SYSTEM['Pluto']
    charon = SYSTEM['Charon']

    mercurySurface = SurfaceLocation(mercury)
    stationarySunOrbit = StationaryLocation(sun)
    sunMercuryL2 = LibrationLocation(mercury, LibrationLocation.L2)
    earthSurface = SurfaceLocation(earth)
    lowEarthOrbit = AltitudeLocation(earth, 200e3) # LEO is 200 km
    highEarthOrbit = AltitudeLocation(earth, 1400e3) # ICO starts at 1400 km
    stationaryEarthOrbit = StationaryLocation(earth)
    sunEarthL1 = LibrationLocation(earth, LibrationLocation.L1)
    sunEarthL3 = LagrangeLocation(earth, LagrangeLocation.L3)
    sunEarthL4 = LagrangeLocation(earth, LagrangeLocation.L4)
    moonSurface = SurfaceLocation(moon)
    grazingMoonOrbit = GrazingLocation(moon)
    earthMoonL1 = LibrationLocation(moon, LibrationLocation.L1)
    earthMoonL2 = LibrationLocation(moon, LibrationLocation.L2)
    earthMoonL3 = LagrangeLocation(moon, LagrangeLocation.L3)
    earthMoonL4 = LagrangeLocation(moon, LagrangeLocation.L4)
    lowMarsOrbit = AltitudeLocation(mars, 200e3)
    marsSurface = SurfaceLocation(mars)
    lowJupiterOrbit = AltitudeLocation(jupiter, 500e3)
    stationaryJupiterOrbit = StationaryLocation(jupiter)
    grazingJupiterOrbit = GrazingLocation(jupiter)
    sunJupiterL4 = LagrangeLocation(jupiter, LagrangeLocation.L4)
    lowIoOrbit = AltitudeLocation(io, 100e3)
    charonSurface = SurfaceLocation(charon)

    AU = earth.distance()

    GEE = earth.gravityAtSurface()

    MINUTE = 60.0
    HOUR = 60*MINUTE
    DAY = 24*HOUR
    WEEK = 7*DAY
    YEAR = earth.orbitalPeriodAroundPrimary()
    MONTH = YEAR/12
    TURN = 1e6 ###

    LIGHT_MINUTE = C*MINUTE
    LIGHT_HOUR = C*HOUR
    LIGHT_DAY = C*DAY
    LIGHT_WEEK = C*WEEK
    LIGHT_YEAR = C*YEAR
    LIGHT_MONTH = C*MONTH
