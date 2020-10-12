from starlette.applications import Starlette
from starlette.responses import JSONResponse
from starlette.routing import Route
import botec
import humanize
from datetime import timedelta
from starlette.templating import Jinja2Templates
import re

NUMBERS = re.compile(r'^\d+$')

WORLD = botec.load()
LOCATIONS = sorted(loc for loc in WORLD.keys() if not NUMBERS.match(loc) and WORLD[loc].radius() > 0)

templates = Jinja2Templates(directory='templates')

def precisedelta(delta, *args, **kwargs):
    td = timedelta(seconds=delta)
    return humanize.precisedelta(td, *args, **kwargs)
templates.env.filters['precisedelta'] = precisedelta


def loc(locstring):
    return WORLD[locstring]
    if '+' in locstring:
        body, altitude = locstring.split('+')
        the_body = WORLD[body]
        if the_body.hasSurface() and the_body.radius() > 0:
            return botec.AltitudeLocation(WORLD[body], int(altitude))
        else:
            return botec.Location(WORLD[body], int(altitude))
    loc = WORLD[locstring]
    if loc.hasSurface():
        return botec.SurfaceLocation(WORLD[locstring])
    return botec.AltitudeLocation(loc, 100000)

async def transfer_data(request):
    origin = loc(request.query_params['origin'])
    destination = loc(request.query_params['destination'])
    
    assert origin is not None, destination is not None
    course = botec.Course(origin, destination)
    return {
        'course': str(course),
        'source': str(origin),
        'destination': str(destination),
        'opportunities': [
            {
                'opportunity': str(opportunity),
                'angle': opportunity.angle(),
                'period': humanize.precisedelta(timedelta(seconds=opportunity.period()), minimum_unit="hours"),
            }
            for opportunity in course.opportunities()
        ],
        'maneuvers': [{
            'burns': len(maneuver.burns()),
            'deltavee': maneuver.deltavee(),
            'description': str(maneuver),
        } for maneuver in course.maneuvers()],
        'deltavee': course.deltavee(),
        'duration': course.duration(),
        'humanized_duration': humanize.precisedelta(
            timedelta(seconds=course.duration()),
        ) if course.duration() else 0,
    }

async def transfer(request):
    origin = request.query_params.get('origin')
    destination = request.query_params.get('destination')
    if origin and destination and origin != destination:
        origin = loc(request.query_params['origin'])
        destination = loc(request.query_params['destination'])
        course = botec.Course(origin, destination)
    else:
        origin = None
        destination = None
        course = None

    return templates.TemplateResponse('index.html', {
        'origin': origin,
        'destination': destination,
        'course': course, 'request': request, 'LOCATIONS': LOCATIONS})


app = Starlette(debug=True, routes=[
    Route('/', transfer, name='transfer'),
])
