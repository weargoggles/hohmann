from starlette.applications import Starlette
from starlette.responses import JSONResponse
from starlette.routing import Route
import botec
import humanize
from datetime import timedelta
from starlette.templating import Jinja2Templates
import re
import sentry_sdk
from sentry_sdk.integrations.asgi import SentryAsgiMiddleware

sentry_sdk.init(
    "https://d0b9b59743564012a6e876525f876a7f@o297975.ingest.sentry.io/5461409",
    traces_sample_rate=1.0
)

NUMBERS = re.compile(r'^\d+$')

WORLD = botec.load()
LOCATIONS = sorted(
        loc for loc in WORLD.keys()
        if not NUMBERS.match(loc)
        and WORLD[loc].radius() > 0
    )

templates = Jinja2Templates(directory='templates')

def precisedelta(delta, *args, **kwargs):
    td = timedelta(seconds=delta)
    return humanize.precisedelta(td, *args, **kwargs)

templates.env.filters['precisedelta'] = precisedelta

def loc(request, name):
    origin_loc = WORLD[request.query_params[name]]
    if bool(request.query_params.get(name+'__is_surface', False)) and origin_loc.period() and origin_loc.hasSurface():
        return botec.SurfaceLocation(origin_loc)
    return origin_loc

async def transfer(request):
    origin = request.query_params.get('origin')
    destination = request.query_params.get('destination')
    if origin and destination and origin in LOCATIONS and destination in LOCATIONS and origin != destination:
        origin_loc = loc(request, 'origin')
        destination_loc = loc(request, 'destination')
        course = botec.Course(origin_loc, destination_loc)
    else:
        origin = None
        destination = None
        course = None

    return templates.TemplateResponse('index.html', {
        'request': request,
        'origin': origin,
        'destination': destination,
        'course': course, 'request': request, 'LOCATIONS': LOCATIONS})


app = Starlette(debug=True, routes=[
    Route('/', transfer, name='transfer'),
])

app = SentryAsgiMiddleware(app)
