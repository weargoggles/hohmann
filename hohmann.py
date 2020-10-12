from starlette.applications import Starlette
from starlette.responses import JSONResponse
from starlette.routing import Route
import botec from botec


async def homepage(request):
    return JSONResponse({'hello': 'world'})


app = Starlette(debug=True, routes=[
    Route('/', homepage),
])
