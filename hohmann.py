from starlette.applications import Starlette
from starlette.responses import JSONResponse
from starlette.routing import Route
import botec.botec as botec


async def homepage(request):
    return JSONResponse({'hello': 'world'})


app = Starlette(debug=True, routes=[
    Route('/', homepage),
])
