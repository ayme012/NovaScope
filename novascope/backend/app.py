"""
NovaScope backend: a small web server (Python standard library only).

It does two jobs:
  1. serves the website in ../frontend  (index.html = star map, timelapse/ = SPHEREx time-lapse)
  2. answers /api/... requests by asking NASA/ESA archives (see spherex.py)

Run (from the project folder):
    python backend/app.py              real data from the archives, opens http://localhost:8000
    python backend/app.py --demo       synthetic data, no internet needed
    python backend/app.py --port 9000 --host 0.0.0.0 --no-browser

API (all GET, all answer JSON except /api/cutout, which answers a FITS file):
    /api/ping                                          {"ok": true, "demo": false}
    /api/search?ra=&dec=&collection=                   {"collection", "count", "rows": [...]}
    /api/cutout?url=&ra=&dec=&size=                    FITS bytes
    /api/gaia?ra=&dec=&r=                              {"epoch": 2016.0, "rows": [...]}
    /api/resolve?name=                                 {"name", "ra", "dec", "pmra", "pmdec", "epoch"}
Errors answer {"error": "<what went wrong>"} with HTTP 400 (bad request) or 502 (archive problem).

To add an endpoint: write a function that takes the query dict and returns a dict (or bytes),
then add it to ROUTES below.
"""
import argparse
import json
import os
import socket
import sys
import threading
import urllib.error
import urllib.parse
import webbrowser
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer

import config
import spherex


# ----------------------------------------------------------------------------- endpoints

def _float(q, key, default=None):
    if key not in q:
        if default is None:
            raise ValueError(f"Missing parameter '{key}'")
        return default
    try:
        return float(q[key])
    except ValueError:
        raise ValueError(f"Parameter '{key}' must be a number, got '{q[key]}'")


def api_ping(q):
    return {"ok": True, "demo": config.DEMO}


def api_search(q):
    return spherex.search_images(_float(q, "ra"), _float(q, "dec"), q.get("collection"))


def api_cutout(q):
    if "url" not in q:
        raise ValueError("Missing parameter 'url'")
    return spherex.get_cutout(q["url"], _float(q, "ra"), _float(q, "dec"), _float(q, "size", 0.1))


def api_gaia(q):
    return spherex.gaia_stars(_float(q, "ra"), _float(q, "dec"), _float(q, "r", 0.1))


def api_resolve(q):
    if not q.get("name"):
        raise ValueError("Missing parameter 'name'")
    return spherex.resolve_name(q["name"])


ROUTES = {
    "/api/ping": api_ping,
    "/api/search": api_search,
    "/api/cutout": api_cutout,
    "/api/gaia": api_gaia,
    "/api/resolve": api_resolve,
}


# ----------------------------------------------------------------------------- web server

class Handler(SimpleHTTPRequestHandler):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, directory=config.FRONTEND_DIR, **kwargs)

    def log_message(self, fmt, *args):
        if args and "/api/" in str(args[0]):           # only log API calls, not every image or script
            sys.stderr.write("  " + (fmt % args) + "\n")

    def end_headers(self):
        # always send the newest page files, so an update shows up without clearing the browser's cache
        if not self.path.startswith("/api/"):
            self.send_header("Cache-Control", "no-cache")
        super().end_headers()

    def reply(self, code, body, content_type):
        self.send_response(code)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Access-Control-Allow-Origin", "*")
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        parsed = urllib.parse.urlparse(self.path)
        route = ROUTES.get(parsed.path)
        if route is None:
            return super().do_GET()                     # a frontend file
        query = {k: v[0] for k, v in urllib.parse.parse_qs(parsed.query).items()}
        try:
            result = route(query)
            if isinstance(result, bytes):
                return self.reply(200, result, "application/fits")
            return self.reply(200, json.dumps(result).encode(), "application/json")
        except (ValueError, LookupError) as e:
            code, message = 400, str(e)
        except urllib.error.HTTPError as e:
            code, message = 502, f"The archive answered HTTP {e.code}. Try again in a minute."
        except (urllib.error.URLError, TimeoutError, OSError) as e:
            code, message = 502, f"Could not reach the archive ({getattr(e, 'reason', e)}). Check your internet connection."
        except Exception as e:                          # keep the server running whatever happens
            code, message = 500, f"{type(e).__name__}: {e}"
        self.reply(code, json.dumps({"error": message}).encode(), "application/json")


VERSION = "2026-11-15 (blue theme)"      # shown when the server starts, so you can tell which copy is running


class Server(ThreadingHTTPServer):
    """On Windows, SO_REUSEADDR lets a second server share a busy port without any error, and the browser
    keeps talking to the old one. Ask for the port exclusively instead, so a busy port is reported."""
    allow_reuse_address = os.name != "nt"

    def server_bind(self):
        if os.name == "nt" and hasattr(socket, "SO_EXCLUSIVEADDRUSE"):
            self.socket.setsockopt(socket.SOL_SOCKET, socket.SO_EXCLUSIVEADDRUSE, 1)
        super().server_bind()


def open_server(host, port):
    """Use the asked-for port; if another program (often an older NovaScope) holds it, say so and try the next ones."""
    for p in range(port, port + 10):
        try:
            return Server((host, p), Handler), p
        except OSError:
            print(f"Port {p} is already in use, probably by another NovaScope that is still running "
                  f"(an older copy?). Close its black window, or end 'Python' in Task Manager.")
    raise SystemExit("No free port found between %d and %d." % (port, port + 9))


def main():
    ap = argparse.ArgumentParser(description="NovaScope backend")
    ap.add_argument("--host", default=config.HOST)
    ap.add_argument("--port", type=int, default=config.PORT)
    ap.add_argument("--demo", action="store_true", help="synthetic data, no internet needed")
    ap.add_argument("--no-browser", action="store_true")
    args = ap.parse_args()
    if args.demo:
        config.DEMO = True

    server, port = open_server(args.host, args.port)
    url = f"http://localhost:{port}/"
    if port != args.port:
        print(f"Using port {port} instead. Open {url} (not :{args.port}, which shows the other copy).")
    print(f"NovaScope {VERSION} running at {url}" + ("  [DEMO MODE: synthetic data]" if config.DEMO else ""))
    print("Press Ctrl+C to stop.")
    if config.OPEN_BROWSER and not args.no_browser and args.host in ("127.0.0.1", "localhost"):
        threading.Timer(0.8, lambda: webbrowser.open(url)).start()
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nStopped.")


if __name__ == "__main__":
    main()
