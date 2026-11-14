"""
The NIS MCP server.

    python -m nismcp serve --decay-data data/jeff311.dec --reaction-data data/jeff4.nisreac \\
        --fission-yields data/jeff311.nfy --sf-yields data/jeff311.sfy              # stdio, for a local client
    python -m nismcp serve --transport http --keys-file keys.txt ...                  # HTTP, with API keys
    python -m nismcp keygen --name alice --keys-file keys.txt                         # a new API key

The data paths can also be given in the environment: NIS_DECAY_DATA, NIS_REACTION_DATA, NIS_FISSION_YIELDS and
NIS_SF_YIELDS.
"""

import argparse
import logging
import os
import sys
from pathlib import Path


def _env_path(name: str):
    value = os.environ.get(name)
    return Path(value) if value else None


def _serve(args) -> int:
    from .data import DataFiles, DataLoader
    from .model import Runner
    from .server import build_server

    if args.decay_data is None:
        print("give --decay-data (or NIS_DECAY_DATA)", file=sys.stderr)
        return 2
    files = DataFiles(decay=args.decay_data, reaction=args.reaction_data, fission_yields=args.fission_yields,
                      sf_yields=args.sf_yields)
    loader = DataLoader(files)
    runner = Runner(loader, max_runs=args.max_runs)

    if args.transport == "stdio":
        server = build_server(loader, runner)
        server.run("stdio")
        return 0

    import uvicorn
    from mcp.server.transport_security import TransportSecuritySettings

    from .auth import ApiKeyMiddleware, KeyStore

    if args.keys_file is None:
        print("the HTTP server needs --keys-file, make a key with: python -m nismcp keygen --name <name> "
              "--keys-file <file>", file=sys.stderr)
        return 2
    store = KeyStore(args.keys_file)
    server = build_server(loader, runner, key_store=store)
    security = None
    if args.allowed_host:
        # behind a reverse proxy the Host header is the public name, so allow it
        hosts = list(args.allowed_host) + ["127.0.0.1:*", "localhost:*"]
        origins = [f"https://{h}" for h in args.allowed_host] + ["http://127.0.0.1:*", "http://localhost:*"]
        security = TransportSecuritySettings(enable_dns_rebinding_protection=True, allowed_hosts=hosts,
                                             allowed_origins=origins)
    app = server.streamable_http_app(host=args.host, transport_security=security)
    # read the data now, so the first request does not wait for it
    loader.start()
    uvicorn.run(ApiKeyMiddleware(app, store), host=args.host, port=args.port, log_level=args.log_level.lower(),
                proxy_headers=True, forwarded_allow_ips="127.0.0.1")
    return 0


def _keygen(args) -> int:
    from .auth import add_key
    try:
        key = add_key(args.keys_file, args.name)
    except ValueError as e:
        print(e, file=sys.stderr)
        return 1
    print(f"API key for {args.name} (it is shown only one time, the file has only its hash):\n\n{key}\n",
          file=sys.stderr)
    return 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="python -m nismcp", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)

    serve = commands.add_parser("serve", help="run the MCP server")
    serve.add_argument("--transport", choices=("stdio", "http"), default="stdio")
    serve.add_argument("--decay-data", "-d", type=Path, default=_env_path("NIS_DECAY_DATA"))
    serve.add_argument("--reaction-data", "-r", type=Path, default=_env_path("NIS_REACTION_DATA"),
                       help="a reaction cache (.nisreac, recommended) or a group-wise ENDF file")
    serve.add_argument("--fission-yields", "-y", type=Path, default=_env_path("NIS_FISSION_YIELDS"))
    serve.add_argument("--sf-yields", type=Path, default=_env_path("NIS_SF_YIELDS"))
    serve.add_argument("--host", default="127.0.0.1", help="HTTP: the address to listen on (default 127.0.0.1, "
                                                           "put a reverse proxy with HTTPS in front)")
    serve.add_argument("--port", type=int, default=8000)
    serve.add_argument("--keys-file", type=Path, default=_env_path("NIS_KEYS_FILE"),
                       help="HTTP: the file of the API keys (see keygen)")
    serve.add_argument("--allowed-host", action="append", default=[],
                       help="HTTP: the public host name of the server (the Host header from the reverse proxy), "
                            "for example nis.example.com; give it more than one time for more names")
    serve.add_argument("--max-runs", type=int, default=20, help="the runs kept in memory for each user")
    serve.add_argument("--log-level", default="INFO", choices=("DEBUG", "INFO", "WARNING", "ERROR"))
    serve.set_defaults(func=_serve)

    keygen = commands.add_parser("keygen", help="make a new API key and add it to the keys file")
    keygen.add_argument("--name", required=True, help="the name of the user, for example alice")
    keygen.add_argument("--keys-file", type=Path, required=True)
    keygen.set_defaults(func=_keygen)

    args = parser.parse_args(argv)
    # the log goes to stderr: over stdio, stdout is the protocol
    logging.basicConfig(level=getattr(args, "log_level", "INFO"), stream=sys.stderr,
                        format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
