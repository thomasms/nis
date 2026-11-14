# NIS MCP server

An [MCP](https://modelcontextprotocol.io) server for NIS. With it, an AI assistant (for example Claude) can
irradiate a material with NIS, get the nuclides and plot them. A user can ask, for example:

> If I irradiate 1 kg of U235 with a 14 MeV fusion flux for 1 day, what is the activity 1 day later?

> Plot the top 3 radionuclides that it makes.

The server calls `pynis` directly. It reads the nuclear data one time at the start and keeps it in memory.

## Tools

| Tool | What it does |
|---|---|
| `irradiate` | Irradiates a material (mass fractions of nuclides or natural elements, and a mass in kg, or atoms) in a spectrum (`14MeV`, `thermal` or a group spectrum) with a flux, then cools it. Gives a `run_id`, the total activity and decay heat and the top nuclides at the end of each step, and the assumptions of the run. |
| `get_results` | The nuclides of a run at one time, ranked by activity, decay heat, atoms or mass, with the half-lives. |
| `plot_nuclides` | A PNG plot of a quantity of some nuclides (or the top N) against time. |
| `list_runs` | The runs of the user (the server keeps the last 20 of each user, in memory). |
| `nuclear_data_info` | The status of the server and the nuclear data files. |

Times are text with a unit (`"30 min"`, `"1 d"`, `"10 y"`) or a number of seconds. The flux is in n/cm²/s. If the
user gives no flux, the server uses 1e14 n/cm²/s and says so in the result.

**Limits of the model.** NIS is a 0-D model: the flux is the same in all the material. There is no
self-shielding and no neutron transport. For a large block, the results are an upper limit near the surface. The
tool results say this, so the assistant can tell the user.

## Install

The server needs Python 3.10 or later, and the `mcp` dependency group (the `dev` group has it too):

```bash
uv sync --group mcp
```

It needs the nuclear data (see the main [README](../README.md)). Use the reaction cache (`.nisreac`), not the ENDF
file: it loads in seconds and needs less memory.

```bash
./build/bin/nis_cache --reaction-data data/jeff4.reac --output data/jeff4.nisreac
python scripts/get_fission_yields.py      # data/jeff311.nfy and data/jeff311.sfy
```

Without fission yields, a fission removes the nuclide but makes no fission products. For U235, give them.

## Use it on your computer (stdio, no keys)

The client (Claude Code or Claude Desktop) starts the server as a process. Only you can use it, so there are no
keys. Give the full paths:

```bash
claude mcp add nis -- /path/to/nis/.venv/bin/python -m nismcp serve \
    --decay-data /path/to/nis/data/jeff311.dec --reaction-data /path/to/nis/data/jeff4.nisreac \
    --fission-yields /path/to/nis/data/jeff311.nfy --sf-yields /path/to/nis/data/jeff311.sfy
```

For Claude Desktop, add this to `claude_desktop_config.json`:

```json
{
  "mcpServers": {
    "nis": {
      "command": "/path/to/nis/.venv/bin/python",
      "args": ["-m", "nismcp", "serve"],
      "cwd": "/path/to/nis",
      "env": {
        "NIS_DECAY_DATA": "/path/to/nis/data/jeff311.dec",
        "NIS_REACTION_DATA": "/path/to/nis/data/jeff4.nisreac",
        "NIS_FISSION_YIELDS": "/path/to/nis/data/jeff311.nfy",
        "NIS_SF_YIELDS": "/path/to/nis/data/jeff311.sfy"
      }
    }
  }
}
```

The data loads in the background, so the client connects at once. A tool waits for the data if it is not ready.

## Host it on a server (HTTP, with API keys)

### 1. Make the API keys

Make one key for each user. The command adds a line to the keys file and shows the key **one time**. The file has
only a hash of each key, so a copy of the file does not give access.

```bash
uv run python -m nismcp keygen --name alice --keys-file /etc/nis/keys.txt
```

Send the key to the user on a safe channel. To remove access, delete the line of the user from the file. The
server reads the file again when it changes, so you do not have to restart it.

### 2. Run the server

```bash
uv run python -m nismcp serve --transport http --host 127.0.0.1 --port 8000 \
    --keys-file /etc/nis/keys.txt --allowed-host nis.example.com \
    --decay-data data/jeff311.dec --reaction-data data/jeff4.nisreac \
    --fission-yields data/jeff311.nfy --sf-yields data/jeff311.sfy
```

- The MCP endpoint is `/mcp`. Each request needs the key, in `Authorization: Bearer <key>` or in `X-API-Key: <key>`.
  A request without a valid key gets 401.
- `/health` needs no key. It gives the status of the data (`loading`, `ready` or the error).
- Each user sees only their own runs.
- The server listens on 127.0.0.1. Put a reverse proxy with HTTPS in front of it, so that the keys are not sent
  in clear text. `--allowed-host` is the public name of the server (the proxy sends it in the `Host` header).

A reverse proxy with [Caddy](https://caddyserver.com), which gets the HTTPS certificate itself
(`/etc/caddy/Caddyfile`):

```
nis.example.com {
    reverse_proxy 127.0.0.1:8000
}
```

A systemd service, so that the server starts again after a reboot or a crash (`/etc/systemd/system/nis-mcp.service`):

```ini
[Unit]
Description=NIS MCP server
After=network.target

[Service]
User=nis
WorkingDirectory=/opt/nis
ExecStart=/opt/nis/.venv/bin/python -m nismcp serve --transport http --port 8000 \
    --keys-file /etc/nis/keys.txt --allowed-host nis.example.com \
    --decay-data data/jeff311.dec --reaction-data data/jeff4.nisreac \
    --fission-yields data/jeff311.nfy --sf-yields data/jeff311.sfy
Restart=on-failure

[Install]
WantedBy=multi-user.target
```

### 3. Connect a client

Claude Code:

```bash
claude mcp add --transport http nis https://nis.example.com/mcp --header "Authorization: Bearer <key>"
```

Claude Desktop (its config file starts local servers only, so use the `mcp-remote` bridge, which needs Node.js):

```json
{
  "mcpServers": {
    "nis": {
      "command": "npx",
      "args": ["mcp-remote", "https://nis.example.com/mcp", "--header", "Authorization: Bearer ${NIS_KEY}"],
      "env": { "NIS_KEY": "<key>" }
    }
  }
}
```

The claude.ai web and mobile apps cannot send a fixed key. They need OAuth, which this server does not have yet.

### Size of the server

The memory is set by the nuclear data, not by the number of users. With the reaction cache, a VM with 2 vCPU,
4 GB of memory and 20 GB of disk is enough for a few requests per minute. The server does one run at a time;
other requests wait.

## Tests

```bash
uv run pytest src/tests/python/mcptests.py
```

The tests use a small fake data set, so they need no nuclear data. They test the times, the materials, the API keys,
the runs (against the analytic solution of Co59 activation), the plots, and the tools in-process and over HTTP with
keys. One test irradiates 1 kg of U235 at 14 MeV with the real data, if it is there.
