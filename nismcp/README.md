# NIS MCP server

An [MCP](https://modelcontextprotocol.io) server for NIS. With it, an AI assistant (for example Claude) can
irradiate a material with NIS, get the nuclides and plot them. A user can ask, for example:

> If I irradiate 1 kg of U235 with a 14 MeV fusion flux for 1 day, what is the activity 1 day later?

> Plot the top 3 radionuclides that it makes.

It runs on your computer (stdio), or on a server over HTTP with API keys and OAuth (GitHub login), so that
claude.ai users can connect. The server calls `pynis` directly. It reads the nuclear data one time at the start and keeps it in memory.

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
keys. The `nismcp` package is not installed in the virtual environment, so Python finds it only from the NIS
folder: `uv run --directory` starts it there, from any folder. Give the full paths:

```bash
claude mcp add nis -- uv run --directory /path/to/nis python -m nismcp serve \
    --decay-data /path/to/nis/data/jeff311.dec --reaction-data /path/to/nis/data/jeff4.nisreac \
    --fission-yields /path/to/nis/data/jeff311.nfy --sf-yields /path/to/nis/data/jeff311.sfy
```

For Claude Desktop, add this to `claude_desktop_config.json` (give the full path of `uv` if Claude Desktop does
not find it, for example `/opt/homebrew/bin/uv`; `which uv` shows it):

```json
{
  "mcpServers": {
    "nis": {
      "command": "uv",
      "args": ["run", "--directory", "/path/to/nis", "python", "-m", "nismcp", "serve"],
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
server reads the file again when it changes, so you do not have to restart it. If the file has an error or is
not there, no key works until you correct it (the server writes the error in its log).

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

The claude.ai web, desktop and mobile apps cannot send a fixed key. For them, turn on OAuth (next section).

## OAuth with the GitHub login (for claude.ai)

With OAuth, a claude.ai user adds the server as a custom connector and logs in with GitHub. Only the GitHub users
in your allow list get in. The server is its own OAuth authorization server (dynamic client registration, PKCE,
access tokens of 1 hour, refresh tokens of 30 days); GitHub only does the login. The API keys still work on the
same server, for Claude Code and scripts.

### 1. Make a GitHub OAuth app

On GitHub: **Settings → Developer settings → OAuth Apps → New OAuth App**:

- **Homepage URL:** `https://nis.example.com`
- **Authorization callback URL:** `https://nis.example.com/oauth/github/callback`

Copy the **Client ID**, and make a **client secret**. The app needs no scopes: the server reads only the GitHub
login of the user.

### 2. Make the allow list

One GitHub login on each line (`#` starts a comment), optionally with the GitHub account id, for example
`/etc/nis/users.txt`:

```
alice
bob 1234567
```

To remove a user, delete the line. The server reads the file again when it changes, and checks it each time a token
is used, so the user loses access at once. If the file has an error or is not there, no one gets in until you
correct it.

A GitHub login can change: a user can rename their account, and then another person can register the old login.
So the server identifies a user by the account id, not by the login. It keeps the id of each login the first time
the user logs in (in the state file), and refuses the login if the id changes. To be sure from the start, give the
id in the allow list (`https://api.github.com/users/<login>` shows it as `id`). If the id of a login changes for a
correct reason, give the new id in the allow list, and remove the login from `github_ids` in the state file.

### 3. Run the server

Give the secret in the environment only (not on the command line, where other users of the machine can see it):

```bash
export NIS_GITHUB_CLIENT_SECRET=<client secret>
uv run python -m nismcp serve --transport http --port 8000 \
    --public-url https://nis.example.com --github-client-id <client id> \
    --allowed-users-file /etc/nis/users.txt --oauth-state-file /var/lib/nis/oauth.json \
    --keys-file /etc/nis/keys.txt \
    --decay-data data/jeff311.dec --reaction-data data/jeff4.nisreac \
    --fission-yields data/jeff311.nfy --sf-yields data/jeff311.sfy
```

- `--public-url` is the HTTPS URL of the server, as the users see it (behind the reverse proxy).
- `--oauth-state-file` keeps the registered clients, the tokens and the account id of each login over a restart,
  so the users do not have to log in again. It has only a hash of each token. Without it, everything is in memory.
- `--keys-file` is optional with OAuth. With it, the API keys work too.

With systemd, put the secret in a file that only root can read, and give it to the service:

```ini
[Service]
# the file has the line NIS_GITHUB_CLIENT_SECRET=... (chmod 600). Keep this comment on its own line: systemd
# reads a # after a value as part of the value.
EnvironmentFile=/etc/nis/oauth.env
```

### 4. Connect claude.ai

In claude.ai: **Settings → Connectors → Add custom connector**. Give the URL `https://nis.example.com/mcp`, and
leave the OAuth client ID and secret empty (the server registers claude.ai itself). Then click **Connect**. The user sees a
page of the NIS server that names the client and the site the login goes back to, and clicks **Log in with
GitHub**. The login must complete in the same browser. GitHub asks the user to log in (and, the first time, to allow the app), and the user comes back to
claude.ai.

The NIS page is necessary for security: the server uses one GitHub app for all the clients, and anyone can register
a client, so without the page a bad client could get access as a user that allowed the app before. The page warns
if the login does not go back to claude.ai. On a Team or Enterprise plan,
an owner adds the connector for the organization, and each user connects with their own GitHub login.

Each user sees only their own runs (the user is `github:<account id>`).

### Size of the server

The memory is mostly the nuclear data. The results of the runs use at most 500 MB for all the users together
(`--max-results-mb`); when they need more, the oldest runs are removed. With the reaction cache, a VM with 2 vCPU,
4 GB of memory and 20 GB of disk is enough for a few requests per minute. The server does one run at a time;
other requests wait.

## Tests

```bash
uv run pytest src/tests/python/mcptests.py
```

The tests use a small fake data set, so they need no nuclear data. They test the times, the materials, the API keys,
the runs (against the analytic solution of Co59 activation), the plots, the tools in-process and over HTTP with
keys, and the OAuth flow as claude.ai does it (registration, PKCE, a fake GitHub login, tokens, refresh, the
allow list, the state file). One test irradiates 1 kg of U235 at 14 MeV with the real data, if it is there.
