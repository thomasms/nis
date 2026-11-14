"""
OAuth for the HTTP server, so that claude.ai (web, desktop and mobile) can connect. Users log in with GitHub, and
only the GitHub users in the allow list get in.

The server is its own OAuth authorization server (the MCP SDK gives the endpoints: metadata, /register for dynamic
client registration, /authorize, /token and /revoke). For the login it sends the user to GitHub:

    claude.ai -> /authorize -> /oauth/consent -> GitHub login -> /oauth/github/callback -> claude.ai -> /token

The consent page shows the user which client asks for access, and where the code goes. It is necessary because
the server uses one GitHub app for all the clients, and anyone can register a client: GitHub remembers that the
user allowed the app, so without the page a bad client could get a code with no question to the user (the MCP
security best practices require it for this case). The page cannot be in a frame, and its form needs a cookie
that a request from another site does not have.

The tokens are random strings. The state file keeps the registered clients and a hash of each token (not the
tokens), so a restart does not log the users out, and a copy of the file gives no access.

The allow list file has one GitHub login on each line (# for comments). The server reads it again when it changes,
and it checks the list each time a token is used, so a user that you remove loses access at once.
"""

import hashlib
import json
import logging
import os
import secrets
import threading
import time
import urllib.parse
from dataclasses import dataclass
import urllib.request
from pathlib import Path

from mcp.server.auth.provider import (AccessToken, AuthorizationCode, AuthorizationParams, AuthorizeError,
                                      RefreshToken, TokenError, construct_redirect_uri)
from mcp.shared.auth import OAuthClientInformationFull, OAuthToken

log = logging.getLogger("nismcp.oauth")

CALLBACK_PATH = "/oauth/github/callback"
CONSENT_PATH = "/oauth/consent"
CSRF_COOKIE = "nis_consent"


def csrf_cookie_name(state: str) -> str:
    """One cookie for each login, so two logins at the same time do not replace the cookie of the other."""
    return f"{CSRF_COOKIE}_{_hash(state)[:16]}"
USER_PREFIX = "github:"

ACCESS_TOKEN_SECONDS = 3600
REFRESH_TOKEN_SECONDS = 30 * 86400
CODE_SECONDS = 300
LOGIN_SECONDS = 600


def _hash(token: str) -> str:
    return hashlib.sha256(token.encode("utf-8")).hexdigest()


class AllowList:
    """The GitHub logins that can use the server, from a file that is read again when it changes."""

    def __init__(self, path: Path):
        self.path = Path(path)
        self._lock = threading.Lock()
        self._stamp = None
        self._logins = set()
        self._reload()
        if not self._logins:
            log.warning("the allow list %s is empty, so no GitHub user can log in", self.path)

    def _reload(self) -> None:
        stat = self.path.stat()
        stamp = (stat.st_mtime_ns, stat.st_size)
        if stamp == self._stamp:
            return
        logins = set()
        for line in self.path.read_text().splitlines():
            line = line.split("#", 1)[0].strip()
            if line:
                logins.add(line.lower())
        self._logins, self._stamp = logins, stamp
        log.info("%d GitHub users in the allow list %s", len(logins), self.path)

    def allows(self, login: str) -> bool:
        with self._lock:
            try:
                self._reload()
            except OSError as e:
                # keep the list we have, a bad edit must not stop the server
                log.error("cannot read the allow list again: %s", e)
            return login.lower() in self._logins


class GitHub:
    """The calls to GitHub: the login page, the code exchange and the user. The tests replace it."""

    AUTHORIZE_URL = "https://github.com/login/oauth/authorize"
    TOKEN_URL = "https://github.com/login/oauth/access_token"
    USER_URL = "https://api.github.com/user"

    def __init__(self, client_id: str, client_secret: str):
        self.client_id = client_id
        self.client_secret = client_secret

    def login_url(self, state: str, redirect_uri: str) -> str:
        # no scope: the public profile (the login) is enough
        query = urllib.parse.urlencode({"client_id": self.client_id, "redirect_uri": redirect_uri, "state": state,
                                        "allow_signup": "false"})
        return f"{self.AUTHORIZE_URL}?{query}"

    def user_login(self, code: str, redirect_uri: str) -> str:
        """The GitHub login of the user that gave the code. Blocks, so call it in a thread."""
        body = urllib.parse.urlencode({"client_id": self.client_id, "client_secret": self.client_secret,
                                       "code": code, "redirect_uri": redirect_uri}).encode()
        request = urllib.request.Request(self.TOKEN_URL, data=body, headers={"Accept": "application/json"})
        with urllib.request.urlopen(request, timeout=20) as response:
            token = json.load(response)
        if "access_token" not in token:
            raise ValueError(f"GitHub refused the code: {token.get('error_description') or token.get('error')}")
        request = urllib.request.Request(self.USER_URL, headers={
            "Accept": "application/vnd.github+json", "Authorization": f"Bearer {token['access_token']}",
            "User-Agent": "nis-mcp"})
        with urllib.request.urlopen(request, timeout=20) as response:
            return json.load(response)["login"]


@dataclass
class _Login:
    """A login in progress: the request of the client, until the user comes back from GitHub."""
    expires: float
    client_id: str
    params: AuthorizationParams
    csrf: str
    approved: bool = False


class GitHubOAuthProvider:
    """The OAuth authorization server of the MCP SDK (OAuthAuthorizationServerProvider), with the GitHub login."""

    def __init__(self, github: GitHub, allow_list: AllowList, public_url: str, state_file: Path = None,
                 key_store=None):
        self.github = github
        # API keys work too (for Claude Code and scripts), see load_access_token
        self.key_store = key_store
        self.allow_list = allow_list
        self.public_url = public_url.rstrip("/")
        self.resource = f"{self.public_url}/mcp"
        self.callback_url = f"{self.public_url}{CALLBACK_PATH}"
        self.state_file = Path(state_file) if state_file else None
        self._lock = threading.RLock()
        self.clients = {}       # client_id -> OAuthClientInformationFull
        self.logins = {}        # state -> _Login
        self.codes = {}         # code -> AuthorizationCode
        self.access = {}        # hash -> AccessToken (with the hash as the token)
        self.refresh = {}       # hash -> RefreshToken (with the hash as the token)
        self._load()

    # ------------------------------------------------------------------ the state file

    def _load(self) -> None:
        if not self.state_file or not self.state_file.exists():
            return
        state = json.loads(self.state_file.read_text())
        now = time.time()
        self.clients = {c["client_id"]: OAuthClientInformationFull.model_validate(c) for c in state.get("clients", [])}
        self.access = {t["token"]: AccessToken.model_validate(t) for t in state.get("access", [])
                       if t.get("expires_at") is None or t["expires_at"] > now}
        self.refresh = {t["token"]: RefreshToken.model_validate(t) for t in state.get("refresh", [])
                        if t.get("expires_at") is None or t["expires_at"] > now}
        log.info("OAuth state: %d clients, %d refresh tokens", len(self.clients), len(self.refresh))

    def _save(self) -> None:
        if not self.state_file:
            return
        state = {"clients": [c.model_dump(mode="json", exclude_none=True) for c in self.clients.values()],
                 "access": [t.model_dump(mode="json") for t in self.access.values()],
                 "refresh": [t.model_dump(mode="json") for t in self.refresh.values()]}
        temporary = self.state_file.with_suffix(".tmp")
        fd = os.open(temporary, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o600)
        with os.fdopen(fd, "w") as f:
            json.dump(state, f)
        os.replace(temporary, self.state_file)

    def _clean(self) -> None:
        now = time.time()
        self.logins = {k: v for k, v in self.logins.items() if v.expires > now}
        self.codes = {k: v for k, v in self.codes.items() if v.expires_at > now}
        self.access = {k: v for k, v in self.access.items() if v.expires_at is None or v.expires_at > now}
        self.refresh = {k: v for k, v in self.refresh.items() if v.expires_at is None or v.expires_at > now}

    # ------------------------------------------------------------------ clients

    async def get_client(self, client_id: str):
        return self.clients.get(client_id)

    async def register_client(self, client_info: OAuthClientInformationFull) -> None:
        with self._lock:
            self.clients[client_info.client_id] = client_info
            self._save()
        log.info("registered the OAuth client %s (%s)", client_info.client_id, client_info.client_name)

    # ------------------------------------------------------------------ the login

    async def authorize(self, client: OAuthClientInformationFull, params: AuthorizationParams) -> str:
        """Keep the request of the client, and send the user to the consent page."""
        state = secrets.token_urlsafe(32)
        with self._lock:
            self._clean()
            self.logins[state] = _Login(time.time() + LOGIN_SECONDS, client.client_id, params,
                                        csrf=secrets.token_urlsafe(32))
        return f"{self.public_url}{CONSENT_PATH}?{urllib.parse.urlencode({'state': state})}"

    def pending(self, state: str):
        """(the login in progress, its client), or None if the state is not known or too old."""
        with self._lock:
            self._clean()
            entry = self.logins.get(state)
        if entry is None or entry.client_id not in self.clients:
            return None
        return entry, self.clients[entry.client_id]

    def decide(self, state: str, allow: bool) -> str:
        """After the consent page: the GitHub login if the user allows, else back to the client with an error."""
        with self._lock:
            entry = self.logins.get(state)
            if entry is None:
                raise ValueError("the login took too long or is not known, start again from your client")
            if not allow:
                del self.logins[state]
            else:
                entry.approved = True
        if not allow:
            return construct_redirect_uri(str(entry.params.redirect_uri), error="access_denied",
                                          state=entry.params.state, error_description="the user did not allow access")
        return self.github.login_url(state, self.callback_url)

    def complete_login(self, state: str, login: str) -> str:
        """After the GitHub login: the URL that sends the user back to the client, with a code or an error."""
        with self._lock:
            self._clean()
            entry = self.logins.pop(state, None)
        if entry is None or not entry.approved:
            raise ValueError("the login took too long or is not known, start again from your client")
        client_id, params = entry.client_id, entry.params
        if not self.allow_list.allows(login):
            log.warning("refused the GitHub user %s: not in the allow list", login)
            return construct_redirect_uri(str(params.redirect_uri), error="access_denied", state=params.state,
                                          error_description=f"the GitHub user {login} is not allowed on this server")
        code = secrets.token_urlsafe(32)
        with self._lock:
            self.codes[code] = AuthorizationCode(
                code=code, scopes=params.scopes or [], expires_at=time.time() + CODE_SECONDS, client_id=client_id,
                code_challenge=params.code_challenge, redirect_uri=params.redirect_uri,
                redirect_uri_provided_explicitly=params.redirect_uri_provided_explicitly, resource=params.resource,
                subject=USER_PREFIX + login.lower())
        log.info("GitHub user %s logged in", login)
        return construct_redirect_uri(str(params.redirect_uri), code=code, state=params.state)

    # ------------------------------------------------------------------ tokens

    def _check_resource(self, resource):
        # RFC 8707: the token is for this server (the MCP endpoint, or the server URL)
        if resource is not None and resource.rstrip("/") not in (self.resource, self.public_url):
            raise TokenError(error="invalid_target", error_description=f"this server is {self.resource}")

    def _issue(self, client_id: str, scopes: list, resource, subject: str) -> OAuthToken:
        access, refresh = secrets.token_urlsafe(32), secrets.token_urlsafe(32)
        now = int(time.time())
        with self._lock:
            self.access[_hash(access)] = AccessToken(token=_hash(access), client_id=client_id, scopes=scopes,
                                                     expires_at=now + ACCESS_TOKEN_SECONDS, resource=resource,
                                                     subject=subject)
            self.refresh[_hash(refresh)] = RefreshToken(token=_hash(refresh), client_id=client_id, scopes=scopes,
                                                        expires_at=now + REFRESH_TOKEN_SECONDS, resource=resource,
                                                        subject=subject)
            self._clean()
            self._save()
        return OAuthToken(access_token=access, expires_in=ACCESS_TOKEN_SECONDS, refresh_token=refresh,
                          scope=" ".join(scopes) or None)

    async def load_authorization_code(self, client, authorization_code: str):
        code = self.codes.get(authorization_code)
        if code is None or code.client_id != client.client_id or code.expires_at < time.time():
            return None
        return code

    async def exchange_authorization_code(self, client, authorization_code: AuthorizationCode) -> OAuthToken:
        with self._lock:
            if self.codes.pop(authorization_code.code, None) is None:
                raise TokenError(error="invalid_grant", error_description="the code was used already")
        self._check_resource(authorization_code.resource)
        return self._issue(client.client_id, authorization_code.scopes, authorization_code.resource,
                           authorization_code.subject)

    async def load_refresh_token(self, client, refresh_token: str):
        token = self.refresh.get(_hash(refresh_token))
        if token is None or token.client_id != client.client_id:
            return None
        if token.expires_at is not None and token.expires_at < time.time():
            return None
        return token

    async def exchange_refresh_token(self, client, refresh_token: RefreshToken, scopes: list) -> OAuthToken:
        login = (refresh_token.subject or "").removeprefix(USER_PREFIX)
        with self._lock:
            # the token is used one time: a new refresh token comes with the new access token
            self.refresh.pop(refresh_token.token, None)
            self._save()
        if not self.allow_list.allows(login):
            raise TokenError(error="invalid_grant", error_description="the user is not allowed on this server")
        return self._issue(client.client_id, scopes or refresh_token.scopes, refresh_token.resource,
                           refresh_token.subject)

    async def load_access_token(self, token: str):
        """The OAuth access token, or an API key of the keys file. None if neither."""
        access = self.access.get(_hash(token))
        if access is None:
            user = self.key_store.user(token) if self.key_store is not None else None
            if user is None:
                return None
            return AccessToken(token="api-key", client_id="api-key", scopes=[], subject=user)
        if access.expires_at is not None and access.expires_at < time.time():
            return None
        # a user removed from the allow list loses access at once
        if not self.allow_list.allows((access.subject or "").removeprefix(USER_PREFIX)):
            return None
        return access

    async def revoke_token(self, token) -> None:
        with self._lock:
            subject, client_id = token.subject, token.client_id
            # revoke the access and refresh tokens of the same login of the client
            self.access = {k: v for k, v in self.access.items()
                           if not (v.subject == subject and v.client_id == client_id)}
            self.refresh = {k: v for k, v in self.refresh.items()
                            if not (v.subject == subject and v.client_id == client_id)}
            self._save()

    async def exchange_identity_assertion(self, client, params) -> OAuthToken:
        raise TokenError(error="unsupported_grant_type", error_description="not supported")


_PAGE_HEADERS = {
    # the page cannot be in a frame of another site (clickjacking), and it loads nothing
    "X-Frame-Options": "DENY",
    "Content-Security-Policy": "default-src 'none'; style-src 'unsafe-inline'; frame-ancestors 'none'",
    "Cache-Control": "no-store",
    "Referrer-Policy": "no-referrer",
}

_STYLE = ("body{font-family:system-ui,sans-serif;max-width:32rem;margin:3rem auto;padding:0 1rem;color:#0b0b0b;"
          "background:#fcfcfb}code{background:#ebeae6;padding:0 .25rem;border-radius:4px;word-break:break-all}"
          "button{font:inherit;padding:.5rem 1rem;margin-right:.5rem;border-radius:6px;border:1px solid #52514e;"
          "cursor:pointer}.allow{background:#2a78d6;color:#fff;border-color:#2a78d6}.warn{color:#b3261e}"
          "@media (prefers-color-scheme:dark){body{background:#1a1a19;color:#f0efec}code{background:#33322f}}")


def _page(status: int, body: str, headers: dict = None):
    from starlette.responses import HTMLResponse
    html = (f"<!doctype html><html lang=en><meta charset=utf-8><meta name=viewport content='width=device-width'>"
            f"<title>NIS login</title><style>{_STYLE}</style><body>{body}</body></html>")
    return HTMLResponse(html, status_code=status, headers={**_PAGE_HEADERS, **(headers or {})})


def consent_handler(provider: GitHubOAuthProvider):
    """The consent page (GET) and the answer of the user (POST), a Starlette handler for CONSENT_PATH."""
    import html
    from starlette.responses import RedirectResponse

    secure = provider.public_url.startswith("https://")

    async def consent(request):
        if request.method == "GET":
            state = request.query_params.get("state", "")
            found = provider.pending(state)
            if found is None:
                return _page(400, "<p>This login is too old or not known. Start again from your client.</p>")
            entry, client = found
            target = urllib.parse.urlparse(str(entry.params.redirect_uri))
            name = html.escape(client.client_name or "An unnamed client")
            host = html.escape(target.netloc or str(entry.params.redirect_uri))
            known = target.hostname in ("claude.ai", "claude.com") or (target.hostname or "").endswith(
                (".claude.ai", ".claude.com", ".anthropic.com"))
            warning = "" if known else ("<p class=warn>This client does not send you back to claude.ai. Allow it "
                                        "only if you started this login yourself and you know the client.</p>")
            body = (f"<h1>Allow access to NIS?</h1><p><b>{name}</b> asks to use the NIS server as you.</p>"
                    f"<p>After the GitHub login, you go back to <code>{host}</code>.</p>{warning}"
                    f"<form method=post action='{CONSENT_PATH}'>"
                    f"<input type=hidden name=state value='{html.escape(state)}'>"
                    f"<input type=hidden name=csrf value='{html.escape(entry.csrf)}'>"
                    f"<button class=allow name=action value=allow>Log in with GitHub</button>"
                    f"<button name=action value=deny>Cancel</button></form>")
            response = _page(200, body)
            # the form needs this cookie: a page of another site cannot send it (SameSite=Strict)
            response.set_cookie(csrf_cookie_name(state), entry.csrf, max_age=LOGIN_SECONDS, path=CONSENT_PATH, secure=secure,
                                httponly=True, samesite="strict")
            return response

        form = urllib.parse.parse_qs((await request.body()).decode("utf-8", "replace"))
        state = form.get("state", [""])[0]
        csrf = form.get("csrf", [""])[0]
        found = provider.pending(state)
        cookie = request.cookies.get(csrf_cookie_name(state), "")
        if (found is None or not csrf or not secrets.compare_digest(csrf, found[0].csrf)
                or not secrets.compare_digest(cookie, found[0].csrf)):
            return _page(400, "<p>This request is not valid. Start again from your client.</p>")
        try:
            target = provider.decide(state, form.get("action", [""])[0] == "allow")
        except ValueError as e:
            return _page(400, f"<p>{html.escape(str(e))}</p>")
        response = RedirectResponse(target, status_code=303)
        response.delete_cookie(csrf_cookie_name(state), path=CONSENT_PATH)
        return response

    return consent


def callback_handler(provider: GitHubOAuthProvider):
    """The page GitHub sends the user back to after the login (a Starlette handler for CALLBACK_PATH)."""
    import html

    import anyio.to_thread
    from starlette.responses import RedirectResponse

    def page(status: int, text: str):
        return _page(status, f"<p>{html.escape(text)}</p>")

    async def callback(request):
        state = request.query_params.get("state")
        code = request.query_params.get("code")
        if not state or not code:
            error = request.query_params.get("error_description") or request.query_params.get("error") or "no code"
            return page(400, f"The GitHub login did not complete: {error}")
        try:
            login = await anyio.to_thread.run_sync(provider.github.user_login, code, provider.callback_url)
            target = provider.complete_login(state, login)
        except Exception as e:  # noqa: BLE001, the user sees the reason
            log.warning("GitHub login failed: %s", e)
            return page(400, f"The GitHub login failed: {e}")
        return RedirectResponse(target, status_code=302)

    return callback


class ApiKeyHeader:
    """ASGI middleware: X-API-Key: <key> becomes Authorization: Bearer <key>, so the OAuth check takes it too."""

    def __init__(self, app):
        self.app = app

    async def __call__(self, scope, receive, send):
        if scope["type"] == "http":
            headers = scope.get("headers", [])
            names = {k.lower() for k, _ in headers}
            if b"x-api-key" in names and b"authorization" not in names:
                key = next(v for k, v in headers if k.lower() == b"x-api-key")
                scope = dict(scope, headers=list(headers) + [(b"authorization", b"Bearer " + key.strip())])
        return await self.app(scope, receive, send)


def auth_settings(public_url: str):
    from mcp.server.auth.settings import AuthSettings, ClientRegistrationOptions, RevocationOptions
    public_url = public_url.rstrip("/")
    return AuthSettings(issuer_url=public_url, resource_server_url=f"{public_url}/mcp",
                        validate_token_resource=False,
                        client_registration_options=ClientRegistrationOptions(enabled=True),
                        revocation_options=RevocationOptions(enabled=True))
