"""
API keys for the HTTP server.

The keys file has one line for each user:

    <name> <sha256 of the key, in hex>

Lines that start with # are comments. The file has only the hash of each key, so a copy of the file does not give
access. Make a key with `python -m nismcp keygen --name <name> --keys-file <file>`: it adds the line and shows the key
one time. To remove access, delete the line of the user. The server reads the file again when it changes, so you do
not have to restart it.

A client sends the key in a header, one of:

    Authorization: Bearer <key>
    X-API-Key: <key>
"""

import hashlib
import hmac
import json
import logging
import os
import re
import secrets
import threading
from pathlib import Path

log = logging.getLogger("nismcp.auth")

_NAME = re.compile(r"^[A-Za-z0-9_.@-]+$")
_HASH = re.compile(r"^[0-9a-f]{64}$")
KEY_PREFIX = "nis_"


def hash_key(key: str) -> str:
    return hashlib.sha256(key.encode("utf-8")).hexdigest()


def new_key() -> str:
    """A new random key, 256 bits."""
    return KEY_PREFIX + secrets.token_urlsafe(32)


def parse_keys(text: str, where: str = "keys file") -> dict:
    """{hash: name} from the text of a keys file."""
    keys = {}
    for number, line in enumerate(text.splitlines(), start=1):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        words = line.split()
        if len(words) != 2 or not _NAME.match(words[0]) or not _HASH.match(words[1]):
            raise ValueError(f"{where}, line {number}: give '<name> <sha256 hex>'")
        name, digest = words
        if digest in keys:
            raise ValueError(f"{where}, line {number}: the key of {name} is also the key of {keys[digest]}")
        keys[digest] = name
    return keys


def add_key(keys_file: Path, name: str) -> str:
    """Make a new key for the user, add its hash to the keys file and return the key."""
    if not _NAME.match(name):
        raise ValueError("the name can have only letters, digits and _ . @ -")
    keys_file = Path(keys_file)
    text = keys_file.read_text() if keys_file.exists() else ""
    if name in parse_keys(text, str(keys_file)).values():
        raise ValueError(f"{name} has a key already, delete the line of {name} first to make a new one")
    key = new_key()
    if text and not text.endswith("\n"):
        text += "\n"
    text += f"{name} {hash_key(key)}\n"
    # only the owner can read the file
    fd = os.open(keys_file, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o600)
    with os.fdopen(fd, "w") as f:
        f.write(text)
    return key


class KeyStore:
    """The keys of the keys file. It reads the file again when the file changes."""

    def __init__(self, keys_file: Path):
        self.path = Path(keys_file)
        self._lock = threading.Lock()
        self._stamp = None
        self._keys = {}
        self._reload()
        if not self._keys:
            log.warning("the keys file %s has no keys, so no one can use the server", self.path)

    def _reload(self) -> None:
        stat = self.path.stat()
        stamp = (stat.st_mtime_ns, stat.st_size)
        if stamp == self._stamp:
            return
        keys = parse_keys(self.path.read_text(), str(self.path))
        self._keys, self._stamp = keys, stamp
        log.info("%d API keys read from %s", len(keys), self.path)

    def user(self, key: str):
        """The name of the user of the key, or None if the key is not in the file."""
        with self._lock:
            try:
                self._reload()
            except (OSError, ValueError) as e:
                # fail closed: if the file is deleted or has an error, no key works until it is correct again
                # (keeping the old keys would keep a user that was just removed)
                log.error("cannot read the keys file, so no API key works until it is correct: %s", e)
                self._keys, self._stamp = {}, None
            keys = self._keys
        digest = hash_key(key)
        for known, name in keys.items():
            if hmac.compare_digest(known, digest):
                return name
        return None


def key_from_headers(headers: dict):
    """The key from the Authorization (Bearer) or X-API-Key header, or None."""
    value = headers.get("authorization", "")
    if value[:7].lower() == "bearer ":
        return value[7:].strip() or None
    return headers.get("x-api-key", "").strip() or None


class ApiKeyMiddleware:
    """ASGI middleware: an HTTP request without a known key gets 401. The name of the user goes in
    scope["state"]["nis_user"], so the tools can keep the runs of each user apart."""

    def __init__(self, app, store: KeyStore, open_paths=("/health",)):
        self.app = app
        self.store = store
        self.open_paths = set(open_paths)

    async def __call__(self, scope, receive, send):
        if scope["type"] != "http" or scope.get("path") in self.open_paths:
            return await self.app(scope, receive, send)
        headers = {k.decode("latin-1").lower(): v.decode("latin-1") for k, v in scope.get("headers", [])}
        key = key_from_headers(headers)
        user = self.store.user(key) if key else None
        if user is None:
            client = scope.get("client") or ("?", 0)
            log.warning("refused a request from %s: %s", client[0], "bad key" if key else "no key")
            body = json.dumps({"error": "unauthorized",
                               "detail": "give a valid API key in 'Authorization: Bearer <key>' or 'X-API-Key'"})
            await send({"type": "http.response.start", "status": 401,
                        "headers": [(b"content-type", b"application/json"),
                                    (b"www-authenticate", b'Bearer realm="nis"')]})
            await send({"type": "http.response.body", "body": body.encode()})
            return
        scope.setdefault("state", {})["nis_user"] = user
        return await self.app(scope, receive, send)
