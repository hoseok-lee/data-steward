"""The catalogue: a SQLite file shared by everyone on the server."""

from __future__ import annotations

import contextlib
import grp
import json
import os
import socket
import sqlite3
from dataclasses import dataclass, field
from pathlib import Path

from . import util

SCHEMA_VERSION = 4

SCHEMA_V1 = """
CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT NOT NULL);
CREATE TABLE datasets (
    id          INTEGER PRIMARY KEY,
    name        TEXT NOT NULL UNIQUE,
    path        TEXT NOT NULL UNIQUE,
    platform    TEXT,
    format      TEXT,
    organ       TEXT,
    species     TEXT,
    study       TEXT,
    description TEXT,
    owner       TEXT,
    steward     TEXT,
    size_bytes  INTEGER,
    n_files     INTEGER,
    added_by    TEXT NOT NULL,
    added_at    TEXT NOT NULL,
    updated_at  TEXT NOT NULL,
    checked_at  TEXT
);
CREATE TABLE dataset_modalities (
    dataset_id INTEGER NOT NULL REFERENCES datasets(id) ON DELETE CASCADE,
    modality   TEXT NOT NULL,
    PRIMARY KEY (dataset_id, modality)
);
CREATE TABLE dataset_tags (
    dataset_id INTEGER NOT NULL REFERENCES datasets(id) ON DELETE CASCADE,
    tag        TEXT NOT NULL,
    PRIMARY KEY (dataset_id, tag)
);
CREATE TABLE links (
    id           INTEGER PRIMARY KEY,
    dataset_id   INTEGER REFERENCES datasets(id) ON DELETE SET NULL,
    dataset_name TEXT NOT NULL,
    user         TEXT NOT NULL,
    link_path    TEXT NOT NULL,
    created_at   TEXT NOT NULL,
    released_at  TEXT
);
CREATE TABLE events (
    id      INTEGER PRIMARY KEY,
    at      TEXT NOT NULL,
    user    TEXT NOT NULL,
    action  TEXT NOT NULL,
    dataset TEXT,
    detail  TEXT
);
"""

# v2: where a dataset came from, and a record of downloads.
MIGRATE_V2 = [
    "ALTER TABLE datasets ADD COLUMN source TEXT",
    "ALTER TABLE datasets ADD COLUMN source_key TEXT",
    "CREATE INDEX datasets_source_key ON datasets(source_key)",
    """CREATE TABLE downloads (
        id          INTEGER PRIMARY KEY,
        source_key  TEXT NOT NULL,
        source      TEXT NOT NULL,
        user        TEXT NOT NULL,
        host        TEXT NOT NULL,
        pid         INTEGER NOT NULL,
        dest        TEXT NOT NULL,
        status      TEXT NOT NULL,
        started_at  TEXT NOT NULL,
        finished_at TEXT,
        dataset_id  INTEGER REFERENCES datasets(id) ON DELETE SET NULL,
        detail      TEXT
    )""",
]

# v3: what is inside a dataset, how datasets relate, file checksums, saved
# queries, and background downloads.
MIGRATE_V3 = [
    "ALTER TABLE datasets ADD COLUMN n_cells INTEGER",
    "ALTER TABLE datasets ADD COLUMN n_features INTEGER",
    "ALTER TABLE datasets ADD COLUMN details TEXT",
    "ALTER TABLE datasets ADD COLUMN inspected_at TEXT",
    "ALTER TABLE datasets ADD COLUMN fingerprint TEXT",
    "ALTER TABLE datasets ADD COLUMN fingerprinted_at TEXT",
    """CREATE TABLE relations (
        child_id  INTEGER NOT NULL REFERENCES datasets(id) ON DELETE CASCADE,
        parent_id INTEGER NOT NULL REFERENCES datasets(id) ON DELETE CASCADE,
        kind      TEXT NOT NULL,
        added_by  TEXT NOT NULL,
        added_at  TEXT NOT NULL,
        PRIMARY KEY (child_id, parent_id)
    )""",
    """CREATE TABLE dataset_files (
        dataset_id INTEGER NOT NULL REFERENCES datasets(id) ON DELETE CASCADE,
        relpath    TEXT NOT NULL,
        size       INTEGER NOT NULL,
        mtime      REAL NOT NULL,
        hash       TEXT NOT NULL,
        PRIMARY KEY (dataset_id, relpath)
    )""",
    "CREATE INDEX dataset_files_hash ON dataset_files(hash, size)",
    """CREATE TABLE saved_queries (
        user       TEXT NOT NULL,
        name       TEXT NOT NULL,
        args       TEXT NOT NULL,
        created_at TEXT NOT NULL,
        PRIMARY KEY (user, name)
    )""",
    "ALTER TABLE downloads ADD COLUMN background INTEGER NOT NULL DEFAULT 0",
    "ALTER TABLE downloads ADD COLUMN notified INTEGER NOT NULL DEFAULT 0",
    "ALTER TABLE downloads ADD COLUMN log TEXT",
    "ALTER TABLE downloads ADD COLUMN total_bytes INTEGER",
]

# v4: datasets vs references (genomes, indices, annotations).
MIGRATE_V4 = [
    "ALTER TABLE datasets ADD COLUMN kind TEXT NOT NULL DEFAULT 'dataset'",
]

# Shared saved queries are stored under this user name.
SHARED = "*"

# "uses" links a dataset to the reference it was processed with.
RELATION_KINDS = ("derived-from", "part-of", "same-as", "uses")

NOW = "datetime('now', 'localtime')"
SEP = "\x1f"

COLS = """d.*,
    (SELECT group_concat(modality, char(31)) FROM dataset_modalities WHERE dataset_id = d.id) AS _modalities,
    (SELECT group_concat(tag, char(31)) FROM dataset_tags WHERE dataset_id = d.id) AS _tags"""

# Columns `set` may change.
EDITABLE = ("name", "platform", "organ", "species", "study", "description", "steward", "source", "source_key",
            "n_cells", "n_features", "details", "inspected_at", "format", "kind")

DATASET_FIELDS = ("name", "path", "platform", "format", "organ", "species", "study", "description",
                  "owner", "steward", "size_bytes", "n_files", "source", "source_key", "n_cells", "n_features",
                  "details", "inspected_at", "kind")


class CatalogError(Exception):
    pass


def _dataset(row: sqlite3.Row) -> dict:
    d = {k: row[k] for k in row.keys() if not k.startswith("_")}
    d["modalities"] = sorted(row["_modalities"].split(SEP)) if row["_modalities"] else []
    d["tags"] = sorted(row["_tags"].split(SEP)) if row["_tags"] else []
    if isinstance(d.get("details"), str):
        d["details"] = json.loads(d["details"])
    return d


@dataclass
class Filter:
    """Values within a field are alternatives (any of several organs);
    modalities and tags must all be present."""
    organs: list[str] = field(default_factory=list)
    species: list[str] = field(default_factory=list)
    modalities: list[str] = field(default_factory=list)
    platforms: list[str] = field(default_factory=list)
    formats: list[str] = field(default_factory=list)
    tags: list[str] = field(default_factory=list)
    owner: str | None = None
    study: str | None = None
    name: str | None = None
    under: str | None = None
    source_key: str | None = None
    source_text: str | None = None
    kinds: list[str] = field(default_factory=list)

    def is_empty(self) -> bool:
        return not any(vars(self).values())


@dataclass
class Changes:
    fields: list[tuple[str, str | None]] = field(default_factory=list)
    set_modalities: list[str] | None = None
    add_modalities: list[str] = field(default_factory=list)
    rm_modalities: list[str] = field(default_factory=list)
    add_tags: list[str] = field(default_factory=list)
    rm_tags: list[str] = field(default_factory=list)

    def is_empty(self) -> bool:
        return not (self.fields or self.set_modalities is not None or self.add_modalities
                    or self.rm_modalities or self.add_tags or self.rm_tags)

    def describe(self) -> str:
        parts = [f"{k}={v if v is not None else '(unset)'}" for k, v in self.fields
                 if k not in ("source_key", "details", "inspected_at")]
        if self.set_modalities is not None:
            parts.append("modalities=" + ",".join(self.set_modalities))
        parts += [f"+modality {m}" for m in self.add_modalities]
        parts += [f"-modality {m}" for m in self.rm_modalities]
        parts += [f"+tag {t}" for t in self.add_tags]
        parts += [f"-tag {t}" for t in self.rm_tags]
        return "; ".join(parts)


class Catalog:
    def __init__(self, path: Path, conn: sqlite3.Connection):
        self.path = path
        self.conn = conn
        self.user = util.current_user()

    @staticmethod
    def _connect(path: Path) -> sqlite3.Connection:
        conn = sqlite3.connect(path, timeout=15, isolation_level=None)
        conn.row_factory = sqlite3.Row
        conn.execute("PRAGMA foreign_keys = ON")
        return conn

    @contextlib.contextmanager
    def tx(self):
        """A write transaction, taken immediately so concurrent writers queue."""
        self.conn.execute("BEGIN IMMEDIATE")
        try:
            yield self.conn
        except BaseException:
            self.conn.execute("ROLLBACK")
            raise
        self.conn.execute("COMMIT")

    def log(self, action: str, dataset: str | None = None, detail: str | None = None) -> None:
        self.conn.execute(f"INSERT INTO events (at, user, action, dataset, detail) VALUES ({NOW}, ?, ?, ?, ?)",
                          (self.user, action, dataset, detail))

    @classmethod
    def create(cls, path: Path, group: str | None = None) -> "Catalog":
        """A new catalogue. With `group`, members of that Unix group can edit it."""
        if path.exists():
            raise CatalogError(f"a catalogue already exists at {path}")
        path.parent.mkdir(parents=True, exist_ok=True)
        cat = cls(path, cls._connect(path))
        with cat.tx() as c:
            for stmt in SCHEMA_V1.split(";"):
                if stmt.strip():
                    c.execute(stmt)
            for stmt in MIGRATE_V2 + MIGRATE_V3 + MIGRATE_V4:
                c.execute(stmt)
            c.execute("INSERT INTO meta VALUES ('schema_version', ?)", (str(SCHEMA_VERSION),))
            cat.log("init", None, group)
        if group:
            try:
                gid = grp.getgrnam(group).gr_gid
            except KeyError:
                raise CatalogError(f"no such group: {group}")
            for p, mode in ((path.parent, 0o2775), (path, 0o664)):
                os.chown(p, -1, gid)
                os.chmod(p, mode)
        return cat

    @classmethod
    def open_readonly(cls, path: Path) -> "Catalog":
        """For the web page: never writes, never upgrades."""
        conn = sqlite3.connect(f"file:{path}?mode=ro", uri=True, timeout=15, isolation_level=None,
                               check_same_thread=False)
        conn.row_factory = sqlite3.Row
        version = int(conn.execute("SELECT value FROM meta WHERE key = 'schema_version'").fetchone()[0])
        if version != SCHEMA_VERSION:
            raise CatalogError(f"the catalogue is at schema {version}; run any stew command once to upgrade it")
        return cls(path, conn)

    @classmethod
    def open(cls, path: Path) -> "Catalog":
        if not path.exists():
            raise CatalogError(
                f"no catalogue at {path}. Create one with `stew init <data folder>`, or point --catalog "
                "(or DATA_STEWARD_CATALOG) at an existing one")
        cat = cls(path, cls._connect(path))
        try:
            row = cat.conn.execute("SELECT value FROM meta WHERE key = 'schema_version'").fetchone()
        except sqlite3.DatabaseError:
            row = None
        if row is None:
            raise CatalogError(f"{path} is not a data-steward catalogue")
        version = int(row[0])
        if version > SCHEMA_VERSION:
            raise CatalogError("the catalogue was written by a newer data-steward; please upgrade")
        if version < SCHEMA_VERSION:
            cat._migrate(version)
        return cat

    def _migrate(self, version: int) -> None:
        try:
            with self.tx() as c:
                if version < 2:
                    for stmt in MIGRATE_V2:
                        c.execute(stmt)
                if version < 3:
                    for stmt in MIGRATE_V3:
                        c.execute(stmt)
                if version < 4:
                    for stmt in MIGRATE_V4:
                        c.execute(stmt)
                c.execute("UPDATE meta SET value = ? WHERE key = 'schema_version'", (str(SCHEMA_VERSION),))
                self.log("upgrade", None, f"schema {version} -> {SCHEMA_VERSION}")
        except sqlite3.OperationalError as e:
            if "readonly" in str(e):
                raise CatalogError(
                    "the catalogue needs a one-time upgrade for this version of stew. "
                    "Ask someone with write access to run any stew command once.")
            raise

    # --- datasets -------------------------------------------------------------

    def insert(self, d: dict) -> int:
        cols = [f for f in DATASET_FIELDS if f in d]
        with self.tx() as c:
            cur = c.execute(
                f"INSERT INTO datasets ({', '.join(cols)}, added_by, added_at, updated_at, checked_at) "
                f"VALUES ({', '.join('?' * len(cols))}, ?, {NOW}, {NOW}, {NOW})",
                [d[k] for k in cols] + [self.user])
            did = cur.lastrowid
            c.executemany("INSERT OR IGNORE INTO dataset_modalities VALUES (?, ?)",
                          [(did, m) for m in d.get("modalities", [])])
            c.executemany("INSERT OR IGNORE INTO dataset_tags VALUES (?, ?)", [(did, t) for t in d.get("tags", [])])
            self.log("add", d["name"], d["path"])
        return did

    def _one(self, cond: str, arg) -> dict | None:
        row = self.conn.execute(f"SELECT {COLS} FROM datasets d WHERE {cond}", (arg,)).fetchone()
        return _dataset(row) if row else None

    def by_id(self, did: int) -> dict | None:
        return self._one("d.id = ?", did)

    def by_path(self, path: Path | str) -> dict | None:
        return self._one("d.path = ?", str(path))

    def names(self) -> list[str]:
        return [r[0] for r in self.conn.execute("SELECT name FROM datasets ORDER BY name")]

    def find(self, key: str) -> dict:
        """A dataset by id, name or path."""
        if key.isdigit():
            d = self.by_id(int(key))
        elif "/" in key:
            p = Path(key)
            d = self.by_path(p.resolve() if p.exists() else p)
        else:
            d = self._one("d.name = ?", key)
        if d:
            return d
        k = key.lower()
        close = sorted((0 if k in n else util.levenshtein(k, n), n) for n in self.names()
                       if k in n or util.levenshtein(k, n) <= max(2, len(k) // 3))
        if close:
            raise CatalogError(f"no dataset '{key}'. Did you mean: {', '.join(n for _, n in close[:5])}?")
        raise CatalogError(f"no dataset '{key}' in the catalogue (see `stew ls`)")

    def name_taken(self, name: str) -> bool:
        return self.conn.execute("SELECT 1 FROM datasets WHERE name = ?", (name,)).fetchone() is not None

    def unique_name(self, base: str, alt: str | None = None) -> str:
        """`base` if free, else `alt`, else base-2, base-3, ..."""
        for n in (base, alt):
            if n and not self.name_taken(n):
                return n
        i = 2
        while self.name_taken(f"{base}-{i}"):
            i += 1
        return f"{base}-{i}"

    def list(self, f: Filter | None = None) -> list[dict]:
        f = f or Filter()
        sql = f"SELECT {COLS} FROM datasets d WHERE 1 = 1"
        args: list = []

        def any_of(col, vals):
            nonlocal sql
            if vals:
                sql += f" AND {col} IN ({', '.join('?' * len(vals))})"
                args.extend(vals)

        any_of("d.organ", f.organs)
        any_of("d.species", f.species)
        any_of("d.platform", f.platforms)
        any_of("d.format", f.formats)
        any_of("d.kind", f.kinds)
        for m in f.modalities:
            sql += " AND EXISTS (SELECT 1 FROM dataset_modalities m WHERE m.dataset_id = d.id AND m.modality = ?)"
            args.append(m)
        for t in f.tags:
            sql += " AND EXISTS (SELECT 1 FROM dataset_tags t WHERE t.dataset_id = d.id AND t.tag = ?)"
            args.append(t)
        if f.owner:
            sql += " AND (d.owner = ? OR d.steward = ?)"
            args += [f.owner, f.owner]
        if f.study:
            sql += " AND d.study = ? COLLATE NOCASE"
            args.append(f.study)
        if f.name:
            sql += " AND d.name LIKE ?"
            args.append(f"%{f.name}%")
        if f.under:
            prefix = f.under.rstrip("/") + "/"
            sql += " AND (d.path = ? OR substr(d.path, 1, length(?)) = ?)"
            args += [f.under, prefix, prefix]
        if f.source_key or f.source_text:
            sql += " AND (d.source_key = ? OR d.source LIKE ?)"
            args += [f.source_key, f"%{f.source_text or f.source_key}%"]
        sql += " ORDER BY d.id"
        return [_dataset(r) for r in self.conn.execute(sql, args)]

    def update(self, d: dict, ch: Changes) -> None:
        with self.tx() as c:
            for col, val in ch.fields:
                assert col in EDITABLE, col
                c.execute(f"UPDATE datasets SET {col} = ? WHERE id = ?", (val, d["id"]))
            if ch.set_modalities is not None:
                c.execute("DELETE FROM dataset_modalities WHERE dataset_id = ?", (d["id"],))
                c.executemany("INSERT OR IGNORE INTO dataset_modalities VALUES (?, ?)",
                              [(d["id"], m) for m in ch.set_modalities])
            c.executemany("INSERT OR IGNORE INTO dataset_modalities VALUES (?, ?)",
                          [(d["id"], m) for m in ch.add_modalities])
            c.executemany("DELETE FROM dataset_modalities WHERE dataset_id = ? AND modality = ?",
                          [(d["id"], m) for m in ch.rm_modalities])
            c.executemany("INSERT OR IGNORE INTO dataset_tags VALUES (?, ?)", [(d["id"], t) for t in ch.add_tags])
            c.executemany("DELETE FROM dataset_tags WHERE dataset_id = ? AND tag = ?",
                          [(d["id"], t) for t in ch.rm_tags])
            c.execute(f"UPDATE datasets SET updated_at = {NOW} WHERE id = ?", (d["id"],))
            self.log("set", d["name"], ch.describe())

    def containing(self, path: Path) -> dict | None:
        """The catalogued dataset at `path`, or the closest one that contains it."""
        p = Path(path)
        for q in [p, *p.parents]:
            d = self.by_path(q)
            if d:
                return d
        return None

    def get_meta(self, key: str) -> str | None:
        r = self.conn.execute("SELECT value FROM meta WHERE key = ?", (key,)).fetchone()
        return r[0] if r else None

    def set_meta(self, key: str, value: str) -> None:
        self.conn.execute("INSERT OR REPLACE INTO meta VALUES (?, ?)", (key, value))

    def move_paths(self, moves: list[tuple[dict, str]], detail: str) -> None:
        """Records new paths after files were moved on disk."""
        with self.tx() as c:
            for d, new in moves:
                c.execute(f"UPDATE datasets SET path = ?, updated_at = {NOW} WHERE id = ?", (new, d["id"]))
            self.log("mv", moves[0][0]["name"], detail)

    def under(self, path: str) -> list[dict]:
        """Datasets at or below a path."""
        return self.list(Filter(under=path))

    def update_stats(self, did: int, owner: str | None, size: int, n_files: int) -> None:
        self.conn.execute(f"UPDATE datasets SET owner = ?, size_bytes = ?, n_files = ?, checked_at = {NOW} "
                          "WHERE id = ?", (owner, size, n_files, did))

    def delete(self, d: dict) -> None:
        with self.tx() as c:
            c.execute("DELETE FROM datasets WHERE id = ?", (d["id"],))
            self.log("rm", d["name"], d["path"])

    # --- links -----------------------------------------------------------------

    def add_link(self, d: dict, link: Path) -> None:
        with self.tx() as c:
            c.execute(f"INSERT INTO links (dataset_id, dataset_name, user, link_path, created_at) "
                      f"VALUES (?, ?, ?, ?, {NOW})", (d["id"], d["name"], self.user, str(link)))
            self.log("get", d["name"], str(link))

    def links(self, user: str | None = None, dataset_id: int | None = None, released: bool = False) -> list[dict]:
        sql, args = "SELECT * FROM links WHERE 1 = 1", []
        if user:
            sql += " AND user = ?"
            args.append(user)
        if dataset_id is not None:
            sql += " AND dataset_id = ?"
            args.append(dataset_id)
        if not released:
            sql += " AND released_at IS NULL"
        return [dict(r) for r in self.conn.execute(sql + " ORDER BY id", args)]

    def active_link(self, link: Path) -> dict | None:
        r = self.conn.execute("SELECT * FROM links WHERE link_path = ? AND released_at IS NULL",
                              (str(link),)).fetchone()
        return dict(r) if r else None

    def release_link(self, link: dict) -> None:
        with self.tx() as c:
            c.execute(f"UPDATE links SET released_at = {NOW} WHERE id = ?", (link["id"],))
            self.log("release", link["dataset_name"], link["link_path"])

    # --- sources and downloads ----------------------------------------------------

    def by_source(self, keys: list[str]) -> list[dict]:
        if not keys:
            return []
        q = f"SELECT {COLS} FROM datasets d WHERE d.source_key IN ({', '.join('?' * len(keys))}) ORDER BY d.id"
        return [_dataset(r) for r in self.conn.execute(q, keys)]

    def mentioning(self, accession: str) -> list[dict]:
        """Datasets whose study, name or path mentions an accession: probably the
        same data, registered before it had a source."""
        a = accession.lower()
        q = (f"SELECT {COLS} FROM datasets d WHERE lower(d.study) = ? OR instr(lower(d.name), ?) > 0 "
             "OR instr(lower(d.path), ?) > 0 ORDER BY d.id")
        return [_dataset(r) for r in self.conn.execute(q, (a, a, a))]

    def downloads(self, status: str | None = None, key: str | None = None, limit: int = 50) -> list[dict]:
        sql, args = "SELECT * FROM downloads WHERE 1 = 1", []
        if status:
            sql += " AND status = ?"
            args.append(status)
        if key:
            sql += " AND source_key = ?"
            args.append(key)
        sql += " ORDER BY id DESC LIMIT ?"
        return [dict(r) for r in self.conn.execute(sql, args + [limit])]

    def active_downloads(self, key: str) -> list[dict]:
        """Downloads of this source still in progress. One whose process has died
        on this machine is marked interrupted instead."""
        alive = []
        for d in self.downloads(status="running", key=key):
            if d["host"] == socket.gethostname() and not _pid_alive(d["pid"]):
                self.finish_download(d["id"], "interrupted")
            else:
                alive.append(d)
        return alive

    def start_download(self, key: str, source: str, dest: Path, check) -> int:
        """Records a download as running. `check()` repeats the duplicate checks
        inside the same transaction, so two people can't start the same one."""
        with self.tx() as c:
            check()
            cur = c.execute(f"INSERT INTO downloads (source_key, source, user, host, pid, dest, status, started_at) "
                            f"VALUES (?, ?, ?, ?, ?, ?, 'running', {NOW})",
                            (key, source, self.user, socket.gethostname(), os.getpid(), str(dest)))
            self.log("download", None, f"started {source} -> {dest}")
            return cur.lastrowid

    def finish_download(self, did: int, status: str, dataset_id: int | None = None, detail: str | None = None):
        self.conn.execute(f"UPDATE downloads SET status = ?, finished_at = {NOW}, dataset_id = ?, detail = ? "
                          "WHERE id = ?", (status, dataset_id, detail, did))

    def set_download_meta(self, did: int, background: bool, log: str | None, total: int | None) -> None:
        self.conn.execute("UPDATE downloads SET background = ?, log = ?, total_bytes = ? WHERE id = ?",
                          (int(background), log, total, did))

    def unnotified(self) -> list[dict]:
        """Background downloads of mine that ended since I last ran stew."""
        rows = [dict(r) for r in self.conn.execute(
            "SELECT * FROM downloads WHERE user = ? AND background = 1 AND notified = 0 "
            "AND status != 'running' ORDER BY id", (self.user,))]
        if rows:
            self.conn.execute(f"UPDATE downloads SET notified = 1 WHERE id IN ({', '.join('?' * len(rows))})",
                              [r["id"] for r in rows])
        return rows

    # --- relations --------------------------------------------------------------------

    def relate(self, child: dict, parent: dict, kind: str) -> None:
        if child["id"] == parent["id"]:
            raise CatalogError("a dataset can't be related to itself")
        if kind not in RELATION_KINDS:
            raise CatalogError(f"kind must be one of {', '.join(RELATION_KINDS)}")
        if kind not in ("same-as", "uses") and child["id"] in {a["id"] for a in self.ancestors(parent)}:
            raise CatalogError(f"{parent['name']} already comes from {child['name']}; that would make a loop")
        with self.tx() as c:
            c.execute(f"INSERT OR REPLACE INTO relations VALUES (?, ?, ?, ?, {NOW})",
                      (child["id"], parent["id"], kind, self.user))
            self.log("relate", child["name"], f"{kind} {parent['name']}")

    def unrelate(self, child: dict, parent: dict) -> bool:
        with self.tx() as c:
            n = c.execute("DELETE FROM relations WHERE (child_id = ? AND parent_id = ?) OR "
                          "(child_id = ? AND parent_id = ?)",
                          (child["id"], parent["id"], parent["id"], child["id"])).rowcount
            if n:
                self.log("unrelate", child["name"], parent["name"])
        return bool(n)

    def relations(self, d: dict) -> tuple[list[tuple[str, dict]], list[tuple[str, dict]]]:
        """(what d comes from, what comes from d), each as (kind, dataset)."""
        up = [(r["kind"], self.by_id(r["parent_id"])) for r in
              self.conn.execute("SELECT * FROM relations WHERE child_id = ? ORDER BY parent_id", (d["id"],))]
        down = [(r["kind"], self.by_id(r["child_id"])) for r in
                self.conn.execute("SELECT * FROM relations WHERE parent_id = ? ORDER BY child_id", (d["id"],))]
        return up, down

    def ancestors(self, d: dict) -> list[dict]:
        seen, todo, out = {d["id"]}, [d["id"]], []
        while todo:
            for r in self.conn.execute("SELECT parent_id FROM relations WHERE child_id = ? "
                                       "AND kind NOT IN ('same-as', 'uses')",
                                       (todo.pop(),)):
                if r[0] not in seen:
                    seen.add(r[0])
                    todo.append(r[0])
                    out.append(self.by_id(r[0]))
        return out

    def all_relations(self) -> list[dict]:
        return [dict(r) for r in self.conn.execute("SELECT child_id, parent_id, kind FROM relations")]

    # --- checksums ----------------------------------------------------------------------

    def files(self, did: int) -> dict[str, dict]:
        return {r["relpath"]: dict(r) for r in
                self.conn.execute("SELECT * FROM dataset_files WHERE dataset_id = ?", (did,))}

    def store_files(self, d: dict, manifest: list[tuple], fingerprint: str) -> None:
        with self.tx() as c:
            c.execute("DELETE FROM dataset_files WHERE dataset_id = ?", (d["id"],))
            c.executemany("INSERT INTO dataset_files VALUES (?, ?, ?, ?, ?)", [(d["id"], *m) for m in manifest])
            c.execute(f"UPDATE datasets SET fingerprint = ?, fingerprinted_at = {NOW} WHERE id = ?",
                      (fingerprint, d["id"]))
            self.log("checksum", d["name"], f"{len(manifest)} file(s), {fingerprint}")

    def duplicate_files(self) -> list[dict]:
        """Files with the same hash and size in more than one dataset."""
        q = """SELECT f.hash, f.size, group_concat(d.name || ':' || f.relpath, char(31)) AS places,
                      count(DISTINCT f.dataset_id) AS n
               FROM dataset_files f JOIN datasets d ON d.id = f.dataset_id
               GROUP BY f.hash, f.size HAVING n > 1 ORDER BY f.size * (n - 1) DESC"""
        return [dict(r) for r in self.conn.execute(q)]

    # --- saved queries ------------------------------------------------------------------

    def save_query(self, name: str, args: list[str], shared: bool) -> None:
        with self.tx() as c:
            c.execute(f"INSERT OR REPLACE INTO saved_queries VALUES (?, ?, ?, {NOW})",
                      (SHARED if shared else self.user, name, json.dumps(args)))

    def saved_query(self, name: str) -> list[str] | None:
        """Mine first, then shared."""
        r = self.conn.execute("SELECT args FROM saved_queries WHERE name = ? AND user IN (?, ?) "
                              "ORDER BY user = ? DESC", (name, self.user, SHARED, SHARED)).fetchone()
        return json.loads(r[0]) if r else None

    def saved_queries(self) -> list[dict]:
        return [dict(r) for r in self.conn.execute(
            "SELECT * FROM saved_queries WHERE user IN (?, ?) ORDER BY user = ?, name", (self.user, SHARED, SHARED))]

    def delete_query(self, name: str, shared: bool) -> bool:
        return self.conn.execute("DELETE FROM saved_queries WHERE user = ? AND name = ?",
                                 (SHARED if shared else self.user, name)).rowcount > 0

    # --- reporting ------------------------------------------------------------------

    def events(self, limit: int = 20) -> list[dict]:
        return [dict(r) for r in self.conn.execute("SELECT * FROM events ORDER BY id DESC LIMIT ?", (limit,))]

    def summary(self) -> dict:
        n, size = self.conn.execute("SELECT count(*), coalesce(sum(size_bytes), 0) FROM datasets").fetchone()

        def counts(sql):
            return [(r[0], r[1]) for r in self.conn.execute(sql)]

        def by(col):
            return counts(f"SELECT coalesce({col}, '(unset)'), count(*) FROM datasets GROUP BY 1 ORDER BY 2 DESC, 1")

        return {
            "datasets": n, "bytes": size,
            "modality": counts("SELECT modality, count(*) FROM dataset_modalities GROUP BY 1 ORDER BY 2 DESC, 1"),
            "organ": by("organ"), "species": by("species"), "platform": by("platform"), "kind": by("kind"),
        }


def _pid_alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True
