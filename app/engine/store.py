"""PostgreSQL-style relational schema on SQLite, so the workbench runs anywhere with no server."""
import sqlite3, threading, time, random
import pandas as pd

SCHEMA = """
CREATE TABLE samples(sample_id INTEGER PRIMARY KEY, name TEXT UNIQUE NOT NULL);
CREATE TABLE sequences(seq_id INTEGER PRIMARY KEY, sample_id INTEGER NOT NULL REFERENCES samples(sample_id),
  name TEXT NOT NULL, kind TEXT NOT NULL CHECK(kind IN ('DNA','PROTEIN')),
  length INTEGER NOT NULL CHECK(length > 0), gc REAL, residues TEXT NOT NULL);
CREATE TABLE jobs(job_id INTEGER PRIMARY KEY, run_id INTEGER NOT NULL, task TEXT NOT NULL,
  seq_id INTEGER REFERENCES sequences(seq_id), priority INTEGER NOT NULL DEFAULT 0,
  status TEXT NOT NULL DEFAULT 'pending' CHECK(status IN ('pending','running','done','failed')));
CREATE TABLE executions(exec_id INTEGER PRIMARY KEY, job_id INTEGER NOT NULL REFERENCES jobs(job_id),
  worker INTEGER, started_ms REAL, finished_ms REAL, cpu_ms REAL);
CREATE TABLE results(result_id INTEGER PRIMARY KEY, job_id INTEGER NOT NULL REFERENCES jobs(job_id),
  seq_id INTEGER REFERENCES sequences(seq_id), position INTEGER, score REAL, detail TEXT);
"""
INDEXES = {
    'idx_seq_sample': "CREATE INDEX idx_seq_sample ON sequences(sample_id)",
    'idx_seq_length': "CREATE INDEX idx_seq_length ON sequences(length)",
    'idx_jobs_run': "CREATE INDEX idx_jobs_run ON jobs(run_id)",
    'idx_exec_job': "CREATE INDEX idx_exec_job ON executions(job_id)",
    'idx_res_seq': "CREATE INDEX idx_res_seq ON results(seq_id)"}

class Store:
    def __init__(self):
        self.lock = threading.Lock(); self._open()
    def _open(self):
        self.con = sqlite3.connect(':memory:', check_same_thread=False)
        self.con.execute("PRAGMA foreign_keys=ON"); self.con.executescript(SCHEMA); self.con.commit()
        self.run_id = 0; self.indexed = False
    def reset(self):
        with self.lock: self.con.close(); self._open()
    def add_sequences(self, sample, seqs, kind_of, gc):
        with self.lock, self.con:
            cur = self.con.execute("INSERT INTO samples(name) VALUES(?)", (sample,)); sid = cur.lastrowid
            self.con.executemany("INSERT INTO sequences(sample_id,name,kind,length,gc,residues) VALUES(?,?,?,?,?,?)",
                                 [(sid, n, kind_of(s), len(s), gc(s), s) for n, s in seqs.items()])
    def seq_ids(self):
        with self.lock: return dict(self.con.execute("SELECT name, seq_id FROM sequences"))
    def df(self, sql, params=()):
        with self.lock: return pd.read_sql_query(sql, self.con, params=params)
    def explain(self, sql):
        with self.lock: return [r[3] for r in self.con.execute("EXPLAIN QUERY PLAN " + sql)]
    def timed(self, sql, reps=20):
        with self.lock:
            best = 1e9
            for _ in range(reps):
                t = time.perf_counter(); self.con.execute(sql).fetchall(); best = min(best, time.perf_counter() - t)
        return best * 1000
    def set_indexes(self, on):
        with self.lock, self.con:
            for n, ddl in INDEXES.items():
                self.con.execute(f"DROP INDEX IF EXISTS {n}")
                if on: self.con.execute(ddl)
            self.indexed = on
    def record_run(self, runs, seq_of):
        """Persist one workflow run in a single transaction: jobs, their executions and their results."""
        with self.lock, self.con:
            self.run_id += 1; ids = {}
            for r in runs:
                sid = seq_of.get(r['id'].split(':', 1)[-1])
                c = self.con.execute("INSERT INTO jobs(run_id,task,seq_id,priority,status) VALUES(?,?,?,?,'done')", (self.run_id, r['id'], sid, r['priority']))
                ids[r['id']] = (c.lastrowid, sid)
                self.con.execute("INSERT INTO executions(job_id,worker,started_ms,finished_ms,cpu_ms) VALUES(?,?,?,?,?)", (c.lastrowid, r['lane'], r['start'], r['end'], r['cpu_ms']))
                res = r['result'] or {}
                if 'score' in res or 'hits' in res:
                    self.con.execute("INSERT INTO results(job_id,seq_id,position,score,detail) VALUES(?,?,?,?,?)",
                                     (c.lastrowid, sid, res.get('position'), res.get('score', res.get('hits')), r['kind']))
            return self.run_id
    def index_experiment(self, n, reps=15):
        """Time a filter and a join on n synthetic rows with and without indexes."""
        rnd = random.Random(5)
        with self.lock, self.con:
            c = self.con
            for t in ('bench_rows', 'bench_samples'): c.execute(f"DROP TABLE IF EXISTS {t}")
            c.execute("CREATE TABLE bench_samples(sample_id INTEGER PRIMARY KEY, name TEXT)")
            c.execute("CREATE TABLE bench_rows(seq_id INTEGER PRIMARY KEY, sample_id INT REFERENCES bench_samples, length INT, gc REAL)")
            c.executemany("INSERT INTO bench_samples VALUES(?,?)", [(i, f"s{i}") for i in range(1000)])
            c.executemany("INSERT INTO bench_rows VALUES(?,?,?,?)", ((i, rnd.randrange(1000), rnd.randrange(200, 5000), rnd.random()) for i in range(n)))
        qs = {'Filter by sample': "SELECT COUNT(*), AVG(gc) FROM bench_rows WHERE sample_id=7",
              'Join with filter': "SELECT s.name, COUNT(*) FROM bench_rows q JOIN bench_samples s USING(sample_id) WHERE q.length>4900 GROUP BY s.name"}
        rows = []
        for ix in (False, True):
            with self.lock, self.con:
                self.con.execute("DROP INDEX IF EXISTS b1"); self.con.execute("DROP INDEX IF EXISTS b2")
                if ix: self.con.execute("CREATE INDEX b1 ON bench_rows(sample_id)"); self.con.execute("CREATE INDEX b2 ON bench_rows(length)")
            for name, q in qs.items(): rows.append(dict(query=name, indexed='With indexes' if ix else 'No indexes', ms=self.timed(q, reps)))
        return pd.DataFrame(rows)
