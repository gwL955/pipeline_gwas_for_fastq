#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""线程安全日志：主日志 + 每样本日志，时间戳 + 耗时。仅标准库。"""

import os
import sys
import threading
from datetime import datetime


class Logger:
    """双通道日志：stdout + 文件（log_file=None 时仅控制台，配合 capture_stdio
    由 tee 统一落盘）；进程内全局锁保证线程安全。"""

    _lock = threading.Lock()

    def __init__(self, log_file=None, prefix="MAIN", quiet=False):
        self._fh = None
        if log_file is not None:
            os.makedirs(os.path.dirname(log_file), exist_ok=True)
            self._fh = open(log_file, "a", encoding="utf-8")
        self._fh_lock = threading.Lock()
        self._start = datetime.now()
        self._prefix = prefix
        self._quiet = quiet

    # ── 基础 ──
    def _ts(self):
        return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    def _elapsed(self):
        e = (datetime.now() - self._start).total_seconds()
        if e < 60:
            return f"{e:.1f}s"
        m, s = divmod(int(e), 60)
        h, m = divmod(m, 60)
        return f"{int(h)}h{int(m)}m{s}s"

    def _emit(self, level, msg):
        tag = f"[{self._prefix}]" if self._prefix != "MAIN" else ""
        line = f"[{self._ts()}] [{level:5s}] [elapsed={self._elapsed()}] {tag} {msg}"
        if not self._quiet:
            with Logger._lock:
                print(line, flush=True)
        if self._fh is not None:
            with self._fh_lock:
                self._fh.write(line + "\n")
                self._fh.flush()
        return line

    def info(self, msg):   self._emit("INFO", msg)
    def step(self, msg):   self._emit("STEP", msg)
    def cmd(self, msg):    self._emit("CMD", msg)
    def out(self, msg):    self._emit("OUT", msg)
    def warn(self, msg):   self._emit("WARN", msg)
    def error(self, msg):  self._emit("ERROR", msg)
    def result(self, msg): self._emit("RESULT", msg)
    def skip(self, msg):   self._emit("INFO", f"[SKIP] {msg}")

    def tail(self, n=20):
        """日志文件末尾 n 行（失败通知用）"""
        if self._fh is None:
            return ""
        try:
            with open(self._fh.name, encoding="utf-8") as f:
                return "".join(f.readlines()[-n:]).rstrip()
        except OSError:
            return ""

    def close(self):
        try:
            if self._fh is not None:
                self._fh.close()
        except OSError:
            pass


class TeeStream:
    """把 stdout/stderr 同步镜像（tee）到文件的流包装——任何裸 print、
    未捕获 traceback、子进程继承输出都会自动落盘，无需 shell 重定向。
    path=None 时先写入内存缓冲，tee_set_log_path() 确定目标后一并落盘
    （运行日志归属批次结果目录，而目标目录要等批次发现后才知道；
    dry-run/启动即退场景从不设路径 → 缓冲丢弃，零落盘）。"""

    def __init__(self, stream, path=None):
        self._stream = stream
        self._buf = []
        self._fh = None
        if path is not None:
            self._open(path)

    def _open(self, path):
        os.makedirs(os.path.dirname(path), exist_ok=True)
        self._fh = open(path, "a", encoding="utf-8", buffering=1)
        if self._buf:
            self._fh.write("".join(self._buf))
            self._buf = []

    def set_path(self, path):
        if self._fh is None:
            self._open(path)

    def write(self, data):
        try:
            self._stream.write(data)
        except (ValueError, OSError):
            pass   # 控制台已关闭（nohup 后台）时只落盘
        if self._fh is not None:
            self._fh.write(data)
        else:
            self._buf.append(data)
        return len(data)

    def flush(self):
        try:
            self._stream.flush()
        except (ValueError, OSError):
            pass
        if self._fh is not None:
            self._fh.flush()

    def isatty(self):
        return False

    def __getattr__(self, name):
        return getattr(self._stream, name)


def capture_stdio(log_path=None):
    """接管进程 stdout/stderr：控制台照常显示，同时镜像到运行日志。
    log_path=None → 缓冲模式，稍后 tee_set_log_path() 落盘"""
    sys.stdout = TeeStream(sys.stdout, log_path)
    sys.stderr = TeeStream(sys.stderr, log_path)


def tee_set_log_path(path):
    """确定运行日志落盘位置（幂等；stdout/stderr 两个流写同一文件）"""
    for s in (sys.stdout, sys.stderr):
        if isinstance(s, TeeStream):
            s.set_path(path)


class SampleLoggerFactory:
    """每样本独立日志：logs/sample_<样本>.log（按批次隔离目录）"""

    def __init__(self):
        self._loggers = {}
        self._lock = threading.Lock()

    def get(self, sample, log_dir):
        with self._lock:
            key = (log_dir, sample)
            if key not in self._loggers:
                path = os.path.join(log_dir, f"sample_{sample}.log")
                self._loggers[key] = Logger(path, prefix=sample)
            return self._loggers[key]

    def close_all(self):
        with self._lock:
            for lg in self._loggers.values():
                lg.close()
            self._loggers.clear()
