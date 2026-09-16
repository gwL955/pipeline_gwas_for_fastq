#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""singularity(apptainer) 命令封装：实时输出采集、超时、返回码、
产物存在性检查（幂等 SKIP）、容器路径映射。仅标准库。"""

import os
import shlex
import shutil
import threading
import subprocess

import config


def nonempty(path):
    return bool(path) and os.path.isfile(str(path)) and os.path.getsize(str(path)) > 0


class Runner:
    def __init__(self, dry_run=False, extra_binds=None):
        self.work = config.WORK_DIR
        self.rt = self._resolve_rt(config.CONTAINER_RT)
        self.dry_run = dry_run
        self.extra_binds = list(extra_binds or [])   # [(host, container)]
        self._env = self._hardened_env()

    def _resolve_rt(self, name):
        """容器运行时解析为绝对路径：不依赖调用方 PATH 的完整性
        （实测：从 PATH 受限环境启动时 `/bin/sh: singularity: not found`
        exit=127，而本机实际装在 /usr/local/bin/singularity）"""
        exe = shutil.which(name)
        if exe:
            return exe
        for d in ("/usr/local/bin", "/usr/local/sbin", "/usr/bin", "/usr/sbin",
                  os.path.expanduser("~/.local/bin"), os.path.expanduser("~/bin")):
            cand = os.path.join(d, name)
            if os.path.isfile(cand) and os.access(cand, os.X_OK):
                return cand
        return name   # 全部找不到则保留原名，让 127 报错可读

    def _hardened_env(self):
        """继承完整环境，但 PATH 缺标准系统目录时补齐
        （精简环境可能连 cat/grep 都找不到；singularity 也会把该 PATH 带进容器）"""
        env = {**os.environ, "LC_ALL": "C"}
        have = env.get("PATH", "")
        missing = [d for d in ("/usr/local/sbin", "/usr/local/bin",
                               "/usr/sbin", "/usr/bin", "/sbin", "/bin")
                   if d not in have.split(":")]
        if missing:
            env["PATH"] = ":".join(missing + ([have] if have else []))
        return env

    # ── 路径映射：宿主路径 → 容器路径 ──
    def cpath(self, path):
        p = str(path)
        if p.startswith(self.work):
            return "/data" + p[len(self.work):]
        for h, c in self.extra_binds:
            if p.startswith(h):
                return c + p[len(h):]
        return p

    # ── 拼装容器命令 ──
    def tool(self, sif_key, args, binds=None):
        bind_args = [f"--bind {shlex.quote(self.work)}:/data"]
        for h, c in list(self.extra_binds) + list(binds or []):
            bind_args.append(f"--bind {shlex.quote(h)}:{c}")
        return f"{self.rt} exec {' '.join(bind_args)} {shlex.quote(config.SIF[sif_key])} {args}"

    # ── 执行 ──
    def run(self, cmd, logger=None, outputs=(), timeout=None, capture=False):
        """执行 shell 命令。outputs 为本命令的宿主侧产物列表：
        全部存在且非空 → [SKIP]；否则执行并校验产物落盘。"""
        lg = logger
        outputs = [str(o) for o in outputs if o]
        if outputs and all(nonempty(o) for o in outputs):
            if lg:
                lg.skip(f"产物已存在，跳过命令: {os.path.basename(outputs[0])}"
                        + (" 等" if len(outputs) > 1 else ""))
            return (0, "", "") if capture else 0
        if self.dry_run:
            if lg:
                lg.cmd(f"[DRY-RUN] {cmd}")
            return (0, "", "") if capture else 0
        if lg:
            lg.cmd(cmd)
        try:
            proc = subprocess.Popen(
                cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                text=True, bufsize=1, env=self._env,
            )
        except Exception as e:   # noqa: BLE001
            if lg:
                lg.error(f"命令启动失败: {e}")
            return (-2, "", "") if capture else -2
        timed_out = {"fired": False}

        def _timeout_kill():
            timed_out["fired"] = True
            try:
                proc.kill()
            except OSError:
                pass

        # 超时用后台定时器强制 kill（无输出的长命令不会触发 wait 的 TimeoutExpired，
        # 因为 stdout 排空本身会阻塞到进程退出——曾致 sleep 类命令超时不生效）
        timer = threading.Timer(timeout, _timeout_kill) \
            if timeout is not None else None
        try:
            if timer is not None:
                timer.start()
            if capture:
                stdout, _ = proc.communicate(
                    timeout=None if timeout is None else timeout)
                rc = proc.returncode
            else:
                out_chunks = []
                for line in proc.stdout:
                    line = line.rstrip()
                    if lg:
                        lg.out(line)
                    out_chunks.append(line)
                proc.wait()
                rc = proc.returncode
                stdout = "\n".join(out_chunks)
            if timed_out["fired"]:
                rc = -1
            if rc != 0 and lg:
                lg.error(f"命令失败 (exit={rc}): {cmd[:300]}")
            if rc == 0 and outputs:
                missing = [o for o in outputs if not nonempty(o)]
                if missing and lg:
                    lg.warn(f"命令成功但产物缺失/为空: {missing}")
            if capture:
                return rc, stdout or "", ""
            return rc
        except subprocess.TimeoutExpired:
            proc.kill()
            if lg:
                lg.error(f"命令超时(>{timeout}s) 已杀死: {cmd[:300]}")
            return (-1, "", "") if capture else -1
        except Exception as e:   # noqa: BLE001
            if lg:
                lg.error(f"命令异常: {e}")
            return (-2, "", "") if capture else -2
        finally:
            if timer is not None:
                timer.cancel()

    def out(self, cmd, logger=None):
        """执行并仅返回 stdout（一次性统计查询用）"""
        rc, o, _ = self.run(cmd, logger=logger, capture=True)
        return o if rc == 0 else ""
