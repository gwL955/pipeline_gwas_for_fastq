#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""multiqc：聚合 qc/ 全部 QC 产物（fastqc/fastp/flagstat/stats/markdup/
hsmetrics/mosdepth/bcftools_stats）。对应笔记《5-测序质量》代码块 7。"""

import glob
import os

import config


def run_multiqc(runner, qc_dir_host, out_dir_host, logger):
    """multiqc 报告文件名带标题前缀（空格→'-'），按 glob 检查实现幂等"""
    existing = glob.glob(os.path.join(out_dir_host, "*multiqc_report.html"))
    if existing and not runner.dry_run:
        logger.skip(f"MultiQC 报告已存在: {os.path.basename(existing[-1])}")
        return True
    rc = runner.run(
        runner.tool("multiqc",
                    f"multiqc {runner.cpath(qc_dir_host)} "
                    f"-o {runner.cpath(out_dir_host)} "
                    f"--title {config.MULTIQC_TITLE!r} --force"),
        logger=logger)
    return rc == 0 and bool(glob.glob(os.path.join(out_dir_host,
                                                   "*multiqc_report.html")))
