#!/usr/bin/env bash
# 快速回归：单元测试 + 设计文档一致性校验（纯标准库，秒级，无容器/网络）
cd "$(dirname "$0")"
python3 -m unittest discover -s tests -v "$@" || exit 1
python3 check_design.py
