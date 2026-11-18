#!/usr/bin/env python3
"""Tests for rotifer.core.session._files.

Run under pytest, or standalone: python test/core/session/test_files.py
"""

import sys

import pytest

from rotifer.core.session._files import (
    SessionLocked, atomic_write_bytes, atomic_write_text, atomic_writer,
    session_lock)


def test_write_bytes_and_text(tmp_path):
    atomic_write_bytes(tmp_path / "a", b"xyz")
    atomic_write_text(tmp_path / "b", "olá")
    assert (tmp_path / "a").read_bytes() == b"xyz"
    assert (tmp_path / "b").read_text(encoding="utf-8") == "olá"
    assert sorted(p.name for p in tmp_path.iterdir()) == ["a", "b"]


def test_error_keeps_old_content_and_leaves_no_temp(tmp_path):
    target = tmp_path / "f"
    target.write_bytes(b"old")
    with pytest.raises(RuntimeError):
        with atomic_writer(target) as handle:
            handle.write(b"new")
            raise RuntimeError("boom")
    assert target.read_bytes() == b"old"
    assert [p.name for p in tmp_path.iterdir()] == ["f"]


def test_error_on_new_file_creates_nothing(tmp_path):
    with pytest.raises(RuntimeError):
        with atomic_writer(tmp_path / "sub" / "f") as handle:
            handle.write(b"x")
            raise RuntimeError
    assert list((tmp_path / "sub").iterdir()) == []


def test_text_mode(tmp_path):
    with atomic_writer(tmp_path / "t", "w") as handle:
        handle.write("hi")
    assert (tmp_path / "t").read_text() == "hi"


def test_lock_contention(tmp_path):
    with session_lock(tmp_path, timeout=1):
        assert (tmp_path / ".lock").exists()
        with pytest.raises(SessionLocked):
            with session_lock(tmp_path, timeout=0.05):
                pass
    with session_lock(tmp_path, timeout=0.05):
        pass


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
