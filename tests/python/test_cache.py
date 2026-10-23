"""User cache precedence is independent of the process working directory."""

import pytest

from beach.tooling._cache import default_periodic_cache_dir


@pytest.mark.parametrize(
    ("beach", "xdg", "expected"),
    [
        ("/project/cache", "/user/cache", "/project/cache/periodic2"),
        ("", "/user/cache", "/user/cache/beach/periodic2"),
        (None, "", "/user/home/.cache/beach/periodic2"),
    ],
)
def test_default_periodic_cache_dir(monkeypatch, tmp_path, beach, xdg, expected):
    monkeypatch.setenv("HOME", "/user/home")
    for key, value in (("BEACH_CACHE_DIR", beach), ("XDG_CACHE_HOME", xdg)):
        if value is None:
            monkeypatch.delenv(key, raising=False)
        else:
            monkeypatch.setenv(key, value)
    monkeypatch.chdir(tmp_path)
    assert default_periodic_cache_dir() == expected
