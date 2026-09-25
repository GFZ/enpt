"""Cross-platform task helpers for pixi."""

import glob
import os
import shutil
import subprocess
import sys
import webbrowser


def clean_build():
    for d in ["build", "dist", ".eggs"]:
        shutil.rmtree(d, ignore_errors=True)
    for p in glob.glob("**/*.egg-info", recursive=True):
        shutil.rmtree(p, ignore_errors=True)
    for p in glob.glob("**/*.egg", recursive=True):
        os.remove(p)


def clean_pyc():
    for pat in ["**/*.pyc", "**/*.pyo", "**/*~"]:
        for p in glob.glob(pat, recursive=True):
            os.remove(p)
    for p in glob.glob("**/__pycache__", recursive=True):
        shutil.rmtree(p, ignore_errors=True)


def clean_test():
    for f in [".coverage", "report.html", "report.xml", "coverage.xml"]:
        if os.path.exists(f):
            os.remove(f)
    for d in ["htmlcov", ".pytest_cache"]:
        shutil.rmtree(d, ignore_errors=True)
    for p in glob.glob(".coverage.*"):
        os.remove(p)


def docs(pkg):
    for f in [os.path.join("docs", pkg + ".rst"), os.path.join("docs", "modules.rst")]:
        if os.path.exists(f):
            os.remove(f)
    subprocess.run(
        ["sphinx-apidoc", pkg, "-o", "docs/", "--private",
         "--doc-project", "Python API reference"],
        check=True,
    )
    shutil.rmtree(os.path.join("docs", "_build"), ignore_errors=True)
    subprocess.run(
        ["sphinx-build", "-b", "html", "docs", os.path.join("docs", "_build", "html")],
        check=True,
    )
    webbrowser.open(os.path.join("docs", "_build", "html", "index.html"))


def dist():
    subprocess.run([sys.executable, "-m", "build"], check=True)
    for f in os.listdir("dist"):
        print(f)


if __name__ == "__main__":
    task = sys.argv[1] if len(sys.argv) > 1 else None
    extra = sys.argv[2] if len(sys.argv) > 2 else None
    tasks = {
        "clean-build": lambda: clean_build(),
        "clean-pyc": lambda: clean_pyc(),
        "clean-test": lambda: clean_test(),
        "docs": lambda: docs(extra) if extra else sys.exit("docs requires package name"),
        "dist": lambda: dist(),
    }
    if task in tasks:
        tasks[task]()
    else:
        print("Unknown task: {}. Available: {}".format(task, ", ".join(tasks)))
        sys.exit(1)
