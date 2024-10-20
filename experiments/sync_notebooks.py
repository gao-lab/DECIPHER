import os
import subprocess
from argparse import ArgumentParser
from pathlib import Path
from loguru import logger


def get_args():
    parser = ArgumentParser()
    parser.add_argument("--dir", default=".", help="The directory to sync")
    parser.add_argument("-d", "--dry-run", action="store_true", default=False)
    parser.add_argument("-f", "--force", action="store_true", default=False)
    parser.add_argument("-i", "--inverse", action="store_true", default=False, help="From .py to .ipynb")
    return parser.parse_args()

    # iter the sub dir and grep all the ipynb files

def notebook2py(args):
    for notebook in Path().rglob("*.ipynb"):
        # check if .py file exists
        py_file = notebook.with_suffix(".py")
        if py_file.exists():
            logger.debug(f"{py_file} already exists")
            # compare the modification time of ipynb and py file
            if notebook.stat().st_mtime > py_file.stat().st_mtime or args.force:
                logger.warning(f"{notebook} updated, re-generate .py file")
                if not args.dry_run:
                    subprocess.run(f"rm {py_file}", shell=True)
            else:
                continue
        cmd = f"jupytext --set-formats ipynb,py:percent {notebook}"
        # print(cmd)
        if not args.dry_run and 'plot' not in cmd:
            subprocess.run(cmd, shell=True)


def py2notebook(args):
    for py in Path().rglob("*.py"):
        # check if .ipynb file exists
        notebook = py.with_suffix(".ipynb")
        if notebook.exists():
            logger.debug(f"{notebook} already exists")
            # compare the modification time of ipynb and py file
            if py.stat().st_mtime > notebook.stat().st_mtime or args.force:
                logger.warning(f"{py} updated, update .ipynb file")
                cmd = f"jupytext --to notebook {py}"
            else:
                continue
        else:
            cmd = f"jupytext --to notebook {py}"
        # print(cmd)
        if not args.dry_run:
            subprocess.run(cmd, shell=True)


if __name__ == "__main__":
    args = get_args()
    os.chdir(args.dir)
    if args.inverse:
        py2notebook(args)
    else:
        notebook2py(args)
