#!/usr/bin/env python
import subprocess
import io
from pathlib import Path
import pandas
import pytest
import os
from pytest import approx
#
R = Path(__file__).resolve().parents[1]
exe = './testdrv'
if os.name == 'nt':
    exe = exe[2:]


def test_fortran():

    def _csv2dat(txt):
        if isinstance(txt, str):
            txt1 = io.StringIO(txt)
            txt2 = io.StringIO(txt)
        else:
            txt1 = txt2 = txt

        comp = pandas.read_csv(txt1, sep='\s+', nrows=102,    index_col=0)
        light = pandas.read_csv(txt2, sep='\s+', skiprows=103, index_col=0)

        return pandas.concat((comp, light), axis=1, join='inner')

    ret = subprocess.check_output([exe, '-v'], cwd=R, universal_newlines=True)
    data = _csv2dat(ret)

    ref = _csv2dat(R/'reference/aur981.basic.out')

    assert data.values == approx(ref.values, rel=0.05)

    return data


if __name__ == '__main__':
    pytest.main(['-xrsv', __file__])
