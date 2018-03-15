#!/usr/bin/env python
import subprocess
import io
from pathlib import Path
import pandas
import numpy as np
#
root = Path(__file__).parents[1]


def test_fortran():

    def _csv2dat(txt):
        if isinstance(txt,str):
            txt1 = io.StringIO(txt)
            txt2 = io.StringIO(txt)
        else:
            txt1=txt2=txt

        comp  = pandas.read_csv(txt1, sep='\s+', nrows=102,    index_col=0)
        light = pandas.read_csv(txt2, sep='\s+', skiprows=103, index_col=0)

        return pandas.concat((comp,light),axis=1, join='inner')


    ret = subprocess.check_output(['./testdrv','-v'], cwd=root, universal_newlines=True)
    data = _csv2dat(ret)

    ref = _csv2dat(root/'reference/aur981.basic.out')

    np.testing.assert_allclose(data,ref)

    return data

if __name__ == '__main__':
    data = test_fortran()