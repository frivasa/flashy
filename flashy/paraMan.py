from flashy.IOutils import os
import functools
import importlib


def setupParaData(rc):
    balv = rc.load_balanced_view()
    dirv = rc.direct_view()
    with dirv.sync_imports():
        import flashy.paraMan
    return balv


def setupParaPost(rc):
    balv = rc.load_balanced_view()
    dirv = rc.direct_view()
    with dirv.sync_imports():
        import flashy.post as fpo
    return balv


def throwHammer(balv, N, func, **kwargs):
    kwargf = functools.partial(func, **kwargs)
    res = balv.map_async(kwargf, range(N))
    return res
