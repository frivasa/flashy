import functools


def setupParaData(rc):
    balv = rc.load_balanced_view()
    dirv = rc.direct_view()
    with dirv.sync_imports():
        pass
    return balv


def setupParaPost(rc):
    balv = rc.load_balanced_view()
    dirv = rc.direct_view()
    with dirv.sync_imports():
        pass
    return balv


def throwHammer(balv, N, func, **kwargs):
    kwargf = functools.partial(func, **kwargs)
    res = balv.map_async(kwargf, range(N))
    return res
