import functools
from .IOutils import log


def setupParaData(rc):
    balv = rc.load_balanced_view()
    return balv


def setupParaPost(rc):
    balv = rc.load_balanced_view()
    return balv


def throwHammer(balv, N, func, timeout=None, propagate_errors=False, **kwargs):
    """Execute func in parallel across N workers.

    Args:
        balv: Load-balanced view from ipp
        N: Number of parallel tasks
        func: Function to execute (takes wedgenum as first arg)
        timeout: Optional timeout in seconds
        propagate_errors: If True, raise first exception; if False, log and continue
        **kwargs: Arguments passed to func

    Returns:
        AsyncResult if timeout is None, otherwise the result of .get()
    """
    kwargf = functools.partial(func, **kwargs)
    res = balv.map_async(kwargf, range(N))

    if timeout is None:
        return res

    try:
        return res.get(timeout=timeout)
    except Exception as e:
        log.critical('Parallel task failed: {}: {}'.format(type(e).__name__, e),
                     exc_info=True)
        if propagate_errors:
            raise
        return None
