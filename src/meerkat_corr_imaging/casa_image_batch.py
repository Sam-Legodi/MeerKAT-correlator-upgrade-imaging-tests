"""Finite CASA imaging entrypoint (Python 2/3)."""
from __future__ import print_function
import os
import sys
import traceback

status = 1
try:
    path = os.environ['MCI_IMAGE_BATCH_SCRIPT']
    globals()['__file__'] = path
    # Preserve the arguments after CASA's -c bootstrap filename.
    index = next((i for i, arg in enumerate(sys.argv)
                  if os.path.basename(arg) == 'casa_image_batch.py'), None)
    if index is not None:
        sys.argv = [path] + sys.argv[index + 1:]
    with open(path, 'rb') as handle:
        exec(compile(handle.read(), path, 'exec'), globals())
    status = 0
except BaseException:
    traceback.print_exc()
finally:
    sys.stdout.flush()
    sys.stderr.flush()
    os._exit(status)
