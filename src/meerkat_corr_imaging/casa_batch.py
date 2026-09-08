"""CASA 5/6 batch bootstrap. Executed inside CASA, not imported by the CLI."""
from __future__ import print_function
import os
import sys
import json
import traceback

result = {'run_id': os.environ['MCI_RUN_ID'], 'ms': os.environ['MCI_CAL_MSFILE'],
          'calibration_applied': False, 'quality': 'not_run', 'tasks': []}
result_path = os.environ['MCI_RESULT_PATH']
code = 1
try:
    helper = os.path.join(os.environ['MCI_CASA_SCRIPT_DIR'], 'calibration_checks.py')
    checks = {'__name__': 'calibration_checks'}
    with open(helper, 'rb') as handle:
        exec(compile(handle.read(), helper, 'exec'), checks)
    # Keep CASA's original task registry and bootstrap state intact.
    script_namespace = dict(globals())
    checks['install_checks'](script_namespace, result)
    script = os.environ['MCI_BATCH_SCRIPT']
    script_namespace['__file__'] = script
    with open(script, 'rb') as handle:
        exec(compile(handle.read(), script, 'exec'), script_namespace)
    if not result['calibration_applied']:
        raise RuntimeError('Script ended without validated calibration application')
    code = 0
except BaseException as exc:
    result['error'] = str(exc)
    traceback.print_exc()
finally:
    try:
        with open(result_path + '.tmp', 'w') as handle:
            json.dump(result, handle, indent=2)
        os.rename(result_path + '.tmp', result_path)
    except BaseException:
        traceback.print_exc()
        code = 1
    sys.stdout.flush()
    sys.stderr.flush()
    # CASA/IPython can catch SystemExit and return to a prompt. All task tools
    # have completed/closed before this batch-only process exit.
    os._exit(code)
