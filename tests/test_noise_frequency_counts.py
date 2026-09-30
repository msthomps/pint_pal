"""Execute production attachment branches with real PINT noise parameters.

Isolate the DM/chromatic branches from enterprise/CHOLMOD imports; this checks
the actual property assignment, rather than only checking the source spelling.
Full Discovery-to-PINT covariance parity is outside this regression test.
"""
import ast
from pathlib import Path
from types import SimpleNamespace

import pint.models as pm
import pytest


SOURCE = Path(__file__).resolve().parents[1] / "src/pint_pal/noise_utils.py"
FUNCTION = next(
    node for node in ast.parse(SOURCE.read_text()).body
    if isinstance(node, ast.FunctionDef) and node.name == "add_noise_to_model"
)


@pytest.mark.parametrize("kind, variable, block, component, parameter", [
    ("dm", "dm_pars", "dm_noise", "PLDMNoise", "TNDMC"),
    ("chrom", "chrom_pars", "chromatic_noise", "PLChromNoise", "TNCHROMC"),
])
@pytest.mark.parametrize("requested", [None, 17, 100])
def test_attached_frequency_count(kind, variable, block, component, parameter, requested):
    branch = next(
        node for node in FUNCTION.body
        if isinstance(node, ast.If)
        and ast.unparse(node.test) == f"len({variable}) > 0"
    )
    attached = {}

    def add_component(value, **kwargs):
        attached[type(value).__name__] = value

    noise = {f"J1234_{kind}_gp_log10_A": -13., f"J1234_{kind}_gp_gamma": 3.}
    namespace = {
        "pm": pm,
        "psr_name": "J1234",
        "noise_dict": noise,
        variable: list(noise),
        "model_kwargs": {block: {} if requested is None else {"Nfreqs": requested}},
        "model": SimpleNamespace(components={"ChromaticCM": object()}, add_component=add_component),
        "log": SimpleNamespace(info=lambda *args: None),
        "_apply_pl_component_logfreq_settings": lambda *args, **kwargs: None,
    }
    exec(compile(ast.Module(body=[branch], type_ignores=[]), str(SOURCE), "exec"), namespace)
    count = getattr(attached[component], parameter)
    assert count.value == (100 if requested is None else requested)
    assert count.frozen
    assert not hasattr(count, "quantitity")
