"""Runner policy tests with real PINT models and a recording fitter.

File-based loading keeps these prenoise tests independent of BBX's unrelated
enterprise/CHOLMOD imports. The repository's lite_utils and dmx_utils are used.
"""
import copy
import importlib.util
from io import StringIO
from pathlib import Path
import sys

from astropy.table import Table
from pint.models import get_model
from pint.models.noise_model import PLDMNoise, PLRedNoise
import pytest


ROOT = Path(__file__).resolve().parents[1] / "src" / "pint_pal"


def load_file(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


runner = load_file("classic_dmx_under_test", ROOT / "BBX" / "classic_dmx.py")
load_file("pint_pal.logger", ROOT / "logger.py")
runner.lu = load_file("classic_dmx_lite_utils", ROOT / "lite_utils.py")
runner.du = load_file("classic_dmx_dmx_utils", ROOT / "dmx_utils.py")


class Toas:
    def __init__(self):
        self.table = Table({"index": [5, 6], "mjd_float": [58001., 58011.]})

    def __getitem__(self, key):
        return self.table[key]


class Fitter:
    def __init__(self, toas, model):
        self.toas = toas
        self.model = copy.deepcopy(model)
        self.iterations = None

    def fit_toas(self, maxiter):
        self.iterations = maxiter


class Config:
    def __init__(self, free=None):
        self.free = free or ["F0", "DM", "DM1", "DM2", "DM3", "NE_SW"]

    def construct_fitter(self, to, mo):
        return Fitter(to, mo)

    def get_free_params(self, fitter):
        return self.free + [p for p in fitter.model.params if p.startswith("DMX_")]

    def get_niter(self):
        return 3


@pytest.fixture
def inputs():
    model = get_model(StringIO("""
PSR J1234+5678
F0 1
PEPOCH 58000
DM 10
DMEPOCH 58000
DM1 0.01
DM2 0.001
DM3 0.0001
NE_SW 7
SWM 1
SWP 2
DMX_0001 0.001
DMXR1_0001 58000
DMXR2_0001 58005
DMX_0002 0.002
DMXR1_0002 58010
DMXR2_0002 58015
"""))
    return dict(model=model, toas=Toas(), timing_config=Config())


def test_policy_and_input_isolation(inputs):
    original = inputs["model"].as_parfile(include_info=False)
    run = runner.run_classic_dmx(
        **inputs, dm_order=2, solar_wind="replace",
        sw_kwargs={"NE_SW": 6.67, "SWM": 1}, sw_free=("SWP",),
        frozen_params=("DM", "DM1", "DM2"), fit=False,
    )
    m = run.prenoise.model
    assert set(m.free_params) == {"F0", "SWP", "DMX_0001", "DMX_0002"}
    assert m.NE_SW.value == 6.67
    assert m.DM1.value == pytest.approx(.01)
    assert "DM3" not in m.params
    assert inputs["model"].as_parfile(include_info=False) == original
    assert inputs["timing_config"].free[-1] == "NE_SW"
    run.prenoise.toas.table["index"][0] = 999
    assert inputs["toas"].table["index"][0] == 5
    assert run.metadata["removed_yaml_entries"] == ["DM3"]


def test_all_noise_removed_and_fit_called(inputs):
    for component in (PLDMNoise(), PLRedNoise()):
        inputs["model"].add_component(component, validate=False)
    run = runner.run_classic_dmx(**inputs, dm_order=0, maxiter=5)
    assert set(run.metadata["removed_noise_components"]) == {"PLDMNoise", "PLRedNoise"}
    assert "PLDMNoise" in inputs["model"].components
    assert run.prenoise.iterations == 5
    assert run.fitted
    assert run.prenoise.model.DM1.value == 0
    assert run.prenoise.model.DM1.frozen
    assert "DM2" not in run.prenoise.model.params


def test_derivatives_added(inputs):
    run = runner.run_classic_dmx(
        **inputs, dm_order=4, unfrozen_params=("DMX_*", "DM[0-9]*"), fit=False,
    )
    assert run.prenoise.model.DM4.value == 0
    assert not run.prenoise.model.DM4.frozen


def test_sw_remove(inputs):
    run = runner.run_classic_dmx(**inputs, solar_wind="remove", dm_order=0, fit=False)
    assert "SolarWindDispersion" not in run.prenoise.model.components
    assert "NE_SW" in run.metadata["removed_yaml_entries"]


def test_degeneracy_and_reference_bin(inputs):
    with pytest.raises(ValueError, match="common-offset degeneracy"):
        runner.run_classic_dmx(**inputs, frozen_params=(), fit=False)
    run = runner.run_classic_dmx(
        **inputs, frozen_params=("DMX_0001",), unfrozen_params=("DM",), fit=False,
    )
    assert not run.prenoise.model.DM.frozen
    assert run.prenoise.model.DMX_0001.frozen


def test_empty_bin_not_a_reference_constraint(inputs):
    inputs["toas"].table = inputs["toas"].table[1:]
    with pytest.raises(ValueError, match="common-offset degeneracy"):
        runner.run_classic_dmx(
            **inputs, frozen_params=("DMX_0001",), unfrozen_params=("DM",), fit=False,
        )
    run = runner.run_classic_dmx(**inputs, fit=False)
    assert run.prenoise.model.DMX_0001.frozen
    assert "DMX_0001" in run.metadata["empty_mask_frozen"]


@pytest.mark.parametrize("kwargs, message", [
    ({"frozen_params": ("DM", "DMX_0001")}, "both frozen and unfrozen"),
    ({"unfrozen_params": ("TYPO",)}, "matches no model"),
    ({"solar_wind": "replace", "sw_kwargs": {"SWM": 0}, "sw_free": ("SWP",)}, "SWM=0"),
    ({"sw_free": ("SWM",)}, "model settings"),
    ({"sw_free": ("SWP",), "unfrozen_params": ("NE_SW",)}, "both frozen and unfrozen"),
])
def test_invalid_parameter_choices(inputs, kwargs, message):
    with pytest.raises(ValueError, match=message):
        runner.run_classic_dmx(**inputs, fit=False, **kwargs)


def test_rebuild_uses_returned_toas(inputs, monkeypatch):
    def setup(model, toas, **kwargs):
        assert kwargs == {"bin_width": 6.5}
        toas.table = toas.table[1:]
        return toas
    monkeypatch.setattr(runner.du, "setup_dmx", setup)
    run = runner.run_classic_dmx(
        **inputs, binning="rebuild", dmx_kwargs={"bin_width": 6.5}, fit=False,
    )
    assert run.metadata["used_toa_indices"] == [6]
    assert len(inputs["toas"].table) == 2


def test_fit_failure_propagates(inputs, monkeypatch):
    def fail(self, maxiter):
        raise RuntimeError("fit failed")
    monkeypatch.setattr(Fitter, "fit_toas", fail)
    with pytest.raises(RuntimeError, match="fit failed"):
        runner.run_classic_dmx(**inputs)
