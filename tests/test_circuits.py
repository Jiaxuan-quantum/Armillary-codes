"""Small detector/observable mapping tests; no Monte Carlo campaign is run."""
import numpy as np
import pytest
import stim

from armillary.models import DetectorModel, _dem_arrays, load_model


def test_complete_fault_effects_and_parity():
    # Distinct logical effects with identical syndrome must not be merged.
    # A detectorless logical fault must remain available to the decoder.
    dem = stim.DetectorErrorModel("""
        error(0.1) D0 L0
        error(0.2) D0 L1
        error(0.03) L0
        error(0.04) D0 D0 L1 L1
        error(0.05) D0 ^ D1 L0
    """)
    h, probabilities, logical = _dem_arrays(dem, 2, 2)
    np.testing.assert_array_equal(h.toarray(), [[1, 1, 0, 0, 1], [0, 0, 0, 0, 1]])
    np.testing.assert_allclose(probabilities, [0.1, 0.2, 0.03, 0.04, 0.05])
    np.testing.assert_array_equal(logical[:, 0], [1, 2, 1, 0, 1])


def test_circuit_model_roundtrip(tmp_path):
    circuit = stim.Circuit("""
        R 0 1
        X_ERROR(0.1) 0
        X_ERROR(0.2) 1
        M 0 1
        DETECTOR rec[-2]
        OBSERVABLE_INCLUDE(0) rec[-1]
    """)
    model = DetectorModel.from_circuit(circuit)
    assert model.H.shape == (1, 2)
    assert any(np.diff(model.H.indptr) == 0)
    model.save(tmp_path)
    restored = load_model(tmp_path)
    np.testing.assert_array_equal(restored.H.toarray(), model.H.toarray())
    np.testing.assert_array_equal(restored.logical_map, model.logical_map)
    for e in ([0, 0], [1, 0], [0, 1], [1, 1]):
        np.testing.assert_array_equal(restored.predict_observables(e),
                                      model.predict_observables(e))
    assert restored.metadata["model_id"] == model.metadata["model_id"]
    restored.save(tmp_path)  # Reusing the same verified model is safe.


def test_model_checksum_rejects_modified_file(tmp_path):
    model = DetectorModel.from_circuit(stim.Circuit("R 0\nX_ERROR(0.1) 0\nM 0\nDETECTOR rec[-1]"))
    model.save(tmp_path)
    with (tmp_path / "probabilities.npy").open("ab") as stream:
        stream.write(b"modified")
    with pytest.raises(ValueError, match="checksum"):
        load_model(tmp_path)


def test_prediction_rejects_nonbinary_vector():
    model = DetectorModel.from_circuit(stim.Circuit("R 0\nX_ERROR(0.1) 0\nM 0\nOBSERVABLE_INCLUDE(0) rec[-1]"))
    with pytest.raises(ValueError, match="binary"):
        model.predict_observables([2])


def test_noiseless_model_has_no_fault_variables():
    model = DetectorModel.from_circuit(stim.Circuit("R 0\nM 0\nDETECTOR rec[-1]"))
    assert model.H.shape == (1, 0)
    assert model.probabilities.shape == (0,)
    assert model.predict_observables([]).shape == (0,)
