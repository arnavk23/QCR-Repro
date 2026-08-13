import sys as _sys

from .config import DEFAULT_ANGLES, GateInstance, GateSetName
from .circuits import build_pool, random_circuit, count_gates
from .database import ComputeGraph, ReductionDatabase, load_or_build_database
from .token_pool import TokenPool
from .unitary import equivalent_up_to_global_phase
from .reducer import reduce_with_database, reduce_random_sampling, ReductionStats

# Module aliases so pickle caches written before the qcr_repro module renames
# (tokenizer -> token_pool, qasm_io -> qasm, unitary_utils -> unitary,
# compute_graph -> database, exact_graph -> exact_database) still unpickle.
_ALIASES = {
    "qcr_repro.compute_graph": "qcr_repro.database",
    "qcr_repro.tokenizer": "qcr_repro.token_pool",
    "qcr_repro.qasm_io": "qcr_repro.qasm",
    "qcr_repro.unitary_utils": "qcr_repro.unitary",
    "qcr_repro.exact_graph": "qcr_repro.exact_database",
}
import importlib as _importlib

for _old, _new in _ALIASES.items():
    try:
        _mod = _sys.modules[_new]
    except KeyError:
        _mod = _importlib.import_module(_new)
    _sys.modules.setdefault(_old, _mod)

__all__ = [
    "DEFAULT_ANGLES",
    "GateInstance",
    "GateSetName",
    "TokenPool",
    "build_pool",
    "random_circuit",
    "count_gates",
    "ComputeGraph",
    "ReductionDatabase",
    "load_or_build_database",
    "equivalent_up_to_global_phase",
    "reduce_with_database",
    "reduce_random_sampling",
    "ReductionStats",
]
