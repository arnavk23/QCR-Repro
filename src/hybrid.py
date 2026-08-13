"""Exact/numeric hybrid lookup for NISQ pools.

Routes fully-Clifford windows (RX/RZ at +/-pi/2, CZ) to the exact symplectic engine and everything else to the numeric database -- deep exact reduction where possible, numeric fallback elsewhere."""

from __future__ import annotations

import math

from .config import GateInstance

_CLIFFORD_ANGLES = (math.pi / 2, -math.pi / 2)


def is_clifford_window(block: list[GateInstance], atol: float = 1e-6) -> bool:
    """True if every gate in block is in the Clifford sub-pool (RX/RZ at +/-pi/2, or CZ)."""
    for gate in block:
        if gate.name == "CZ":
            continue
        if gate.name in ("RX", "RZ"):
            if gate.theta is None:
                return False
            if not any(abs(gate.theta - a) <= atol for a in _CLIFFORD_ANGLES):
                return False
        else:
            return False
    return True


class HybridDatabase:
    """Routes lookups to an exact (Clifford) or numeric database per window.

Forwards attribute access to the numeric database so it drops into reduce_circuit; an optional RfGate wraps the numeric fallback."""

    def __init__(self, numeric_db, exact_db, numeric_gate=None, atol: float = 1e-6):
        self._numeric = numeric_db
        self._exact = exact_db
        self._gate = numeric_gate
        self._atol = atol
        # statistics (exposed for the benchmark report)
        self.exact_lookups = 0
        self.numeric_lookups = 0

    def __getattr__(self, name):
        return getattr(self._numeric, name)

    def _route(self, block: list[GateInstance]) -> bool:
        return is_clifford_window(block, self._atol)

    def _numeric_try(self, block, cost: bool):
        self.numeric_lookups += 1
        if self._gate is None:
            if cost:
                return self._numeric.try_reduce_cost(block)
            return self._numeric.try_reduce(block)
        if self._gate.decide(block):
            self._gate.lookups_attempted += 1
            if cost:
                candidate = self._numeric.try_reduce_cost(block)
            else:
                candidate = self._numeric.try_reduce(block)
            self._gate.observe(block, candidate is not None)
            return candidate
        self._gate.lookups_skipped += 1
        return None

    def try_reduce(self, block: list[GateInstance]) -> list[GateInstance] | None:
        if self._route(block):
            candidate = self._exact.try_reduce(block)
            if candidate is not None:
                self.exact_lookups += 1
                return candidate
            # Clifford window the exact graph cannot improve: try the numeric
            # DB (it may hold the same word at a different depth boundary).
        return self._numeric_try(block, cost=False)

    def try_reduce_cost(self, block: list[GateInstance]) -> list[GateInstance] | None:
        if self._route(block):
            candidate = self._exact.try_reduce_cost(block)
            if candidate is not None:
                self.exact_lookups += 1
                return candidate
        return self._numeric_try(block, cost=True)

    def try_reduce_escape(self, block, rng, slack: int = 3, prefer=None):
        if self._route(block):
            candidate = self._exact.try_reduce_escape(block, rng, slack, prefer)
            if candidate is not None:
                self.exact_lookups += 1
                return candidate
        self.numeric_lookups += 1
        return self._numeric.try_reduce_escape(block, rng, slack, prefer)
