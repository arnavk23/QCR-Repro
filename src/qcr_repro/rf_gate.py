"""V3-style RF-gated database lookup (Rosenhahn et al., variant V3).

Skips DB lookups on blocks predicted irreducible: an exact memo cache handles repeated windows; a lazily-trained RandomForest (sklearn, optional [ml] extra) scores novel ones. Skipping only reallocates the time budget, never changes the circuit. Without sklearn it degrades to memo-only mode."""

from __future__ import annotations

import math
from collections import OrderedDict

from .config import GateInstance

_CLIFFORD = (math.pi / 2, -math.pi / 2)
_ANGLE_TOL = 1e-6

# Gate types counted as features (anything else buckets into "other").
_FEATURE_TYPES = ("RX", "RY", "RZ", "RXX", "CZ")


def _signature(block: list[GateInstance]) -> tuple:
    """Exact memo key for a block (angle-rounded, order-preserving)."""
    return tuple(
        (g.name, tuple(sorted(g.qubits)), None if g.theta is None else round(float(g.theta), 8))
        for g in block
    )


def _features(block: list[GateInstance]) -> tuple[float, ...]:
    """Cheap feature vector describing a block for the classifier."""
    n = len(block)
    twq = 0
    wires: set[int] = set()
    counts = {t: 0 for t in _FEATURE_TYPES}
    other = 0
    non_clifford = 0
    for g in block:
        if len(g.qubits) == 2:
            twq += 1
        wires.update(g.qubits)
        if g.name in counts:
            counts[g.name] += 1
        else:
            other += 1
        if g.theta is not None and not any(abs(g.theta - a) <= _ANGLE_TOL for a in _CLIFFORD):
            non_clifford += 1
    return (
        float(n),
        float(twq),
        float(len(wires)),
        *[float(counts[t]) for t in _FEATURE_TYPES],
        float(other),
        float(non_clifford),
    )


class RfGate:
    """Online RF-gated lookup decision maker.

Call decide() before a lookup and observe(outcome) after; trains on first-time outcomes only."""

    def __init__(
        self,
        seed: int = 0,
        min_train_samples: int = 512,
        retrain_every: int = 4096,
        threshold: float = 0.2,
        n_estimators: int = 32,
        max_memo: int = 262144,
    ) -> None:
        """threshold = minimum predicted reducibility probability to attempt a novel block.

Lower is more conservative. 0.2 is the least-bad setting measured on the exhaustive sweep (aggressive gating causes premature convergence)."""
        self._seed = seed
        self._min_train = min_train_samples
        self._retrain_every = retrain_every
        self._threshold = threshold
        self._n_estimators = n_estimators
        self._max_memo = max_memo
        self._memo: dict[tuple, bool] = {}
        self._features: list[tuple[float, ...]] = []
        self._labels: list[int] = []
        self._model = None
        self._samples_since_train = 0
        # statistics (exposed for the benchmark report)
        self.lookups_attempted = 0
        self.lookups_skipped = 0
        self.reductions_found = 0

    # ------------------------------------------------------------------ #
    # public interface
    # ------------------------------------------------------------------ #

    def decide(self, block: list[GateInstance]) -> bool:
        """True if the DB lookup for ``block`` should be attempted."""
        sig = _signature(block)
        known = self._memo.get(sig)
        if known is not None:
            return known
        if self._model is None or len(self._features) < self._min_train:
            return True
        return self._predict(block) >= self._threshold

    def observe(self, block: list[GateInstance], reduced: bool) -> None:
        """Record the true outcome of a lookup on ``block``."""
        sig = _signature(block)
        if sig in self._memo:
            return
        if len(self._memo) >= self._max_memo:
            # bounded memo: drop the oldest half (dicts preserve insertion order)
            for _ in range(len(self._memo) // 2):
                self._memo.pop(next(iter(self._memo)))
        self._memo[sig] = reduced
        self._features.append(_features(block))
        self._labels.append(1 if reduced else 0)
        self._samples_since_train += 1
        if reduced:
            self.reductions_found += 1
        if (self._model is None and len(self._features) >= self._min_train) or (
            self._model is not None and self._samples_since_train >= self._retrain_every
        ):
            self._train()

    # ------------------------------------------------------------------ #
    # internals
    # ------------------------------------------------------------------ #

    def _train(self) -> None:
        try:
            import numpy as np
            from sklearn.ensemble import RandomForestClassifier
        except ImportError:
            return
        if len(self._labels) < self._min_train:
            return
        x = np.asarray(self._features, dtype=float)
        y = np.asarray(self._labels, dtype=int)
        model = RandomForestClassifier(
            n_estimators=self._n_estimators,
            max_depth=8,
            min_samples_leaf=8,
            n_jobs=1,
            random_state=self._seed,
            class_weight="balanced",
        )
        model.fit(x, y)
        self._model = model
        self._samples_since_train = 0

    def _predict(self, block: list[GateInstance]) -> float:
        model = self._model
        if model is None:
            return 1.0
        try:
            import numpy as np

            x = np.asarray([list(_features(block))], dtype=float)
            p = model.predict_proba(x)[0]
            classes = list(model.classes_)
            if 1 in classes:
                return float(p[classes.index(1)])
            return 0.5
        except Exception:
            return 1.0  # on any prediction failure, fall back to attempting


class RfGatedDatabase:
    """Proxy around ReductionDatabase gating try_reduce / try_reduce_cost.

Forwarded attributes make it drop-in for reduce_circuit; escape resampling stays ungated."""

    def __init__(self, db, gate: RfGate):
        self._db = db
        self._gate = gate

    def __getattr__(self, name):
        return getattr(self._db, name)

    def _go(self, block: list[GateInstance]) -> bool:
        if self._gate.decide(block):
            self._gate.lookups_attempted += 1
            return True
        self._gate.lookups_skipped += 1
        return False

    def try_reduce(self, block: list[GateInstance]) -> list[GateInstance] | None:
        if not self._go(block):
            return None
        candidate = self._db.try_reduce(block)
        self._gate.observe(block, candidate is not None and len(candidate) < len(block))
        return candidate

    def try_reduce_cost(self, block: list[GateInstance]) -> list[GateInstance] | None:
        if not self._go(block):
            return None
        candidate = self._db.try_reduce_cost(block)
        self._gate.observe(block, candidate is not None)
        return candidate

    def try_reduce_escape(self, block, rng, slack: int = 3, prefer=None):
        return self._db.try_reduce_escape(block, rng, slack, prefer)
