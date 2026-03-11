from __future__ import annotations

from collections import deque
from dataclasses import dataclass

import numpy as np

from .gates import embedded_gate_matrix
from .tokenizer import TokenPool
from .unitary_utils import unitary_key


@dataclass(frozen=True)
class GraphNode:
    key: tuple[tuple[float, float], ...]
    depth: int
    token_chain: tuple[int, ...]


class ComputeGraphBuilder:
    def __init__(self, token_pool: TokenPool, max_depth: int) -> None:
        self.token_pool = token_pool
        self.max_depth = max_depth
        self.num_qubits = token_pool.num_qubits
        self.dim = 2**self.num_qubits

    def build(self) -> dict[tuple[tuple[float, float], ...], GraphNode]:
        root_unitary = np.eye(self.dim, dtype=complex)
        root_key = unitary_key(root_unitary)
        best: dict[tuple[tuple[float, float], ...], GraphNode] = {
            root_key: GraphNode(key=root_key, depth=0, token_chain=()),
        }

        queue: deque[tuple[np.ndarray, tuple[int, ...], int]] = deque()
        queue.append((root_unitary, (), 0))

        cached_gate_mats = {
            token: embedded_gate_matrix(self.num_qubits, self.token_pool.gate_for_token(token))
            for token in self.token_pool.tokens()
        }

        while queue:
            current_u, current_chain, depth = queue.popleft()
            if depth >= self.max_depth:
                continue

            for token, gate_u in cached_gate_mats.items():
                next_u = gate_u @ current_u
                next_chain = current_chain + (token,)
                key = unitary_key(next_u)

                existing = best.get(key)
                if existing is not None and len(existing.token_chain) <= len(next_chain):
                    continue

                node = GraphNode(key=key, depth=depth + 1, token_chain=next_chain)
                best[key] = node
                queue.append((next_u, next_chain, depth + 1))

        return best
