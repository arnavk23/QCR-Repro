from __future__ import annotations

from .config import GateInstance


def collect_wire_blocks(gates: list[GateInstance], max_wires: int = 3) -> list[list[int]]:
    """Partition gate indices into maximal <=max_wires-wire blocks, respecting
    the true per-wire dependency DAG rather than list adjacency.

    Standard block-collection over a circuit DAG (cf. Qiskit's 2Q-block
    collector), generalized to k<=max_wires wire footprints: a block grows by
    absorbing any newly-scanned gate whose wires overlap it, as long as the
    combined wire footprint stays within max_wires; a gate that would push an
    open block over the limit forces that block closed instead.

    Safety invariant used by :func:`compact_by_blocks`: every block's first
    member is always merged toward the *most recently opened* touched block
    (highest id).  Since ids are assigned in scan order, that block's anchor
    (its own first-member index) is >= every other touched block's anchor, so
    sorting the returned blocks by first-member index yields a valid
    topological order of the circuit's dependency DAG -- see module docstring
    reasoning replicated in check_dag_compact.py.  Returns a list of blocks,
    each a list of original indices in ascending order; every index appears
    in exactly one block.
    """
    block_members: dict[int, list[int]] = {}
    block_wires: dict[int, set[int]] = {}
    open_for_wire: dict[int, int] = {}
    next_id = 0

    def new_block(i: int, wires: set[int]) -> None:
        nonlocal next_id
        bid = next_id
        next_id += 1
        block_members[bid] = [i]
        block_wires[bid] = set(wires)
        for w in wires:
            open_for_wire[w] = bid

    def close(bid: int) -> None:
        for w in block_wires[bid]:
            if open_for_wire.get(w) == bid:
                del open_for_wire[w]

    for i, gate in enumerate(gates):
        wires = set(gate.qubits)
        touched = sorted({open_for_wire[w] for w in wires if w in open_for_wire})

        if not touched:
            new_block(i, wires)
            continue

        # Merging must always prefer the touched block with the largest id
        # (most recently opened): its anchor is >= every other touched
        # block's anchor, which is what keeps the final sort order valid.
        anchor_block = touched[-1]
        union = set(block_wires[anchor_block]) | wires
        if len(union) > max_wires:
            for b in touched:
                close(b)
            new_block(i, wires)
            continue

        chosen = [anchor_block]
        for b in sorted((t for t in touched if t != anchor_block), reverse=True):
            candidate = union | block_wires[b]
            if len(candidate) <= max_wires:
                union = candidate
                chosen.append(b)
            else:
                close(b)

        target = anchor_block
        members = list(block_members[target])
        for b in chosen:
            if b == target:
                continue
            members.extend(block_members[b])
            del block_members[b]
            del block_wires[b]
        members.append(i)
        members.sort()
        block_members[target] = members
        block_wires[target] = union
        for w in union:
            open_for_wire[w] = target

    # Sort by block *id* (creation order), not by members[0]: a block that
    # absorbs an earlier-anchored block keeps its own (later) creation
    # order for output purposes -- that id was exactly what earlier
    # merge/close decisions relied on staying ">= every touched block's id"
    # (see anchor_block = touched[-1] above).  Sorting by the post-merge
    # min member index instead would silently invalidate that invariant.
    return [block_members[bid] for bid in sorted(block_members)]


def compact_by_blocks(gates: list[GateInstance], max_wires: int = 3) -> list[GateInstance]:
    """Reorder ``gates`` so each <=max_wires-wire block becomes contiguous.

    A valid topological reordering of the circuit's dependency DAG (see
    :func:`collect_wire_blocks`), so it preserves the overall unitary exactly
    -- any two gates that swap position act on disjoint wires and therefore
    commute exactly.  This exposes reducible windows to the sweep that were
    previously scattered across unrelated-wire gates in the input order,
    without relying on the stochastic transport_shuffle to find them.
    """
    blocks = collect_wire_blocks(gates, max_wires)
    return [gates[i] for members in blocks for i in members]


def block_size_stats(gates: list[GateInstance], max_wires: int = 3) -> dict:
    """Diagnostic: block-length distribution for a compaction run."""
    blocks = collect_wire_blocks(gates, max_wires)
    sizes = sorted(len(b) for b in blocks)
    n = len(sizes)
    return {
        "num_blocks": n,
        "max_block_len": sizes[-1] if sizes else 0,
        "mean_block_len": (sum(sizes) / n) if n else 0.0,
        "median_block_len": sizes[n // 2] if n else 0,
    }
