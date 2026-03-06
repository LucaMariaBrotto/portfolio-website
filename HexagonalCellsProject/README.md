# Dijkstra Hex Grid Pathfinding

**C implementation of optimal pathfinding on a hexagonal grid with dynamic traversal costs and aerial routes.**

> Final project for *Algorithms and Data Structures* — A.Y. 2025/2026  
> Author: **Luca Maria Brotto**

---

## Overview

This project models a planetary surface as a rectangular grid of hexagonal tiles for the fictional transport company **Movhex**. Each tile has a traversal cost, and vehicles move between adjacent hexagons. The program computes optimal routes using **Dijkstra's algorithm** with a binary min-heap.

Additional **aerial routes** can dynamically connect any two non-adjacent tiles. Costs can change at runtime across radial regions of the map. A **result cache** avoids redundant Dijkstra runs, since `travel_cost` queries vastly outnumber map modifications in real workloads.

---

## Map Model

- The grid is a flat array of `hex` structs, indexed as `row * col_max + col`.
- Coordinates are given as `(column, row)`, both 0-indexed, left-to-right and bottom-to-top.
- Each tile has up to **6 ground neighbours**, computed at init time using offset coordinates (even/odd row alternation).
- Each tile has a **ground exit cost** in `[1, 100]`. A cost of `0` means the tile cannot be departed from (impassable), but can still be entered.
- Moving from tile A to an adjacent tile B costs the **exit cost of A**.
- **Aerial routes** are directed, stored as a singly linked list per tile (`testaporto`/`codaporto`). Each tile supports at most **5 outgoing aerial routes**. When an aerial route is used, the ground exit cost of its source tile is ignored.

---

## Commands

### `init <cols> <rows>`
Allocates and initialises the grid. All tiles start with cost `1`, no aerial routes.  
Re-initialising frees all previously allocated memory first.  
→ Responds `OK`.

---

### `change_cost <x> <y> <v> <radius>`
Applies a radial cost modification centred on `(x, y)` using a **BFS** up to `radius` steps.

For each tile at BFS depth `d < radius`:
```
delta = floor(v × (radius - d) / radius)
new_cost = clamp(old_cost + delta, 0, 100)
```

Both the tile's ground cost and all its outgoing aerial route costs are updated.  
`v` is an integer in `[−10, 10]`.  
→ Responds `KO` if coordinates are invalid or `radius ≤ 0`, otherwise `OK`.  
→ Invalidates the Dijkstra result cache.

---

### `toggle_air_route <x1> <y1> <x2> <y2>`
→ Responds `KO` immediately if either `(x1,y1)` or `(x2,y2)` are outside the grid bounds.

Otherwise, adds or removes a directed aerial route from `(x1,y1)` to `(x2,y2)`:

- **If adding**: the new route's cost is computed as `(sum of existing aerial route costs + ground exit cost of the tile) / num_volo`, where `num_volo` is the 1-based index of the new route (already incremented before the division). This is equivalent to a floor average over all routes including the tile's own ground cost.  
  → `OK` if fewer than 5 outgoing routes exist; otherwise `KO`.
- **If removing**: the route is deleted and the `num_volo` indices of remaining routes are decremented to stay contiguous.  
  → `OK`.

→ Invalidates the Dijkstra result cache.

---

### `travel_cost <xp> <yp> <xd> <yd>`
Returns the minimum cost to travel from `(xp,yp)` to `(xd,yd)`.

- Checks the **result cache** first; returns immediately on a hit.
- On a cache miss, runs **Dijkstra** with a binary min-heap over both ground and aerial edges.
- The exit cost of the **destination tile is not counted**.
- If source equals destination, cost is `0`.
- Returns `-1` if coordinates are invalid or no path exists.

---

## Architecture

The entire project is contained in a single file: **main.c**

| Category | Functions |
|----------|-----------|
| **Data Structures** | `hex`, `rotta` |
| **Grid Management** | `CalcoloCelVet`, `RiempiVicini`, `init`, `pulizia_vettore_porti` |
| **Cost Modification** | `change_cost` |
| **Air Route Management** | `calcolo_costo_volo`, `crea_porto`, `rimuovi_porto`, `toggle_air_route` |
| **Pathfinding** | `heap_push_or_decrease`, `heap_extract_min`, `heap_padre`, `heap_figlio`, `swapMinHeapNode`, `dijkstra` |
| **Cache & I/O** | `travel_cost`, `main` |

---

## Build & Run

**Requirements:** GCC (or any C99-compatible compiler), `libm` for `floor()`.

# Compile
gcc -O2 -o movhex main.c -lm

# Run with input file
./movhex < input.txt

# Or interactively
./movhex


---

## Example

```bash
init 100 100              → OK
change_cost 10 20 -10 5   → OK    (makes a region impassable)
travel_cost 0 0 20 0      → 20    (ground path cost)
toggle_air_route 0 0 20 0 → OK    (add aerial shortcut)
travel_cost 0 0 20 0      → 1     (aerial route is cheaper)
toggle_air_route 0 0 20 0 → OK    (remove aerial route)
travel_cost 0 0 20 0      → 20    (back to ground cost)
change_cost 200 20 -10 5  → KO    (invalid coordinates)
travel_cost 10 20 11 20   → -1    (impassable cell)
```
---

## Known Issues

- `toggle_air_route` prints `"KO \n"` (with a trailing space) instead of `"KO\n"` in error cases. This is inconsistent with all other commands and may cause failures with automated test systems that perform exact string matching.

## Performance Notes

The workload profile assumed by the specification (many `travel_cost` queries, few modifications) drives three key optimizations:

- **Result cache**: all `travel_cost` results are stored in a dynamic array. On a cache hit, the answer is returned in O(1) without touching the graph. The cache is fully invalidated on any `change_cost` or `toggle_air_route` call, and its capacity doubles automatically when full (starting at 500 entries).

- **BFS for `change_cost`**: BFS naturally visits cells in order of increasing distance, allowing early termination at `radius` depth and reusing a single pre-allocated queue.

- **Lazy BFS reset**: a global monotonic counter (`contatore_visitato`) marks visited cells without resetting the entire `visitato[]` array between calls, saving O(V) work per invocation.
