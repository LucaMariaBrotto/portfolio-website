#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define INT_MAX 1000000

/*
 * Represents a single air route departing from a hex cell.
 * Routes are stored as a singly linked list on each cell.
 * route_num is the 1-based index of this route in the list,
 * used to compute the average flight cost across all routes.
 */
typedef struct Route {
    int dest_row;           // Row of the destination cell
    int dest_col;           // Column of the destination cell
    int8_t flight_cost;     // Cost of this air route
    int8_t route_num;       // Progressive index of this route in the list (1-based)
    struct Route *next;     // Pointer to the next route in the list
} route;

/*
 * Represents a single hexagonal cell of the map.
 * cost: traversal cost (range 1-100).
 *       A value of 0 means the cell cannot be left (impassable),
 *       but it can still be entered/visited.
 * neighbors[]: indices (in the flat array) of adjacent hex cells (up to 6).
 * route_head/route_tail: head and tail of the linked list of outgoing air routes,
 *                        kept to allow O(1) insertion at the tail.
 */
typedef struct {
    int8_t cost;
    int neighbors[6];
    int8_t num_neighbors;
    route *route_head;  // Head of the outgoing air routes list
    route *route_tail;  // Tail of the outgoing air routes list
} hex;

/*
 * Stores a cached result of a Dijkstra query.
 * Identified by source (r1,c1) and destination (r2,c2);
 * dist_cost holds the previously computed minimum cost.
 */
typedef struct PathCache {
    int r1, r2, c1, c2;
    int dist_cost;
} PathCache;

/* ---- Global variables ---- */

int V;                  // Total number of vertices (cells) in the map
int *minHeap = NULL;    // Min-heap array storing vertex indices
int *minHeapPos = NULL; // minHeapPos[v] = position of vertex v in minHeap, or -1 if absent
int minHeapSize = 0;    // Current number of elements in the heap

/* Dijkstra result cache */
int cache_count = 0;        // Number of valid entries currently in the cache
PathCache *cache;           // Array of cached (src, dst, cost) triples
int cache_capacity = 500;   // Current allocated capacity of the cache

/*
 * Cache validity flag.
 * Set to 0 whenever the map is modified (change_cost or toggle_air_route),
 * forcing the cache to be invalidated and rebuilt on the next travel_cost query.
 */
short cache_valid = 1;

int *dist = NULL; // dist[v] = current shortest distance from source to vertex v (used by Dijkstra)

/* -----------------------------------------------------------------------
 * cell_index
 * Converts a 2D (row, col) position into a flat 1D array index.
 * Used throughout the code to access the hex array by coordinates.
 * ----------------------------------------------------------------------- */
int cell_index(int r, int c, int col_max)
{
    return r * col_max + c;
}

/* -----------------------------------------------------------------------
 * fill_neighbors
 * Populates the neighbour list of every cell in the hex grid.
 * Hex grids use offset coordinates: neighbour offsets differ depending
 * on whether the row is even or odd, following the standard offset layout.
 * Only neighbours that fall within the grid bounds are recorded.
 * ----------------------------------------------------------------------- */
void fill_neighbors(hex *grid, int row_max, int col_max)
{
    for (int i = 0; i < row_max; i++)
    {
        for (int j = 0; j < col_max; j++)
        {
            int idx = cell_index(i, j, col_max);
            grid[idx].num_neighbors = 0;

            // Determine whether the current row is even or odd,
            // since hex neighbour offsets differ between the two cases.
            short is_even = (i % 2 == 0) ? 1 : 0;

            short row_offset[6];
            short col_offset[6];

            if (!is_even)
            {
                // Odd row neighbour offsets
                short col_offset2[6] = {1, 1, 0, -1, 0, 1};
                short row_offset2[6] = {0, 1, 1, 0, -1, -1};
                memcpy(row_offset, row_offset2, sizeof(row_offset2));
                memcpy(col_offset, col_offset2, sizeof(col_offset2));
            }
            else
            {
                // Even row neighbour offsets
                short col_offset2[6] = {1, 0, -1, -1, -1, 0};
                short row_offset2[6] = {0, 1, 1, 0, -1, -1};
                memcpy(row_offset, row_offset2, sizeof(row_offset2));
                memcpy(col_offset, col_offset2, sizeof(col_offset2));
            }

            // Check each of the 6 possible hex directions
            for (short k = 0; k < 6; k++)
            {
                int new_row = i + row_offset[k];
                int new_col = j + col_offset[k];

                // Add the neighbour only if it lies within the grid
                if (new_col >= 0 && new_col < col_max && new_row >= 0 && new_row < row_max)
                {
                    int neighbor_idx = cell_index(new_row, new_col, col_max);
                    grid[idx].neighbors[grid[idx].num_neighbors] = neighbor_idx;
                    grid[idx].num_neighbors += 1;
                }
            }
        }
    }
}

/* -----------------------------------------------------------------------
 * init
 * Initialises a freshly allocated hex grid:
 *   - sets every cell's traversal cost to 1 (minimum passable cost)
 *   - clears air route lists (route_head/route_tail = NULL)
 *   - calls fill_neighbors to build the adjacency structure
 * Prints "OK" on success.
 * ----------------------------------------------------------------------- */
void init(hex *grid, int row_max, int col_max)
{
    for (int i = 0; i < row_max; i++)
    {
        for (int j = 0; j < col_max; j++)
        {
            int idx = cell_index(i, j, col_max);
            grid[idx].cost = 1;
            grid[idx].route_head = NULL;
            grid[idx].route_tail = NULL;
            grid[idx].num_neighbors = 0;
        }
    }
    fill_neighbors(grid, row_max, col_max);
    printf("OK\n");
}

// Monotonically increasing counter used to mark visited cells during BFS
// without needing to reset the entire visited[] array between calls.
int visit_counter = 1;
int *visited = NULL; // visited[i] == visit_counter means cell i was visited in the current BFS
int *queue = NULL;   // Reusable BFS queue (allocated once at init, size = total cells)

/* -----------------------------------------------------------------------
 * change_cost
 * Applies a radial cost modification centred on cell (r, c).
 * The delta applied to each cell decreases linearly with its BFS depth:
 *   delta = floor(v * (radius - depth) / radius)
 * so the effect is strongest at the centre and fades to zero at the border.
 * Both cell traversal costs and outgoing air route costs are updated.
 * Valid cost range is [0, 100]; values outside are clamped.
 * Invalidates the Dijkstra cache after every successful modification.
 * Prints "OK" on success, "KO" if parameters are out of range.
 *
 * Parameters:
 *   v      - cost variation, in [-10, 10]
 *   radius - radius of effect in hex steps (must be > 0)
 * ----------------------------------------------------------------------- */
void change_cost(hex *grid, int c, int r, int row_max, int col_max, short v, int radius)
{
    if (c < 0 || r >= row_max || r < 0 || c >= col_max || radius <= 0 || v < -10 || v > 10)
    {
        printf("KO\n");
        return;
    }
    cache_valid = 0; // Invalidate the Dijkstra cache

    int head = 0;
    int queue_idx = 0;

    // Seed the BFS with the starting cell
    int start_cell = cell_index(r, c, col_max);
    queue[queue_idx] = start_cell;
    queue_idx = 1;
    visited[start_cell] = visit_counter;

    int depth = 0; // Current BFS depth (= distance from the centre cell)

    while (head < queue_idx) // BFS loop: process until the queue is empty
    {
        int level_size = queue_idx - head;

        // Process all nodes at the current BFS depth level
        for (int d = 0; d < level_size; d++)
        {
            int current = queue[head];
            head++;

            if (depth < radius)
            {
                // Compute the proportional delta for this depth level
                double cost_factor = (double)(radius - depth) / (double)radius;
                if (cost_factor <= 0.0) cost_factor = 0.0; // Guard against floating-point underflow
                short delta = (short)floor(v * cost_factor);

                // Apply delta to the cell's traversal cost, clamping to [0, 100]
                short new_cost = grid[current].cost + delta;
                if (new_cost < 1)
                    grid[current].cost = 0;
                else if (new_cost > 100)
                    grid[current].cost = 100;
                else
                    grid[current].cost = new_cost;

                // Apply the same delta to every outgoing air route of this cell
                route *current_route = grid[current].route_head;
                while (current_route != NULL)
                {
                    short new_flight_cost = current_route->flight_cost + delta;
                    if (new_flight_cost < 1)
                        current_route->flight_cost = 0;
                    else if (new_flight_cost > 100)
                        current_route->flight_cost = 100;
                    else
                        current_route->flight_cost = new_flight_cost;
                    current_route = current_route->next;
                }

                // Enqueue unvisited neighbours for the next depth level
                for (short i = 0; i < grid[current].num_neighbors; i++)
                {
                    int neighbor_idx = (int)grid[current].neighbors[i];
                    if (visited[neighbor_idx] != visit_counter)
                    {
                        visited[neighbor_idx] = visit_counter;
                        queue[queue_idx++] = neighbor_idx;
                    }
                }
            }
        }
        depth++;
    }

    // Advance the global BFS counter so visited[] is effectively reset for next call
    visit_counter++;
    printf("OK\n");
}

/* -----------------------------------------------------------------------
 * compute_flight_cost
 * Computes the cost to assign to a newly added air route from the given cell.
 * Formula: (sum of existing route costs + cell traversal cost) / total routes
 * This distributes the cell's own cost evenly among all its routes.
 * NOTE: route_tail must already point to the new (last) node before calling this.
 * ----------------------------------------------------------------------- */
short compute_flight_cost(hex *cell)
{
    short total = 0;
    route *curr = cell->route_head;

    // Sum up the costs of all currently stored routes
    while (curr != NULL)
    {
        total += curr->flight_cost;
        curr = curr->next;
    }
    // Average over total number of routes, including the cell's own traversal cost
    return ((total + cell->cost) / (cell->route_tail->route_num));
}

/* -----------------------------------------------------------------------
 * create_route
 * Adds a new outgoing air route from 'cell' to destination (r2, c2).
 * Handles two cases:
 *   1. First route ever added to this cell (route_head is NULL).
 *   2. Subsequent route: appended at the tail of the existing list.
 * After insertion the route cost is computed via compute_flight_cost().
 * ----------------------------------------------------------------------- */
void create_route(hex *cell, int c2, int r2)
{
    if (!cell->route_head)
    {
        // First route: initialise both head and tail pointers
        cell->route_head = malloc(sizeof(route));
        cell->route_head->dest_row = r2;
        cell->route_head->dest_col = c2;
        cell->route_head->next = NULL;
        cell->route_head->route_num = 1;
        cell->route_tail = cell->route_head;
        cell->route_head->flight_cost = 0;
        cell->route_head->flight_cost = compute_flight_cost(cell);
    }
    else
    {
        // Subsequent route: append to the tail of the list
        route *temp = malloc(sizeof(route));
        cell->route_tail->next = temp;
        temp->dest_row = r2;
        temp->dest_col = c2;
        temp->next = NULL;
        temp->route_num = cell->route_tail->route_num + 1; // Progressive index
        temp->flight_cost = 0;
        temp->flight_cost = compute_flight_cost(cell);
        cell->route_tail = temp; // Advance the tail pointer
    }
}

/* -----------------------------------------------------------------------
 * remove_route
 * Removes a specific route node from the linked list of 'cell' and frees it.
 * After removal, updates the head/tail pointers if necessary, then
 * decrements route_num for all remaining nodes to keep indices contiguous.
 * ----------------------------------------------------------------------- */
void remove_route(hex *cell, route *to_remove)
{
    route *curr = cell->route_head;
    route *prev = NULL;

    // Walk the list to find the node to remove
    while (curr != to_remove)
    {
        prev = curr;
        curr = curr->next;
    }

    if (prev == NULL)
    {
        // Removing the head node
        cell->route_head = curr->next;
        if (to_remove == cell->route_tail)
            cell->route_tail = cell->route_head; // List is now empty
    }
    else
    {
        prev->next = curr->next;
        if (to_remove == cell->route_tail)
            cell->route_tail = prev; // Update tail if we removed the last node
    }

    free(to_remove);

    // Renumber all remaining routes so route_num stays contiguous
    route *curr2 = cell->route_head;
    while (curr2 != NULL)
    {
        curr2->route_num -= 1;
        curr2 = curr2->next;
    }
}

/* -----------------------------------------------------------------------
 * toggle_air_route
 * Toggles an air route between cell (r1,c1) and cell (r2,c2):
 *   - If the route does NOT exist and the cell has fewer than 5 routes,
 *     it is created.
 *   - If the route already EXISTS, it is removed.
 *   - If the cell already has 5 routes and the route doesn't exist,
 *     prints "KO" (maximum routes per cell reached).
 * Validates that both cells are within grid bounds.
 * Invalidates the Dijkstra cache after any successful modification.
 * ----------------------------------------------------------------------- */
void toggle_air_route(hex *grid, int c1, int r1, int c2, int r2, int row_max, int col_max)
{
    if (c1 < 0 || r1 >= row_max || r1 < 0 || c1 >= col_max ||
        c2 < 0 || r2 >= row_max || r2 < 0 || c2 >= col_max)
    {
        printf("KO \n");
        return;
    }
    cache_valid = 0; // Invalidate the Dijkstra cache

    int cell_idx = cell_index(r1, c1, col_max);
    route *curr = grid[cell_idx].route_head;

    // No routes yet: just create the first one
    if (curr == NULL)
    {
        create_route(&grid[cell_idx], c2, r2);
        printf("OK\n");
        return;
    }

    // Search for an existing route to (r2, c2)
    while (curr->next != NULL && (curr->dest_row != r2 || curr->dest_col != c2))
        curr = curr->next;

    if (curr->dest_row == r2 && curr->dest_col == c2)
    {
        // Route found: remove it (toggle off)
        remove_route(&grid[cell_idx], curr);
        printf("OK\n");
    }
    else if (curr->next == NULL && grid[cell_idx].route_tail->route_num < 5)
    {
        // Route not found and limit not reached: add it (toggle on)
        create_route(&grid[cell_idx], c2, r2);
        printf("OK\n");
    }
    else
    {
        // Route not found but the 5-route limit is already reached
        printf("KO \n");
    }
}

/* ---- Min-heap operations (used by Dijkstra) ---- */

/* -----------------------------------------------------------------------
 * swapMinHeapNode
 * Swaps two elements in the min-heap array and updates their position map.
 * ----------------------------------------------------------------------- */
inline void swapMinHeapNode(int i, int j)
{
    int a = minHeap[i];
    int b = minHeap[j];
    minHeap[i] = b;
    minHeap[j] = a;
    minHeapPos[b] = i;
    minHeapPos[a] = j;
}

/* -----------------------------------------------------------------------
 * heap_sift_down
 * Restores the min-heap property by moving element at index i downward
 * until it is smaller than both its children.
 * ----------------------------------------------------------------------- */
void heap_sift_down(int i)
{
    while (1)
    {
        int l = 2 * i + 1, r = 2 * i + 2, smallest = i;
        if (l < minHeapSize && dist[minHeap[l]] < dist[minHeap[smallest]]) smallest = l;
        if (r < minHeapSize && dist[minHeap[r]] < dist[minHeap[smallest]]) smallest = r;
        if (smallest != i)
        {
            swapMinHeapNode(i, smallest);
            i = smallest;
        }
        else
            break;
    }
}

/* -----------------------------------------------------------------------
 * heap_sift_up
 * Restores the min-heap property by moving element at index i upward
 * until it is greater than or equal to its parent.
 * Used after inserting a new element or decreasing an existing key.
 * ----------------------------------------------------------------------- */
void heap_sift_up(int i)
{
    while (i > 0)
    {
        int p = (i - 1) / 2; // Parent index
        if (dist[minHeap[i]] < dist[minHeap[p]])
        {
            swapMinHeapNode(i, p);
            i = p;
        }
        else
            break;
    }
}

/* -----------------------------------------------------------------------
 * cleanup
 * Frees all dynamically allocated memory associated with the current map:
 *   - air route linked lists attached to every cell
 *   - the hex cell array itself
 *   - the min-heap arrays and the Dijkstra cache
 *   - the BFS auxiliary arrays (visited, queue)
 *   - the dist array used by Dijkstra
 * Called before re-initialising the map or at program exit.
 * ----------------------------------------------------------------------- */
void cleanup(hex *grid, int row_max, int col_max)
{
    // Free the air route list of every cell
    for (int i = 0; i < row_max * col_max; i++)
    {
        route *curr = grid[i].route_head;
        route *temp;
        while (curr != NULL)
        {
            temp = curr;
            curr = curr->next;
            free(temp);
        }
    }

    // Free heap arrays and the Dijkstra result cache
    if (minHeap)
    {
        free(minHeap);
        minHeap = NULL;
    }
    if (minHeapPos)
    {
        free(minHeapPos);
        minHeapPos = NULL;
        free(cache);
    }

    free(grid);
    free(dist);
    free(visited);
    free(queue);
}

/* -----------------------------------------------------------------------
 * heap_push_or_decrease
 * Inserts vertex v into the min-heap if it is not already present,
 * or triggers a sift-up if its distance has been decreased (decrease-key).
 * In both cases, heap_sift_up restores the heap property.
 * ----------------------------------------------------------------------- */
void heap_push_or_decrease(int v)
{
    if (minHeapPos[v] == -1)
    {
        // Vertex not yet in heap: append at the end and sift up
        minHeap[minHeapSize] = v;
        minHeapPos[v] = minHeapSize++;
        heap_sift_up(minHeapPos[v]);
    }
    else
    {
        // Vertex already in heap: its key was decreased, sift up from current position
        heap_sift_up(minHeapPos[v]);
    }
}

/* -----------------------------------------------------------------------
 * heap_extract_min
 * Removes and returns the vertex with the smallest dist[] value (the root).
 * The last element is moved to the root and sifted down to restore the heap.
 * Returns -1 if the heap is empty.
 * ----------------------------------------------------------------------- */
int heap_extract_min()
{
    if (minHeapSize == 0) return -1;

    int root = minHeap[0];
    minHeapSize--;

    if (minHeapSize > 0)
    {
        // Move last element to root position and restore heap order
        minHeap[0] = minHeap[minHeapSize];
        minHeapPos[minHeap[0]] = 0;
        heap_sift_down(0);
    }

    minHeapPos[root] = -1; // Mark as no longer in heap
    return root;
}

/* -----------------------------------------------------------------------
 * dijkstra
 * Computes the minimum travel cost from cell (r1,c1) to cell (r2,c2)
 * using Dijkstra's algorithm with a binary min-heap.
 * Both ground movement (through adjacent hex cells) and air routes
 * (stored per cell as linked lists) are considered as graph edges.
 * A cell with cost == 0 cannot be departed from (edges from it are ignored).
 * An air route with flight_cost == 0 is treated as unusable.
 * Returns the minimum cost, or -1 if the destination is unreachable.
 * ----------------------------------------------------------------------- */
int dijkstra(hex *grid, int r1, int c1, int r2, int c2, int row_max, int col_max)
{
    int src = cell_index(r1, c1, col_max);
    int dest = cell_index(r2, c2, col_max);

    // Initialise all distances to infinity and clear heap positions
    for (int v = 0; v < V; ++v)
    {
        dist[v] = INT_MAX;
        minHeapPos[v] = -1;
    }
    minHeapSize = 0;

    dist[src] = 0;
    heap_push_or_decrease(src);

    while (minHeapSize > 0)
    {
        int u = heap_extract_min();
        if (u == dest) break; // Early exit: shortest path to destination found

        // Relax ground edges (hex adjacency), skipping impassable cells
        if (grid[u].cost != 0)
        {
            for (int i = 0; i < grid[u].num_neighbors; i++)
            {
                int v = grid[u].neighbors[i];
                if (dist[u] != INT_MAX && dist[u] + grid[u].cost < dist[v])
                {
                    dist[v] = dist[u] + grid[u].cost;
                    heap_push_or_decrease(v);
                }
            }
        }

        // Relax air route edges departing from u
        route *current_route = grid[u].route_head;
        while (current_route != NULL)
        {
            int v = cell_index(current_route->dest_row, current_route->dest_col, col_max);
            int flight_cost = current_route->flight_cost;
            if (flight_cost != 0) // Skip disabled air routes
            {
                if (dist[u] != INT_MAX && dist[u] + flight_cost < dist[v])
                {
                    dist[v] = dist[u] + flight_cost;
                    heap_push_or_decrease(v);
                }
            }
            current_route = current_route->next;
        }
    }

    // Return the result, or -1 if the destination was never reached
    return (dist[dest] < INT_MAX) ? dist[dest] : -1;
}

/* -----------------------------------------------------------------------
 * travel_cost
 * Public interface for shortest-path queries between (r1,c1) and (r2,c2).
 * Implements a result cache: if the same query was computed since the last
 * map modification, the cached value is returned immediately without
 * re-running Dijkstra. The cache is invalidated (and reset) whenever
 * change_cost or toggle_air_route modifies the map.
 * If the cache is full, its capacity is doubled via realloc.
 * Prints the minimum cost, or -1 if the destination is unreachable
 * or the coordinates are out of bounds.
 * ----------------------------------------------------------------------- */
void travel_cost(hex *grid, int r1, int c1, int r2, int c2, int row_max, int col_max)
{
    if (r1 < 0 || r1 >= row_max || c1 < 0 || c1 >= col_max ||
        r2 < 0 || r2 >= row_max || c2 < 0 || c2 >= col_max)
    {
        printf("-1\n");
        return;
    }

    // Search the cache for a matching (src, dst) entry
    int i = 0;
    while (i != cache_count && (cache[i].r1 != r1 || cache[i].c1 != c1 ||
                                 cache[i].r2 != r2 || cache[i].c2 != c2))
        i++;

    if (!cache_valid)
    {
        // Map was modified since the last query: discard all cached results
        cache_count = 0;
    }
    else
    {
        if (i != cache_count)
        {
            // Cache hit: return the stored result directly
            printf("%d\n", cache[i].dist_cost);
            return;
        }
        else if (cache_count == cache_capacity)
        {
            // Cache is full: double its capacity
            cache_capacity *= 2;
            cache = realloc(cache, cache_capacity * sizeof(PathCache));
        }
    }

    // Cache miss (or cache was just invalidated): run Dijkstra and store result
    int cost = dijkstra(grid, r1, c1, r2, c2, row_max, col_max);

    cache[cache_count].r1 = r1;
    cache[cache_count].c1 = c1;
    cache[cache_count].r2 = r2;
    cache[cache_count].c2 = c2;
    cache[cache_count].dist_cost = cost;
    cache_count++;
    cache_valid = 1; // Cache is now valid again

    if (cost != -1)
        printf("%d\n", cost);
    else
        printf("-1\n");
}

/* -----------------------------------------------------------------------
 * main
 * Entry point. Reads commands from stdin in a loop until EOF:
 *   init <cols> <rows>                        - Allocates and initialises the hex grid
 *   change_cost <c> <r> <v> <radius>          - Radial cost modification
 *   travel_cost <c1> <r1> <c2> <r2>           - Shortest path query
 *   toggle_air_route <c1> <r1> <c2> <r2>      - Toggle air route
 * All dynamic memory is freed before the program exits.
 * ----------------------------------------------------------------------- */
int main(void)
{
    int col_max, row_max;
    hex *grid = NULL;
    int result;

    char command[64];
    while (scanf("%s", command) != EOF)
    {
        if (strcmp(command, "init") == 0)
        {
            // Free any previously allocated grid before re-initialising
            if (grid)
                cleanup(grid, row_max, col_max);

            result = scanf("%d %d", &col_max, &row_max);

            // Allocate all global data structures sized to the new grid
            grid = malloc(row_max * col_max * sizeof(hex));
            dist = malloc(row_max * col_max * sizeof(int));
            visited = calloc(row_max * col_max, sizeof(int));
            queue = malloc(row_max * col_max * sizeof(int));
            cache = malloc(500 * sizeof(PathCache));
            cache_capacity = 500;

            V = row_max * col_max;
            minHeap = malloc(V * sizeof(int));
            minHeapPos = malloc(V * sizeof(int));
            for (int i = 0; i < V; ++i) minHeapPos[i] = -1;
            minHeapSize = 0;

            init(grid, row_max, col_max);
            cache_valid = 0; // No cached results yet after a fresh init
        }
        else if (strcmp(command, "change_cost") == 0)
        {
            int r, c, v, radius;
            // Note: input order is col then row (c before r)
            result = scanf("%d %d %d %d ", &c, &r, &v, &radius);
            change_cost(grid, c, r, row_max, col_max, v, radius);
        }
        else if (strcmp(command, "travel_cost") == 0)
        {
            int r1, c1, r2, c2;
            // Note: input order is col then row for both endpoints
            result = scanf("%d %d %d %d", &c1, &r1, &c2, &r2);
            travel_cost(grid, r1, c1, r2, c2, row_max, col_max);
        }
        else if (strcmp(command, "toggle_air_route") == 0)
        {
            int c1, r1, c2, r2;
            result = scanf("%d %d %d %d", &c1, &r1, &c2, &r2);
            toggle_air_route(grid, c1, r1, c2, r2, row_max, col_max);
        }
    }

    // Free all memory before exit
    cleanup(grid, row_max, col_max);
    (void)result;
    return 0;
}
