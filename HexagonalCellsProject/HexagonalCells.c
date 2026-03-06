#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define INT_MAX 1000000

/*
 * Represents a single air route departing from a hex cell.
 * Routes are stored as a singly linked list on each cell.
 * num_volo is the 1-based index of this route in the list,
 * used to compute the average flight cost across all routes.
 */
typedef struct Rotta {
    int riga_arr;       // Row of the destination cell
    int col_arr;        // Column of the destination cell
    int8_t costo_volo;  // Cost of this air route
    int8_t num_volo;    // Progressive index of this route in the list (1-based)
    struct Rotta *next; // Pointer to the next route in the list
} rotta;

/*
 * Represents a single hexagonal cell of the map.
 * costo: traversal cost (range 1-100).
 *        A value of 0 means the cell cannot be left (impassable),
 *        but it can still be entered/visited.
 * vicini[]: indices (in the flat array) of adjacent hex cells (up to 6).
 * testaporto/codaporto: head and tail of the linked list of outgoing air routes,
 *                       kept to allow O(1) insertion at the tail.
 */
typedef struct {
    int8_t costo;
    int vicini[6];
    int8_t num_vicini;
    rotta *testaporto; // Head of the outgoing air routes list
    rotta *codaporto;  // Tail of the outgoing air routes list
} hex;

/*
 * Stores a cached result of a Dijkstra query.
 * Identified by source (r1,c1) and destination (r2,c2);
 * costo_dist holds the previously computed minimum cost.
 */
typedef struct costo_tratta {
    int r1, r2, c1, c2;
    int costo_dist;
} costo_tratta;

/* ---- Global variables ---- */

int V;               // Total number of vertices (cells) in the map
int *minHeap = NULL;    // Min-heap array storing vertex indices
int *minHeapPos = NULL; // minHeapPos[v] = position of vertex v in minHeap, or -1 if absent
int minHeapSize = 0;    // Current number of elements in the heap

/* Dijkstra result cache */
int cont_cache = 0;     // Number of valid entries currently in the cache
costo_tratta *cache;    // Array of cached (src, dst, cost) triples
int cache_spazio = 500; // Current allocated capacity of the cache

/*
 * Cache validity flag.
 * Set to 0 whenever the map is modified (change_cost or toggle_air_route),
 * forcing the cache to be invalidated and rebuilt on the next travel_cost query.
 */
short valido = 1;

int *dist = NULL; // dist[v] = current shortest distance from source to vertex v (used by Dijkstra)

/* -----------------------------------------------------------------------
 * CalcoloCelVet
 * Converts a 2D (row, col) position into a flat 1D array index.
 * Used throughout the code to access the hex array by coordinates.
 * ----------------------------------------------------------------------- */
int CalcoloCelVet(int r, int c, int col_max)
{
    return r * col_max + c;
}

/* -----------------------------------------------------------------------
 * RiempiVicini
 * Populates the neighbour list of every cell in the hex grid.
 * Hex grids use offset coordinates: neighbour offsets differ depending
 * on whether the row is even or odd, following the standard offset layout.
 * Only neighbours that fall within the grid bounds are recorded.
 * ----------------------------------------------------------------------- */
void RiempiVicini(hex *vettore, int riga_max, int col_max)
{
    for (int i = 0; i < riga_max; i++)
    {
        for (int j = 0; j < col_max; j++)
        {
            int CelVet = CalcoloCelVet(i, j, col_max);
            vettore[CelVet].num_vicini = 0;

            // Determine whether the current row is even or odd,
            // since hex neighbour offsets differ between the two cases.
            short pari = (i % 2 == 0) ? 1 : 0;

            short distriga[6];
            short distcol[6];

            if (!pari)
            {
                // Odd row neighbour offsets
                short distcol2[6] = {1, 1, 0, -1, 0, 1};
                short distriga2[6] = {0, 1, 1, 0, -1, -1};
                memcpy(distriga, distriga2, sizeof(distriga2));
                memcpy(distcol, distcol2, sizeof(distcol2));
            }
            else
            {
                // Even row neighbour offsets
                short distcol2[6] = {1, 0, -1, -1, -1, 0};
                short distriga2[6] = {0, 1, 1, 0, -1, -1};
                memcpy(distriga, distriga2, sizeof(distriga2));
                memcpy(distcol, distcol2, sizeof(distcol2));
            }

            // Check each of the 6 possible hex directions
            for (short k = 0; k < 6; k++)
            {
                int nuovoi = i + distriga[k];
                int nuovoj = j + distcol[k];

                // Add the neighbour only if it lies within the grid
                if (nuovoj >= 0 && nuovoj < col_max && nuovoi >= 0 && nuovoi < riga_max)
                {
                    int CelVetVic = CalcoloCelVet(nuovoi, nuovoj, col_max);
                    vettore[CelVet].vicini[vettore[CelVet].num_vicini] = CelVetVic;
                    vettore[CelVet].num_vicini += 1;
                }
            }
        }
    }
}

/* -----------------------------------------------------------------------
 * init
 * Initialises a freshly allocated hex grid:
 *   - sets every cell's traversal cost to 1 (minimum passable cost)
 *   - clears air route lists (testaporto/codaporto = NULL)
 *   - calls RiempiVicini to build the adjacency structure
 * Prints "OK" on success.
 * ----------------------------------------------------------------------- */
void init(hex *vettore, int riga_max, int col_max)
{
    for (int i = 0; i < riga_max; i++)
    {
        for (int j = 0; j < col_max; j++)
        {
            int CelVet = CalcoloCelVet(i, j, col_max);
            vettore[CelVet].costo = 1;
            vettore[CelVet].testaporto = NULL;
            vettore[CelVet].codaporto = NULL;
            vettore[CelVet].num_vicini = 0;
        }
    }
    RiempiVicini(vettore, riga_max, col_max);
    printf("OK\n");
}

// Monotonically increasing counter used to mark visited cells during BFS
// without needing to reset the entire visitato[] array between calls.
int contatore_visitato = 1;
int *visitato = NULL; // visitato[i] == contatore_visitato means cell i was visited in the current BFS
int *coda = NULL;     // Reusable BFS queue (allocated once at init, size = total cells)

/* -----------------------------------------------------------------------
 * change_cost
 * Applies a radial cost modification centred on cell (r, c).
 * The delta applied to each cell decreases linearly with its BFS depth:
 *   delta = floor(v * (raggio - depth) / raggio)
 * so the effect is strongest at the centre and fades to zero at the border.
 * Both cell traversal costs and outgoing air route costs are updated.
 * Valid cost range is [0, 100]; values outside are clamped.
 * Invalidates the Dijkstra cache after every successful modification.
 * Prints "OK" on success, "KO" if parameters are out of range.
 *
 * Parameters:
 *   v      - cost variation, in [-10, 10]
 *   raggio - radius of effect in hex steps (must be > 0)
 * ----------------------------------------------------------------------- */
void change_cost(hex *vettore, int c, int r, int riga_max, int col_max, short v, int raggio)
{
    if (c < 0 || r >= riga_max || r < 0 || c >= col_max || raggio <= 0 || v < -10 || v > 10)
    {
        printf("KO\n");
        return;
    }
    valido = 0; // Invalidate the Dijkstra cache

    int testa = 0;
    int coda_idx = 0;

    // Seed the BFS with the starting cell
    int cella_di_partenza = CalcoloCelVet(r, c, col_max);
    coda[coda_idx] = cella_di_partenza;
    coda_idx = 1;
    visitato[cella_di_partenza] = contatore_visitato;

    int profondita = 0; // Current BFS depth (= distance from the centre cell)

    while (testa < coda_idx) // BFS loop: process until the queue is empty
    {
        int dimensione_coda = coda_idx - testa;

        // Process all nodes at the current BFS depth level
        for (int d = 0; d < dimensione_coda; d++)
        {
            int corrente = coda[testa];
            testa++;

            if (profondita < raggio)
            {
                // Compute the proportional delta for this depth level
                double costof = (double)(raggio - profondita) / (double)raggio;
                if (costof <= 0.0) costof = 0.0; // Guard against floating-point underflow
                short delta = (short)floor(v * costof);

                // Apply delta to the cell's traversal cost, clamping to [0, 100]
                short nuovo_costo = vettore[corrente].costo + delta;
                if (nuovo_costo < 1)
                    vettore[corrente].costo = 0;
                else if (nuovo_costo > 100)
                    vettore[corrente].costo = 100;
                else
                    vettore[corrente].costo = nuovo_costo;

                // Apply the same delta to every outgoing air route of this cell
                rotta *current_route = vettore[corrente].testaporto;
                while (current_route != NULL)
                {
                    short nuovo_costo_volo = current_route->costo_volo + delta;
                    if (nuovo_costo_volo < 1)
                        current_route->costo_volo = 0;
                    else if (nuovo_costo_volo > 100)
                        current_route->costo_volo = 100;
                    else
                        current_route->costo_volo = nuovo_costo_volo;
                    current_route = current_route->next;
                }

                // Enqueue unvisited neighbours for the next depth level
                for (short i = 0; i < vettore[corrente].num_vicini; i++)
                {
                    int idx_vicino = (int)vettore[corrente].vicini[i];
                    if (visitato[idx_vicino] != contatore_visitato)
                    {
                        visitato[idx_vicino] = contatore_visitato;
                        coda[coda_idx++] = idx_vicino;
                    }
                }
            }
        }
        profondita++;
    }

    // Advance the global BFS counter so visitato[] is effectively reset for next call
    contatore_visitato++;
    printf("OK\n");
}

/* -----------------------------------------------------------------------
 * calcolo_costo_volo
 * Computes the cost to assign to a newly added air route from the given cell.
 * Formula: (sum of existing route costs + cell traversal cost) / total routes
 * This distributes the cell's own cost evenly among all its routes.
 * NOTE: codaporto must already point to the new (last) node before calling this.
 * ----------------------------------------------------------------------- */
short calcolo_costo_volo(hex *cella)
{
    short complessivo = 0;
    rotta *appoggio = cella->testaporto;

    // Sum up the costs of all currently stored routes
    while (appoggio != NULL)
    {
        complessivo += appoggio->costo_volo;
        appoggio = appoggio->next;
    }
    // Average over total number of routes, including the cell's own traversal cost
    return ((complessivo + cella->costo) / (cella->codaporto->num_volo));
}

/* -----------------------------------------------------------------------
 * crea_porto
 * Adds a new outgoing air route from 'cella' to destination (r2, c2).
 * Handles two cases:
 *   1. First route ever added to this cell (testaporto is NULL).
 *   2. Subsequent route: appended at the tail of the existing list.
 * After insertion the route cost is computed via calcolo_costo_volo().
 * ----------------------------------------------------------------------- */
void crea_porto(hex *cella, int c2, int r2)
{
    if (!cella->testaporto)
    {
        // First route: initialise both head and tail pointers
        cella->testaporto = malloc(sizeof(rotta));
        cella->testaporto->riga_arr = r2;
        cella->testaporto->col_arr = c2;
        cella->testaporto->next = NULL;
        cella->testaporto->num_volo = 1;
        cella->codaporto = cella->testaporto;
        cella->testaporto->costo_volo = 0;
        cella->testaporto->costo_volo = calcolo_costo_volo(cella);
    }
    else
    {
        // Subsequent route: append to the tail of the list
        rotta *appoggio = malloc(sizeof(rotta));
        cella->codaporto->next = appoggio;
        appoggio->riga_arr = r2;
        appoggio->col_arr = c2;
        appoggio->next = NULL;
        appoggio->num_volo = cella->codaporto->num_volo + 1; // Progressive index
        appoggio->costo_volo = 0;
        appoggio->costo_volo = calcolo_costo_volo(cella);
        cella->codaporto = appoggio; // Advance the tail pointer
    }
}

/* -----------------------------------------------------------------------
 * rimuovi_porto
 * Removes a specific route node from the linked list of 'cella' and frees it.
 * After removal, updates the head/tail pointers if necessary, then
 * decrements num_volo for all remaining nodes to keep indices contiguous.
 * ----------------------------------------------------------------------- */
void rimuovi_porto(hex *cella, rotta *dacancellare)
{
    rotta *appoggio = cella->testaporto;
    rotta *prec = NULL;

    // Walk the list to find the node to remove
    while (appoggio != dacancellare)
    {
        prec = appoggio;
        appoggio = appoggio->next;
    }

    if (prec == NULL)
    {
        // Removing the head node
        cella->testaporto = appoggio->next;
        if (dacancellare == cella->codaporto)
            cella->codaporto = cella->testaporto; // List is now empty
    }
    else
    {
        prec->next = appoggio->next;
        if (dacancellare == cella->codaporto)
            cella->codaporto = prec; // Update tail if we removed the last node
    }

    free(dacancellare);

    // Renumber all remaining routes so num_volo stays contiguous
    rotta *appoggio2 = cella->testaporto;
    while (appoggio2 != NULL)
    {
        appoggio2->num_volo -= 1;
        appoggio2 = appoggio2->next;
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
void toggle_air_route(hex *vettore, int c1, int r1, int c2, int r2, int riga_max, int col_max)
{
    if (c1 < 0 || r1 >= riga_max || r1 < 0 || c1 >= col_max ||
        c2 < 0 || r2 >= riga_max || r2 < 0 || c2 >= col_max)
    {
        printf("KO \n");
        return;
    }
    valido = 0; // Invalidate the Dijkstra cache

    int num_cella = CalcoloCelVet(r1, c1, col_max);
    rotta *appoggio = vettore[num_cella].testaporto;

    // No routes yet: just create the first one
    if (appoggio == NULL)
    {
        crea_porto(&vettore[num_cella], c2, r2);
        printf("OK\n");
        return;
    }

    // Search for an existing route to (r2, c2)
    while (appoggio->next != NULL && (appoggio->riga_arr != r2 || appoggio->col_arr != c2))
        appoggio = appoggio->next;

    if (appoggio->riga_arr == r2 && appoggio->col_arr == c2)
    {
        // Route found: remove it (toggle off)
        rimuovi_porto(&vettore[num_cella], appoggio);
        printf("OK\n");
    }
    else if (appoggio->next == NULL && vettore[num_cella].codaporto->num_volo < 5)
    {
        // Route not found and limit not reached: add it (toggle on)
        crea_porto(&vettore[num_cella], c2, r2);
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
 * heap_figlio (sift-down)
 * Restores the min-heap property by moving element at index i downward
 * until it is smaller than both its children.
 * ----------------------------------------------------------------------- */
void heap_figlio(int i)
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
 * heap_padre (sift-up)
 * Restores the min-heap property by moving element at index i upward
 * until it is greater than or equal to its parent.
 * Used after inserting a new element or decreasing an existing key.
 * ----------------------------------------------------------------------- */
void heap_padre(int i)
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
 * pulizia_vettore_porti
 * Frees all dynamically allocated memory associated with the current map:
 *   - air route linked lists attached to every cell
 *   - the hex cell array itself
 *   - the min-heap arrays and the Dijkstra cache
 *   - the BFS auxiliary arrays (visitato, coda)
 *   - the dist array used by Dijkstra
 * Called before re-initialising the map or at program exit.
 * ----------------------------------------------------------------------- */
void pulizia_vettore_porti(hex *vettore, int riga_max, int col_max)
{
    // Free the air route list of every cell
    for (int i = 0; i < riga_max * col_max; i++)
    {
        rotta *appoggio = vettore[i].testaporto;
        rotta *temp;
        while (appoggio != NULL)
        {
            temp = appoggio;
            appoggio = appoggio->next;
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

    free(vettore);
    free(dist);
    free(visitato);
    free(coda);
}

/* -----------------------------------------------------------------------
 * heap_push_or_decrease
 * Inserts vertex v into the min-heap if it is not already present,
 * or triggers a sift-up if its distance has been decreased (decrease-key).
 * In both cases, heap_padre restores the heap property.
 * ----------------------------------------------------------------------- */
void heap_push_or_decrease(int v)
{
    if (minHeapPos[v] == -1)
    {
        // Vertex not yet in heap: append at the end and sift up
        minHeap[minHeapSize] = v;
        minHeapPos[v] = minHeapSize++;
        heap_padre(minHeapPos[v]);
    }
    else
    {
        // Vertex already in heap: its key was decreased, sift up from current position
        heap_padre(minHeapPos[v]);
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
        heap_figlio(0);
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
 * A cell with costo == 0 cannot be departed from (edges from it are ignored).
 * An air route with costo_volo == 0 is treated as unusable.
 * Returns the minimum cost, or -1 if the destination is unreachable.
 * ----------------------------------------------------------------------- */
int dijkstra(hex *vettore, int r1, int c1, int r2, int c2, int riga_max, int col_max)
{
    int src = CalcoloCelVet(r1, c1, col_max);
    int dest = CalcoloCelVet(r2, c2, col_max);

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
        if (vettore[u].costo != 0)
        {
            for (int i = 0; i < vettore[u].num_vicini; i++)
            {
                int v = vettore[u].vicini[i];
                if (dist[u] != INT_MAX && dist[u] + vettore[u].costo < dist[v])
                {
                    dist[v] = dist[u] + vettore[u].costo;
                    heap_push_or_decrease(v);
                }
            }
        }

        // Relax air route edges departing from u
        rotta *current_route = vettore[u].testaporto;
        while (current_route != NULL)
        {
            int v = CalcoloCelVet(current_route->riga_arr, current_route->col_arr, col_max);
            int costo_volo = current_route->costo_volo;
            if (costo_volo != 0) // Skip disabled air routes
            {
                if (dist[u] != INT_MAX && dist[u] + costo_volo < dist[v])
                {
                    dist[v] = dist[u] + costo_volo;
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
void travel_cost(hex *vettore, int r1, int c1, int r2, int c2, int riga_max, int col_max)
{
    if (r1 < 0 || r1 >= riga_max || c1 < 0 || c1 >= col_max ||
        r2 < 0 || r2 >= riga_max || c2 < 0 || c2 >= col_max)
    {
        printf("-1\n");
        return;
    }

    // Search the cache for a matching (src, dst) entry
    int i = 0;
    while (i != cont_cache && (cache[i].r1 != r1 || cache[i].c1 != c1 ||
                                cache[i].r2 != r2 || cache[i].c2 != c2))
        i++;

    if (!valido)
    {
        // Map was modified since the last query: discard all cached results
        cont_cache = 0;
    }
    else
    {
        if (i != cont_cache)
        {
            // Cache hit: return the stored result directly
            printf("%d\n", cache[i].costo_dist);
            return;
        }
        else if (cont_cache == cache_spazio)
        {
            // Cache is full: double its capacity
            cache_spazio *= 2;
            cache = realloc(cache, cache_spazio * sizeof(costo_tratta));
        }
    }

    // Cache miss (or cache was just invalidated): run Dijkstra and store result
    int costo = dijkstra(vettore, r1, c1, r2, c2, riga_max, col_max);

    cache[cont_cache].r1 = r1;
    cache[cont_cache].c1 = c1;
    cache[cont_cache].r2 = r2;
    cache[cont_cache].c2 = c2;
    cache[cont_cache].costo_dist = costo;
    cont_cache++;
    valido = 1; // Cache is now valid again

    if (costo != -1)
        printf("%d\n", costo);
    else
        printf("-1\n");
}

/* -----------------------------------------------------------------------
 * main
 * Entry point. Reads commands from stdin in a loop until EOF:
 *   init <cols> <rows>         - Allocates and initialises the hex grid
 *   change_cost <c> <r> <v> <raggio> - Radial cost modification
 *   travel_cost <c1> <r1> <c2> <r2>  - Shortest path query
 *   toggle_air_route <c1> <r1> <c2> <r2> - Toggle air route
 * All dynamic memory is freed before the program exits.
 * ----------------------------------------------------------------------- */
int main(void)
{
    int col_max, riga_max;
    hex *vettore = NULL;
    int caso;

    char comando[64];
    while (scanf("%s", comando) != EOF)
    {
        if (strcmp(comando, "init") == 0)
        {
            // Free any previously allocated grid before re-initialising
            if (vettore)
                pulizia_vettore_porti(vettore, riga_max, col_max);

            caso = scanf("%d %d", &col_max, &riga_max);

            // Allocate all global data structures sized to the new grid
            vettore = malloc(riga_max * col_max * sizeof(hex));
            dist = malloc(riga_max * col_max * sizeof(int));
            visitato = calloc(riga_max * col_max, sizeof(int));
            coda = malloc(riga_max * col_max * sizeof(int));
            cache = malloc(500 * sizeof(costo_tratta));
            cache_spazio = 500;

            V = riga_max * col_max;
            minHeap = malloc(V * sizeof(int));
            minHeapPos = malloc(V * sizeof(int));
            for (int i = 0; i < V; ++i) minHeapPos[i] = -1;
            minHeapSize = 0;

            init(vettore, riga_max, col_max);
            valido = 0; // No cached results yet after a fresh init
        }
        else if (strcmp(comando, "change_cost") == 0)
        {
            int r, c, v, raggio;
            // Note: input order is col then row (c before r)
            caso = scanf("%d %d %d %d ", &c, &r, &v, &raggio);
            change_cost(vettore, c, r, riga_max, col_max, v, raggio);
        }
        else if (strcmp(comando, "travel_cost") == 0)
        {
            int r1, c1, r2, c2;
            // Note: input order is col then row for both endpoints
            caso = scanf("%d %d %d %d", &c1, &r1, &c2, &r2);
            travel_cost(vettore, r1, c1, r2, c2, riga_max, col_max);
        }
        else if (strcmp(comando, "toggle_air_route") == 0)
        {
            int c1, r1, c2, r2;
            caso = scanf("%d %d %d %d", &c1, &r1, &c2, &r2);
            toggle_air_route(vettore, c1, r1, c2, r2, riga_max, col_max);
        }
    }

    // Free all memory before exit
    pulizia_vettore_porti(vettore, riga_max, col_max);
    (void)caso;
    return 0;
}