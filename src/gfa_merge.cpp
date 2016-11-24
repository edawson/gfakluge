/**
 * Merges a set of GFA files by:
 * 1. Coordinate their ID spaces (assume each one is an independent subgraph)
 * 2. Add both to a single GFAKluge instance. This will provide a sort.
 * 3. Push out the new GFA files to stdout.
 */

