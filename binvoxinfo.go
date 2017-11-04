//
// Written by Maxim Khitrov (November 2017)
//

package main

import (
	"bufio"
	"fmt"
	"io"
	"math"
	"math/bits"
	"os"
	"runtime"
	"sync"
	"sync/atomic"
)

const (
	maxDim = 1280          // Maximum grid dimensions (maxDim^3 <= MaxInt)
	sig    = "#binvox 1\n" // Binvox file signature
)

func main() {
	if len(os.Args) != 2 {
		fmt.Fprintf(os.Stderr, "usage: %s input.binvox\n", os.Args[0])
		os.Exit(2)
	}
	m, err := LoadModel(os.Args[1])
	if err != nil {
		fmt.Fprintf(os.Stderr, "failed to load model: %v\n", err)
		os.Exit(1)
	}

	dist := m.Dist()
	m.VerifyEdges()
	e0, e1, e2, e3 := m.edgeN[0], m.edgeN[1], m.edgeN[2], m.edgeN[3]
	nv := int64(m.count)
	ne := int64(len(m.edge)) - e0
	//m.Save("model.binvox")

	fmt.Printf("grid:        %d x %d x %d = %d\n", m.d, m.h, m.w, len(m.grid))
	fmt.Printf("translation: %v\n", m.tx)
	fmt.Printf("scale:       %f\n", m.scale)
	fmt.Printf("voxels:      %d (%.3f%%)\n", nv, pct(nv, int64(len(m.grid))))
	fmt.Printf("  orphans:   %d (%.1f%%)\n", e0, pct(e0, nv))
	fmt.Printf("  size:      %f (%f scaled)\n", m.vs, m.vs*m.scale)
	fmt.Printf("edges:       %d (%.1f%%)\n", ne, pct(ne, nv))
	fmt.Printf("  1d:        %d (%.1f%%)\n", e1, pct(e1, ne))
	fmt.Printf("  2d:        %d (%.1f%%)\n", e2, pct(e2, ne))
	fmt.Printf("  3d:        %d (%.1f%%)\n", e3, pct(e3, ne))
	fmt.Printf("  removed:   %d (%.1f%%)\n", m.edgeR, pct(m.edgeR, ne+m.edgeR))
	fmt.Printf("  avg epv:   %f\n", m.epv)
	fmt.Printf("  distance:  %f (%f scaled)\n", dist, dist*m.scale)
}

// Edge represents adjacency between two voxels i and j using their linear
// indices in the model grid. The format is:
//
//	i(40) | dx(8) | dz(8) | dy(8)
//
// Where i is the lesser of i and j, and dx,dz,dy are int8 xyz deltas from i to
// j. The current implementation requires all deltas to be in the range [-1,1].
// Fields dy and dz are swapped to match binvox format.
type Edge uint64

// NewEdge creates a new edge from voxel index and xyz deltas.
func NewEdge(i, dx, dy, dz int) Edge {
	return Edge(i)<<24 | Edge(byte(dx))<<16 | Edge(byte(dz))<<8 | Edge(byte(dy))
}

// Info decodes edge information and returns individual fields.
func (e Edge) Info() (i, dx, dy, dz int) {
	return int(e >> 24), int(int8(e >> 16)), int(int8(e)), int(int8(e >> 8))
}

// Dim returns the number of non-zero delta values in e between 0 and 3:
//
//	0: orphan voxel with an edge back to itself (distance = 0)
//	1: face-adjacent voxels (distance = voxelSize)
//	2: edge-adjacent voxels (distance = sqrt(2*voxelSize^2))
//	3: corner-adjacent voxels (distance = sqrt(3*voxelSize^2))
func (e Edge) Dim() int {
	e &= 0x010101
	return bits.OnesCount8(byte(e>>14 | e>>7 | e))
}

// Model is a 3D volumetric representation of an object. Voxels have x,y,z
// coordinates, which are converted to a linear index in the grid bitmap
// according to the binvox file format.
type Model struct {
	grid    []bool            // Voxel bitmap
	edge    map[Edge]struct{} // Adjacency map
	edgeN   [4]int64          // Edge count by dimension
	edgeR   int64             // Number of removed edges by simplification
	d, h, w int               // Grid depth, height, and width
	wh      int               // Cached width * height
	count   int               // Voxel count
	tx      [3]float64        // Model translation
	scale   float64           // Model scale
	vs      float64           // Voxel size
	epv     float64           // Average number of edges per voxel
}

// LoadModel loads the model from a binvox file.
func LoadModel(name string) (*Model, error) {
	f, err := os.Open(name)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	// Verify signature
	r := bufio.NewReader(f)
	if _, err = fmt.Fscanf(r, sig); err != nil {
		return nil, err
	}

	// Read grid dimensions
	var d, h, w uint16
	if _, err = fmt.Fscanf(r, "dim %d %d %d\n", &d, &h, &w); err != nil {
		return nil, err
	}
	if d > maxDim || h > maxDim || w > maxDim {
		return nil, fmt.Errorf("invalid dimensions: %d x %d x %d", d, h, w)
	}
	m := &Model{
		grid: make([]bool, int(d)*int(h)*int(w)),
		d:    int(d),
		h:    int(h),
		w:    int(w),
		wh:   int(w) * int(h),
		vs:   1.0 / float64(maxUint16(maxUint16(d, h), w)),
	}

	// Read rest of header
	_, err = fmt.Fscanf(r, "translate %f %f %f\n", &m.tx[0], &m.tx[1], &m.tx[2])
	if err != nil {
		return nil, err
	}
	if _, err = fmt.Fscanf(r, "scale %f\n", &m.scale); err != nil {
		return nil, err
	}
	if _, err = fmt.Fscanf(r, "data\n"); err != nil {
		return nil, err
	}

	// Read grid data (format is bool | count, 1 byte each)
	for i, b := 0, byte(0); i < len(m.grid); {
		if b, err = r.ReadByte(); err != nil {
			return nil, err
		}
		set := b != 0
		if b, err = r.ReadByte(); err != nil {
			return nil, err
		}
		if i += int(b); set {
			for j := i - int(b); j < i; j++ {
				m.grid[j] = true
			}
			m.count += int(b)
		}
	}
	if _, err = r.ReadByte(); err != io.EOF {
		if err == nil {
			err = fmt.Errorf("unexpected data past end of grid")
		}
		return nil, err
	}
	return m, nil
}

// Save writes the model to a binvox file. The input and output files should be
// identical except for normalization of grid encoding.
func (m *Model) Save(name string) error {
	f, err := os.Create(name)
	if err != nil {
		return err
	}
	defer f.Close()

	// Write header
	w := bufio.NewWriter(f)
	w.WriteString(sig)
	fmt.Fprintf(w, "dim %d %d %d\n", m.d, m.h, m.w)
	fmt.Fprintf(w, "translate %.4f %.4f %.4f\n", m.tx[0], m.tx[1], m.tx[2])
	fmt.Fprintf(w, "scale %.4f\n", m.scale)
	w.WriteString("data\n")

	// Write grid data
	b, n := false, byte(0)
	writeRun := func() {
		if n > 0 {
			if b {
				w.WriteByte(1)
			} else {
				w.WriteByte(0)
			}
			w.WriteByte(n)
		}
	}
	for _, v := range m.grid {
		if b == v && n < math.MaxUint8 {
			n++
		} else {
			writeRun()
			b = v
			n = 1
		}
	}
	writeRun()
	return w.Flush()
}

// Dist returns the total distance of all the edges in the model. This is only
// useful if the model was thinned until each voxel was left with two neighbors
// on average (epv ~2). Multiply this distance by the model scale to convert it
// to original units.
func (m *Model) Dist() float64 {
	if m.edge == nil {
		m.findEdges()
	}
	vs2 := m.vs * m.vs
	return float64(m.edgeN[1])*m.vs + // Face-adjacent
		float64(m.edgeN[2])*math.Sqrt(2*vs2) + // Edge-adjacent
		float64(m.edgeN[3])*math.Sqrt(3*vs2) // Corner-adjacent
}

// VerifyEdges reconstructs the grid from the adjacency map and compares it to
// the original. Errors are reported to stdout.
func (m *Model) VerifyEdges() {
	grid := make([]bool, len(m.grid))
	for e := range m.edge {
		i, dx, dy, dz := e.Info()
		j := i + m.idx(dx, dy, dz)
		grid[i], grid[j] = true, true
	}
	const errLimit = 100
	errCount := uint32(0)
	m.splitGrid(func(orig []bool, off int) {
		for i, v := range orig {
			if i += off; grid[i] != v {
				fmt.Printf("grid mismatch at %d: should be %t\n", i, v)
				if atomic.AddUint32(&errCount, 1) >= errLimit {
					break
				}
			}
		}
	})
	if errCount >= errLimit {
		fmt.Println("too many errors, comparison incomplete")
	}
}

// findEdges constructs the adjacency map for all voxels in the grid.
func (m *Model) findEdges() {
	m.edge = nil
	m.edgeN = [4]int64{}
	m.edgeR = 0

	// Start a goroutine for combining results
	type result struct {
		edge    map[Edge]struct{} // Set of all valid edges
		total   int64             // Sum of all possible edges for each voxel
		removed int64             // Number of simplified edges
	}
	ch, epv := make(chan result), make(chan float64)
	go func() {
		defer close(epv)
		c := result{edge: make(map[Edge]struct{})}
		for r := range ch {
			for e := range r.edge {
				c.edge[e] = struct{}{}
				m.edgeN[e.Dim()]++
			}
			c.total += r.total
			m.edgeR += r.removed
		}
		if int64(len(c.edge)) != m.edgeN[0]+m.edgeN[1]+m.edgeN[2]+m.edgeN[3] {
			panic("edge count mismatch") // Should never happen
		}
		m.edge = c.edge
		epv <- float64(c.total) / float64(m.count)
	}()

	// Visit each voxel and create an edge to each of its neighbors, subject to
	// edge normalization and simplification rules. All neighbors must be
	// evaluated to correctly detect orphans and calculate the average number of
	// edges per voxel.
	m.splitGrid(func(grid []bool, off int) {
		r := result{edge: make(map[Edge]struct{})}
		for i, v := range grid {
			if !v {
				continue
			}
			i += off
			ix, iy, iz := m.xyz(i)
			edgeCount := 0
			for dx := -1; dx <= 1; dx++ {
				x := ix + dx
				if x < 0 || m.d <= x {
					continue
				}
				x *= m.wh
				for dz := -1; dz <= 1; dz++ {
					z := iz + dz
					if z < 0 || m.h <= z {
						continue
					}
					z *= m.w
					for dy := -1; dy <= 1; dy++ {
						y := iy + dy
						j := x + z + y // x and z are premultiplied
						if y < 0 || m.w <= y || !m.grid[j] || i == j {
							continue
						}
						if edgeCount++; i < j {
							if e := NewEdge(i, dx, dy, dz); m.isSimple(e) {
								r.edge[e] = struct{}{}
							} else {
								r.removed++
							}
						}
					}
				}
			}
			if r.total += int64(edgeCount); edgeCount == 0 {
				// Orphans have a single edge back to themselves
				r.edge[NewEdge(i, 0, 0, 0)] = struct{}{}
			}
		}
		ch <- r
	})
	close(ch)
	m.epv = <-epv
}

// deltaMask defines all valid combinations of masking out dx, dz, and dy fields
// in order to find simpler edges.
var deltaMask = [...]Edge{
	^Edge(0x0000ff),
	^Edge(0x00ff00),
	^Edge(0x00ffff),
	^Edge(0xff0000),
	^Edge(0xff00ff),
	^Edge(0xffff00),
}

// isSimple returns true if edge e cannot be represented by multiple edges of
// lower dimensions. For example, three voxels in an L-shape will have two
// face-adjacent (1D) edges and one edge-adjacent (2D) edge. The latter can be
// removed in favor of the other two. This algorithm looks at all 2D and 3D
// edges, and tries to find alternatives by zeroing-out all possible
// combinations of dx, dz, and dy fields.
func (m *Model) isSimple(e Edge) bool {
	if e.Dim() <= 1 {
		return true
	}
	orphan := e & ^Edge(0xffffff)
	for _, mask := range deltaMask {
		if g := e & mask; g != e && g != orphan {
			// g may violate i < j rule, but that doesn't matter here
			if i, dx, dy, dz := g.Info(); m.grid[i+m.idx(dx, dy, dz)] {
				return false
			}
		}
	}
	return true
}

// idx converts xyz coordinates to a linear index.
func (m *Model) idx(x, y, z int) int {
	return x*m.wh + z*m.w + y
}

// xyz converts a linear index to xyz coordinates.
func (m *Model) xyz(i int) (x, y, z int) {
	x = i / m.wh
	i -= x * m.wh
	z = i / m.w
	y = i - z*m.w
	return
}

// splitGrid parallelizes grid processing by splitting it into multiple equal
// parts and invoking f on each part in a separate goroutine.
func (m *Model) splitGrid(f func(grid []bool, off int)) {
	c := runtime.NumCPU()
	if c < 2 || len(m.grid) < 1000 {
		f(m.grid, 0)
		return
	}
	n := len(m.grid) / c
	if n < 500 {
		n = 500
	}
	var wg sync.WaitGroup
	for off, grid := 0, m.grid; len(grid) > 0; off += n {
		if n+n > len(grid) {
			n = len(grid)
		}
		wg.Add(1)
		go func(g []bool, off int) {
			defer wg.Done()
			f(g, off)
		}(grid[:n], off)
		grid = grid[n:]
	}
	wg.Wait()
}

// pct returns a/b as a percentage.
func pct(a, b int64) float64 {
	return float64(a) / float64(b) * 100.0
}

// maxUint16 returns the larger of a or b.
func maxUint16(a, b uint16) uint16 {
	if a > b {
		return a
	}
	return b
}
