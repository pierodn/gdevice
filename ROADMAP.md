# Gdevice — Project Analysis & Improvement Plan

> Saved for resumption. Derived from codebase review (excluding README marketing copy).
> Last updated: 2026-06-13

## What This Project Is

**Gdevice** is a **header-only C++ real-time 3D toolkit** plus a **terrain-focused demo** (`walker`). It is not a full game engine; it is a modular rendering/terrain stack you can compose into your own engine.

### Solution layout

| Project | Role |
|---------|------|
| **bin2c** | Build-time tool: embeds shaders/assets as C char arrays |
| **gdevice** | Header-only library (`include/`) |
| **walker** | Fly-through demo: infinite procedural terrain |

### Core subsystems (`include/`)

| Subsystem | Purpose |
|-----------|---------|
| `gl/` | OpenGL wrappers (VBO, textures, programs, embedded shaders) |
| `type/` | GLSL-like math, OOP scene graph, terrain clipmaps, ECS stub |
| `os/` | Win32 window, input, timer, threads |
| `io/` | Image loading |
| `ut/` | Diagnostics |

### Technical identity (from code)

| Area | What the code does |
|------|-------------------|
| **Terrain LOD** | Nested clipmap-style grids (`Clipmap`, `Tile`); toroidal scroll as camera moves; multiple LOD rings (`CLIPMAPS_COUNT`) |
| **Procedural terrain** | GPU compute generates height quartets, gradients, colors, mixmaps (`generate_terrain.glsl`) |
| **Rendering** | Tessellated terrain + sky; toggles for wireframe, diffuse/specular, Fresnel, indirect, scattering, PBR, shadows (`controls.h`) |
| **Math / CPU–GPU parity** | GLSL-like types in C++ (`vec3`, `mat4`, `mix`, `clamp`) |
| **Scene management** | OOP mixin scene graph today; ECS started but incomplete (`entity.h`) |
| **Platform** | Windows + OpenGL 4.x + Win32 + MSBuild (Win32 only) |
| **License** | GPL v3 |

The `walker` demo: first-person camera over infinitely streamed, procedurally generated terrain with experimental lighting and atmosphere.

---

## Comparison With Similar Projects

### Positioning

- **Not competing with** Godot, Unity, Unreal (no editor, scripting, physics, networking, asset pipeline, cross-platform).
- **Closest peers:** GeoClipmap / CDLOD terrain engines, Outerra-style planetary terrain, research renderers with procedural worlds.
- **Partial overlap:** bgfx, Magnum, Filament, Ogre3D (rendering middleware — gdevice is more opinionated and terrain-centric).

### Strengths vs. competitors

- Very low dependency surface (OpenGL + Win32 + embedded shaders via `bin2c`)
- Coherent CPU/GPU math; good for research and teaching
- Non-trivial terrain + lighting in a compact codebase
- Rich debug toggles (F1–F12+) for fast rendering experiments

### Weaknesses vs. competitors

- Windows-only, legacy MSBuild/Win32 toolchain
- Synchronous tile GPU generation on main path (`Tile::Update()` → `glDispatchCompute`) — frame hitches during LOD scrolls
- OOP mixin architecture with many TODOs; ECS barely started
- Hardcoded materials/uniforms in `Tile::Render()`
- GPL v3 limits commercial embedding vs. MIT/Apache competitors
- No automated tests beyond `walker` integration demo
- Committed binaries in `/binaries` (see `BUILD.md`)

---

## Improvement Plan (Priority Order)

Execute in this order: **performance blockers → architecture → developer experience → differentiation**.

### P0 — Performance and stability blockers

| # | Initiative | Why first | Success metric |
|---|------------|-----------|----------------|
| **1** | **Async terrain tile streaming** | Tile invalidation triggers synchronous compute on render thread. Internal notes cite 60→20 FPS drops during LOD transitions. | Stable ~60 FPS while moving through LOD transitions |
| **2** | **Tile state machine + back-pressure** | Prerequisite for async: `pending → generating → ready → evicted` | No hitch; bounded GPU memory |
| **3** | **Memory ownership cleanup** | Raw `new`/`delete` on VBOs; incomplete destructors in `node.h`; unused `Impostor` stub | No leaks in long `walker` sessions |

**Key files:** `include/type/scene/oop/terrain/tile.h`, `clipmap.h`, `heightmap.h`, `include/type/scene/oop/node.h`, `include/os/thread.h`, `include/os/mutex.h`

---

### P1 — Architecture for extensibility

| # | Initiative | Why | Notes |
|---|------------|-----|-------|
| **4** | **ECS re-platforming** | OOP mixins fight scalability; ECS stub exists with TODOs | Migrate `Heightmap → Clipmap → Tile` to entities + systems (`LODSystem`, `TerrainGenSystem`, `TerrainRenderSystem`); preserve `walker` visual parity |
| **5** | **Data-driven material system** | `Tile::Render()` hardcodes uniforms and debug literals | `Material` struct: program, texture slots, uniform block; F-key modes become material/shader variants |
| **6** | **Impostor / far-LOD optimization** | `Impostor` declared in `node.h` but not wired | GPU time reduction at far clip distances |

**Key files:** `include/type/scene/ecs/entity.h`, `include/type/scene/oop/node.h`, `include/type/scene/oop/terrain/tile.h`, `tests/walker/source/controls.h`

---

### P2 — Developer experience and credibility

| # | Initiative | Why | Notes |
|---|------------|-----|-------|
| **7** | **CI + repo hygiene** | Committed binaries; artifact upload partially in CI | Green CI on `trunk`; release artifacts from CI, not VCS |
| **8** | **Unit tests for terrain core** | `walker` cannot guard clipmap scroll math or invalidation logic | Test `scrollTiles`, LOD blending math, transforms without GPU |
| **9** | **CMake + cross-platform shell** | MSBuild/Win32 locks out contributors | CMake on Windows first; then abstract `os/` (GLFW/SDL) |

**Key files:** `.github/workflows/msbuild.yml`, `BUILD.md`, `gdevice.sln`

---

### P3 — Differentiation (after foundation)

| # | Initiative | Rationale |
|---|------------|-----------|
| **10** | **Configurable terrain pipeline** | Move `TILE_RESOLUTION`, `CLIPMAPS_COUNT`, noise params out of headers/shaders into runtime config |
| **11** | **Static + procedural hybrid** | Import heightmaps via `io/image.h`; broaden beyond pure noise |
| **12** | **Profiler overlay** | Validate PBR/GI/scattering and async streaming costs |
| **13** | **Minimal scene editor / parameter panel** | Only after ECS + materials stabilize |
| **14** | **License strategy review** | LGPL or dual-licensing if middleware adoption is a goal |

---

## Suggested execution sequence

```
P0  Async streaming + state machine + memory cleanup     (start here)
P1  ECS → materials → impostor LOD
P2  CI cleanup → unit tests → CMake (Windows first)
P3  Config pipeline, heightmap import, profiler, editor, license
```

---

## Resume checklist (tomorrow)

When picking this up, choose one track:

- [ ] **Track A — Performance:** Start P0 #1 (async tile streaming). Read `tile.h`, `clipmap.h`, `heightmap.h`, `os/thread.h`.
- [ ] **Track B — Architecture:** Start P1 #4 (ECS). Read `entity.h`, `ECS_test.h`, `node.h`, map current mixin graph.
- [ ] **Track C — DX:** Start P2 #7 (remove committed binaries, fix CI artifacts). Read `BUILD.md`, `.github/workflows/msbuild.yml`.

Recommended default: **Track A (P0 #1)** — highest user-visible impact.

---

## Related internal docs

- `.github/copilot-instructions.md` — build commands, conventions, extended implementation roadmap
- `BUILD.md` — Windows build instructions and known repo issues
