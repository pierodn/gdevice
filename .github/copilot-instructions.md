# Copilot Instructions for gdevice

## Build & Test

### Build Commands
- **Full build (Release)**: From Developer Command Prompt or PowerShell:
  ```
  nuget restore .\gdevice.sln
  msbuild /m /p:Configuration=Release /p:Platform=Win32 .\gdevice.sln
  ```
- **Debug build**: Replace `Release` with `Debug` in the msbuild command
- **Key constraint**: Always specify `/p:Platform=Win32` — the solution only targets 32-bit builds

### Running Tests
- **walker** test application:
  1. Build walker project (listed in gdevice.sln with bin2c as a dependency)
  2. Run from repository root or with working directory `tests\walker`:
     ```
     .\binaries\release\walker.exe
     ```
  3. Walker uses assets from `tests/walker/assets/` and built shaders

### Dependencies
- **NuGet**: Run `nuget restore` before building
- **Windows SDK & MSBuild**: Required for C++ compilation; ensure appropriate version for target Visual Studio version

## Architecture

### Project Layout
The solution contains three projects:
- **bin2c**: Build-time tool that converts binary files to C header files (embeds data as char arrays)
- **gdevice**: Header-only library with core abstractions organized by subsystem
- **walker**: Test/demo application showing terrain rendering with hierarchical clipmaps and scene hierarchy

### Core Subsystems (under `include/`)
Each subsystem is self-contained and organized by concern:

- **`gl/`** — OpenGL abstraction layer
  - VBO management (VertexBuffer, IndexBuffer, Buffer)
  - Texture management, Shader programs, Links
  - glpp.h contains wrapper functions for common VBO operations (create, bind, deallocate)
  - Note: Header-only; no compiled library

- **`type/`** — Core data types and algorithms
  - Scene graph: `type/scene/oop/` (OOP-style with node hierarchy, transforms, scene management)
  - ECS system: `type/scene/ecs/` (Entity Component System framework for more data-driven design)
  - Terrain: `type/scene/oop/terrain/` (Heightmap, Clipmap, Tile with LOD support)
  - Utilities: CPU types (RGBA, int8-64, float packing), GLSL type support

- **`os/`** — Platform abstraction
  - Window management, Timer, Keyboard input, Thread, Mutex, Compiler detection

- **`io/`** — File I/O abstractions
  - Model/mesh loading, Image loading

- **`ut/`** — Utility & diagnostics
  - Logging, Diagnostics macros (DEBUG_CRITICAL, etc.)

### Key Design Patterns

**Header-Only Library**
- All implementation in headers (no .cpp files in include/)
- Enables inline optimization and template instantiation
- bin2c and walker only projects that actually compile C++ code

**Shader Embedding**
- Shaders are compiled to binary, then embedded using bin2c tool (converts to char arrays in .h files)
- .gitignore ignores `*.glsl.h` generated files
- Build process: binary shader → bin2c → char array → linked into executable

**Manual Dependency Management**
- walker explicitly depends on bin2c (see gdevice.sln ProjectDependencies)
- Build gdevice project before walker to ensure bin2c output is available

**Dual Scene Systems**
- OOP approach: Node-based hierarchy, explicit transforms (type/scene/oop/)
- ECS approach: Data-driven entity/component/system pattern (type/scene/ecs/)
- Projects can use either or both depending on use case

## Key Conventions

### Naming
- **Namespaces**: `GL::`, `GL::VBO::` for graphics abstractions
- **Type shortcuts**: Common GLSL types (vec2, vec3, vec4, dmat3, dvec3, rgba) used directly (no std:: prefix)
- **Typedefs**: 
  - Fixed-width integers: `int8` through `int64`, `uint8` through `uint64`
  - Byte type alias: `byte` (unsigned char)
  - Color: `rgba` (vec<byte,4>)

### Directory/File Organization
- `include/` — Exported headers
- `tests/walker/` — Test/demo executable and its assets
- `tests/walker/assets/` — Shader sources and test data
- `binaries/` — Build output (excluded from VCS, generated at build time)

### Platform Specifics
- Windows-only codebase (Win32 API)
- All configuration and file paths use backslashes in command line examples
- Compiler detection via `os/compiler.h`

## Workflow Notes

- **Building bin2c first**: walker depends on bin2c; ensure bin2c is built before rebuilding walker to regenerate shader header files
- **Working directory**: walker test expects to run from repo root or `tests\walker/` directory
- **CI/CD**: GitHub Actions workflow (msbuild.yml) automates builds and uploads walker artifacts; currently archives walker executables and test directory

## Shader Compilation Workflow

### How Shaders are Embedded
Shaders follow a two-stage compilation process:

1. **Shader → Binary**: GLSL source compiled to intermediate binary (occurs during walker build)
2. **Binary → C Header**: bin2c tool converts binary to C char arrays (embedded in walker executable)

### Adding or Modifying Shaders
- **Shader sources**: Located in `tests/walker/source/shaders/` (GLSL files)
- **Generated headers**: `*.glsl.h` files (auto-generated, git-ignored per .gitignore)
- **Workflow**:
  1. Create/edit `.glsl` file in `tests/walker/source/shaders/`
  2. Build walker project (triggers bin2c tool automatically)
  3. `bin2c` converts shader binary to `filename.glsl.h`
  4. Include the generated `.glsl.h` in your C++ source
- **Configuration**: Shader generation programs are hardcoded in heightmap.h:
  ```cpp
  generator.Build(generate_terrain_glsl);   // from generate_terrain.glsl.h
  renderer.Build(render_terrain_glsl);      // from render_terrain.glsl.h
  ```
  Update these includes if adding new shaders

### Shader Structure (from walker)
Shaders are multi-stage programs defined as strings:
- **Format**: Each shader type prefixed in source (e.g., `VERTEX_SHADER:\n#version 430\n...`)
- **Parsing**: `Program::BuildShaders()` splits multi-type shader source and compiles separately

---

## Scene Hierarchy & Transform System

### Scene Graph Structure
Scenes are built from composable node types connected via inheritance mixins:

```
Scene (Node + Transform + Parent + SceneState)
  ├─ Heightmap (Node + Parent + Child + Transform)
  │   ├─ Clipmap (Node + Transform + Parent + Child)  // LOD level 0
  │   │   ├─ Tile (Node + Transform + Child + Geometry + Updatable + Renderable)
  │   │   ├─ Tile
  │   │   └─ ... (CLIPMAP_SIZE × CLIPMAP_SIZE tiles per LOD)
  │   │
  │   ├─ Clipmap (LOD level 1)
  │   │   └─ ... (coarser tiles, same layout)
  │   └─ ... (up to CLIPMAPS_COUNT LODs)
```

### Node Types (Mixins)
Defined in `type/scene/oop/node.h`:
- **Node**: Base class, virtual destructor
- **Transform**: Holds `mat3 transform` (position, rotation, scale)
- **Parent**: Owns vector of `Child*` pointers
- **Child**: Back-pointer to parent
- **Geometry**: Stores vertex/index buffers (VBO, IBO)
- **Renderable**: Virtual `Render(NodeState&, SceneState&)` method
- **Updatable**: Virtual `Update()` method; marks videomem_invalidated, hostmem_invalidated flags
- **Impostor**: Pre-cached impostor texture + direction/distance (TODO: not fully implemented)

### Transform System
- **Per-node transforms**: Each node has position, rotation, scale in `Transform::transform` (mat3 with vec3 fields)
- **Stack-based hierarchy**: Scene traversal maintains `ModelViewMatrixStack` and accumulates `ModelViewMatrix`
- **Matrix composition**: Child transforms composed via `nodeState.ModelViewMatrix *= TransformationMatrix(node.transform)`

### Traversal & Rendering
Scene traversal happens in `Scene::Traverse()`:
1. Push current ModelViewMatrix to stack if node has Transform
2. Compose new ModelViewMatrix from parent + local transform
3. Check Updatable: if `videomem_invalidated`, call `Update()` and regenerate GPU data
4. Check Renderable: if not skipping render, call `Render(nodeState, sceneState)`
5. Recurse into children (Parent nodes)

---

## Terrain & Clipmap System

### Overview
Terrain is generated procedurally using a hierarchical **Clipmap** structure for infinite terrain:
- **Dynamic LODs**: Multiple resolution levels centered on camera
- **Tiles**: Each clipmap divided into tiles (default 33×33 vertices per tile)
- **Compute Shaders**: Tile data (height, gradient, color) generated via GPU compute dispatch

### Clipmap Parameters (in `parameters.h`)
```cpp
TILE_RESOLUTION = 33         // Must be 2^N+1 (e.g., 17, 33, 65, 129)
CLIPMAPS_COUNT = 5           // Number of LOD levels
CLIPMAP_WINDOW = 10          // Width of active tile region per clipmap
CLIPMAP_SIZE = 18            // Total grid size (includes padding and buffering)
```

### Clipmap Structure
- **LOD levels**: Each clipmap represents one resolution level (LOD 0 = finest, LOD N = coarsest)
- **Tile grid**: Each clipmap is CLIPMAP_SIZE × CLIPMAP_SIZE tiles
- **Tiling size**: Tile size = 2^LOD (so LOD 0 tiles are smallest, LOD 2 tiles are 4× larger)
- **Dynamic centering**: `Clipmap::move(location)` repositions tile grid around camera, marking invalidated tiles

### Tile System
- **Tile data**: Each tile stores 4 textures (in VertexBuffer):
  - `quartets`: Height/gradient data (float4 per vertex)
  - `gradients`: Slope/derivative data (float4)
  - `colors`: Diffuse color/material data (float4)
  - `mixmaps`: Texture blend weights (float4)
- **GPU generation**: Tile::Update() dispatches compute shader to regenerate these textures
- **CPU-side tracking**: 
  - `videomem_invalidated`: GPU data needs regeneration (set when tile moves into view)
  - `hostmem_invalidated`: CPU metadata needs sync (set after GPU update)

### LOD Distance Calculation
```cpp
GetVisibilityDistance() = 4 * (1 << children.size())
```
If you have 5 clipmaps, visibility = 4 × 2^5 = 128 (world units)

### Heightmap Movement
- `Heightmap::moveAt(camera)`: Called each frame with camera position
- Updates all clipmaps and marks tiles as invalidated based on scrolling
- Returns Z value (height at camera location)

---

## Rendering Pipeline

### Initialization (walker.cpp OnOpen)
```cpp
GL::Initialize();              // OpenGL context setup
heightmap.setTileResolution(TILE_RESOLUTION);  // Configure tile grid
heightmap.setLODs(CLIPMAPS_COUNT);             // Create LOD clipmaps
heightmap.Initialize();        // Build compute/render shaders
scene.children.push_back(&heightmap);
scene.Initialize();            // Build sky dome shader
```

### Per-Frame Render Loop
1. **Update camera & scene state**
   - Update camera position/rotation from keyboard input
   - Call `heightmap.moveAt(camera)` to update tile positions
   
2. **GPU tile generation** (in Scene::Traverse)
   - For each invalidated tile, call `Tile::Update()`
   - Dispatch compute shader: `glDispatchCompute(width/2+1, height/2+1, 1)`
   - Outputs: height, gradient, color, mixmap textures

3. **Scene traversal** (Scene::Traverse)
   - Recursively traverse scene graph (heightmap → all clipmaps → all tiles)
   - Maintain transform matrix stack
   - For each renderable tile, call `Tile::Render(nodeState, sceneState)`

4. **Tile rendering** (Tile::Render)
   - Bind tile's VBO (vertex buffer with position data)
   - Bind tile's texture units (quartets, gradients, colors, mixmaps)
   - Set uniforms from Controls (rendering flags: wireframe, debug mode, lighting options)
   - Draw indexed: `glDrawElements()` using shared IBO across all tiles
   - Include scene state (sun position, camera matrix)

### Matrix Transformations
- **Projection**: Set once per frame in nodeState
- **View**: Camera matrix (inverse of camera transform)
- **Model**: Accumulated from scene hierarchy (heightmap → clipmap → tile)
- **Final vertex**: `gl_Position = Projection * View * Model * vertex`

### Render State Management
- **Shared IBO**: All tiles use same IndexBuffer (pre-computed grid topology)
- **Per-tile VBO**: Each tile has unique VertexBuffer with position data
- **Texture units**: Quartets (0), Gradients (1), Colors (2), Mixmaps (3)
- **Program switching**: Render shader set per tile (shared across all tiles in practice)

---

## Testing Patterns & Controls

### Walker Test Application
Walker is an interactive terrain viewer demonstrating the complete rendering pipeline:
- **Run**: `walker.exe` (from `binaries/release/` or `binaries/debug/`)
- **Window**: 800×600 by default (configurable in walker.cpp)

### Interactive Controls (F1-F12, H, G, C, U, T, V)
Defined in `controls.h`:
- **F1**: Wireframe toggle
- **F2**: Debug mode (Color, HeightBlend, Normal, Light)
- **F3**: Diffuse lighting
- **F4**: Specular lighting
- **F5**: Fresnel shading
- **F6**: Sky lighting
- **F7**: Indirect lighting
- **F8**: Light scattering
- **F9**: Tessellation
- **F10**: Bump mapping
- **F11**: Shadow mapping
- **F12**: PBR mode
- **H**: Heatmap
- **G**: Gamma correction
- **C**: Contrast adjustment
- **U**: Color saturation
- **T**: Color tint
- **V**: Vignetting

### Keyboard Navigation (Camera Control)
- **WASD**: Move camera forward/left/back/right
- **Mouse/arrows**: Rotate view
- **Shift**: Speed boost (multiplied by BOOST_FACTOR = 10)
- **Page Up/Down**: Altitude control

### Adding New Test Scenarios
1. Create new test class inheriting from `Program` and `Listener` (like Walker)
2. Implement `OnOpen()`: Initialize scene and shaders
3. Implement render loop: Update transforms, traverse scene, render
4. Add to gdevice.sln and configure as test project

---

## Asset Management

### Directory Structure
```
tests/walker/
├── assets/              # Textures, models, procedural data
├── source/
│   ├── walker.cpp       # Main test executable
│   ├── parameters.h     # Terrain/rendering configuration
│   ├── controls.h       # Input key bindings
│   ├── shaders/
│   │   ├── generate_terrain.glsl   # Compute shader for height/gradient generation
│   │   ├── render_terrain.glsl     # Fragment/vertex shader for rendering
│   │   ├── render_sky.glsl         # Sky dome shader
│   │   └── (auto-generates .glsl.h files)
│   ├── textures.h       # Texture format/binding constants
│   ├── ui/
│   └── (test-specific code)
```

### Adding Textures
- Place image files in `tests/walker/assets/`
- Use `io/image.h` to load: `Image::Load("path/to/texture.png")`
- Bind to texture units via `GL::Texturing::bind(unit, texture)`
- Note: Current code uses procedural textures; file I/O partially integrated

### Adding Models/Meshes
- Place model files (OBJ, etc.) in `tests/walker/assets/`
- Use `io/model.h` to load: `Model::Load("path/to/model.obj")`
- Create Geometry node with loaded VBO/IBO
- Attach to scene as child node

### Shader Assets
- GLSL source in `tests/walker/source/shaders/`
- Must be valid for target GPU (walker assumes OpenGL 4.3+)
- Reference in heightmap.h or create new Program object
- Build system automatically runs bin2c on shader binaries

---

## Memory & Performance Considerations

### GPU Memory
- **Per-tile VBO**: Contains 4 textures (quartets, gradients, colors, mixmaps)
  - Size: `(TILE_RESOLUTION * TILE_RESOLUTION * 4 floats * 4 textures) * 4 bytes`
  - Example: 33×33 tile = 17,424 floats per texture = ~272 KB per texture, ~1 MB per tile
- **Shared IBO**: Single IndexBuffer used by all tiles (size: TILE_RESOLUTION²)
- **Total VRAM** (5 LODs, 18×18 tile grid): ~100+ MB

### CPU Memory
- **Scene hierarchy**: O(clipmaps × tiles per clipmap) nodes
  - 5 LODs × 18×18 tiles = 1,620 Tile objects (lightweight, just pointers + IDs)
- **Transform stacks**: Per-frame allocation in NodeState (reserve 4096 matrices)

### Optimization Opportunities (TODO in code)
- **Material system**: Currently hardcoded; should be data-driven Material* per node
- **VAO caching**: Vertex attribute state should be cached in VAO objects
- **Impostor system**: Pre-render terrain to impostor texture when far from camera (reduces geometry)
- **Fixed-point terrain**: Code has commented-out fixed-point implementation for CPU terrain generation
- **Compute shader optimization**: Tile generation dispatch could batch multiple tiles

### Invalidation Strategy
- **Lazy evaluation**: Tiles only regenerate when `videomem_invalidated = true`
- **Hierarchical updates**: Parent clipmaps update child tiles on camera movement
- **Streaming**: Coarser LODs update less frequently than fine LODs

### Threading (Partial)
- **Mutex support**: `os/mutex.h` and `os/thread.h` available
- **Current usage**: Minimal; mostly synchronous render thread
- **Future**: Could parallelize tile generation across threads or stream LODs asynchronously

---

## Common Issues

- **Linking fails**: Verify Windows SDK and C++ build tools match Visual Studio version
- **Walker crashes**: Ensure you're running from correct working directory (repo root or tests\walker) so asset paths resolve
- **NuGet restore fails**: Let Visual Studio open the solution and handle restoration automatically, or ensure NuGet CLI is in PATH
- **Shaders don't update**: After editing `.glsl` files, rebuild walker to regenerate `.glsl.h` headers (clean and rebuild if caching issues)
- **Tiles not visible**: Check camera position vs. visibility distance; verify LOD clipmaps initialized via setLODs()
- **Compute shader errors**: Ensure target GPU supports OpenGL 4.3+ and image texture units (glBindImageTexture)

---

## Implementation Roadmap

This roadmap aligns gdevice with modern game engine architecture while maintaining C++03 compatibility and Windows foundation stability. C++17 migration is deferred and not a blocker.

### Phase 1: Architecture Refactor (Windows Foundation)

**1. ECS Scene Graph Re-platforming** (XLarge effort)
   - **Goal**: Replace OOP node mixins (Transform, Parent, Child, Geometry, Renderable, Updatable) with pure ECS architecture
   - **Scope**: 
     - Remove mixin-based Node system from `type/scene/oop/node.h`
     - Implement ECS framework: entities (IDs), components (data), systems (logic)
     - Migrate Heightmap/Clipmap/Tile hierarchy to ECS entities with domain-specific systems (TileUpdate, TerrainRender, LODManagement)
     - Preserve parent-child hierarchy as component data (Parent/Child components)
   - **Rationale**: ECS enables flexibility, data-driven design, better cache locality, and parallelizable systems (aligns with UE, Godot, Unity patterns)
   - **Testing**: Verify walker renders identically after refactor; all 18 debug modes (F1-F12) functional

**2. Async Terrain Tile Streaming** (XLarge effort)
   - **Goal**: Move tile GPU generation to worker threads; prevent frame stalls during LOD transitions
   - **Scope**:
     - Implement work queue using existing `os/thread.h` + `os/mutex.h` (C++03 compatible)
     - Tile state machine: pending → ready → rendered
     - Worker threads dispatch compute shaders, handle GPU readback
     - Main render thread queues work, polls for ready tiles, renders without blocking
   - **Rationale**: Current synchronous generation causes frame drops (60→20 FPS); competitors use async LOD streaming
   - **Success metric**: Maintain 60 FPS during camera movement through all LOD transitions

**3. Data-Driven Material System** (Large effort)
   - **Goal**: Replace hardcoded shader uniforms with flexible material definitions
   - **Scope**:
     - Define Material struct: shader program, texture bindings, uniform values
     - Remove hardcoded uniforms from `Tile::Render()` (lighting, diffuse, specular, fresnel, etc.)
     - Allow per-tile material assignment (enables varied terrain types, blending)
     - Keep 18 debug modes (F1-F12) as shader variations/material overrides
   - **Rationale**: Current hardcoded approach is inflexible; modern engines use data-driven materials
   - **Dependency**: None; ECS refactor not required but will simplify component design

**4. Impostor LOD System** (Medium effort)
   - **Goal**: Pre-render distant terrain to impostor textures; reduce geometry at far LODs
   - **Scope**:
     - Render high-detail terrain to off-screen texture at coarse LOD positions
     - Use impostor quad instead of full tile geometry when far from camera
     - Implement impostor distance calculation per clipmap
   - **Rationale**: Reduces draw calls and vertex processing; common optimization in terrain engines
   - **Success metric**: 30%+ frame time reduction on low-end GPUs

### Phase 2: Portability & Infrastructure (Later)

**5. CMake Build System + Cross-Platform** (Medium + Large efforts, bundled)
   - **Goal**: Replace MSBuild .vcproj with CMake; enable Linux/macOS builds
   - **Scope**: 
     - Migrate bin2c, gdevice, walker to CMakeLists.txt
     - Automate shader .glsl → .glsl.h generation via custom_command
     - Abstract Win32 APIs (window, threading, timer) → portable libraries (GLFW, pthreads, std::chrono equivalent)
   - **Rationale**: MSBuild is Windows-only; CMake enables cross-platform CI/CD; game engines use CMake or equivalent
   - **Deferred reason**: Requires platform abstraction layer refactor; better after Windows foundation is solid
   - **Timeline**: Start after Phase 1 complete

**6. Unit Testing Framework** (Medium effort, deferred)
   - **Goal**: Add unit tests for terrain generation, transforms, rendering paths
   - **Scope**: Integrate catch2/gtest; test terrain LOD calculations, tile invalidation logic, matrix transforms
   - **Rationale**: Walker is integration test only; unit tests improve quality parity with industry
   - **Deferred reason**: Core functionality validation more critical than test suite initially

### Deferred (Low Priority)

- **Editor tool** (XLarge): Scene tree visualization, parameter tweaking — build after core systems stable
- **Profiler** (Medium): GPU/CPU profiling overlay
- **Shader graph editor** (Large): Visual shader authoring
- **CI/CD automation** (Small): GitHub Actions for builds (tied to CMake phase)
- **Voxel support** (XLarge): Alongside heightmap
- **Networking** (XLarge): Multiplayer terrain sync
- **Spatial audio** (Large): LOD-based audio system
- **C++17 migration** (Not needed): Deferred indefinitely; C++03 + RAII sufficient for all improvements

### Architecture Decision: RAII Without Smart Pointers

Rather than retrofitting smart pointers (unique_ptr, shared_ptr), we fix memory management by redesigning ownership patterns:
- Explicit destructors with proper cleanup (Clipmap, Tile, Heightmap)
- Stack allocation where possible (Scene graph nodes)
- Careful move semantics using C++03 patterns (swap idiom, explicit ownership transfer)
- This keeps the codebase clean, C++03 compatible, and eliminates smart pointer overhead
