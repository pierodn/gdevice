# CONTEXT.md - Gdevice Project Context

## Vision & Goals

Gdevice aims to occupy a unique position in the 3D engine landscape:

**The Problem**: The 3D engine world is divided between:
- Proprietary in-house engines (powerful but exclusive to their studios)
- Proprietary general-purpose engines (impressive but may struggle with helping developers to reach production maturity across diverse use cases)

**The Solution**: Gdevice provides a **middle ground** - a lightweight, open-source ecosystem of modular, reusable components that:
- Gives developers a **flexible foundation** for building their own engine stack
- Avoids the "start from zero" problem
- Can serve as a **production-ready engine** itself

### Target Audience
- **Researchers** - experimenting with rendering techniques, terrain algorithms
- **Indie developers** - building games without reinventing core systems
- **AAA studios** - prototyping or integrating specific subsystems
- **Open-source community** - collaborating on a shared foundation

### Core Philosophy
- **Simplicity over complexity** - Can a 3D engine solution build in seconds?
- **Minimalism** - What is the simplest design that still works?
- **DRY principle** - Can C++ source code be written in a DRY style?
- **Performance focus** - Highly optimized for Windows with SSE and multi-threading
- **Zero-dependency deployment** - Shaders and assets embedded in executables

## Technical Philosophy

### Why C++03?
The choice of C++03 is deliberate and multi-faceted:

1. **Readability**: A restricted subset of C++ forces cleaner, more understandable code
2. **Backward Compatibility**: Supports older compilers and platforms (Visual Studio 2008, Windows XP+)
3. **Minimal Dependencies**: Only STL from C++03 is used; everything else is written from scratch
4. **Performance**: Avoids runtime overhead from modern features
5. **Experimental**: Demonstrates that production-quality code can be written with minimal language features

### Why Header-Only?
The header-only library design was chosen to:
- Enable **faster iteration** (no separate compilation step)
- Allow **easy integration** into other projects (just include headers)
- Support **template-heavy code** naturally
- Avoid **binary compatibility** issues
- Test the viability of this approach at scale (result: no problems found, codebase remains small)

### Why Windows-Only?
- **Focus**: Allows deep optimization for a single platform
- **Simplicity**: No cross-platform abstraction overhead
- **Target audience**: Many game developers still target Windows primarily
- **Future**: Cross-platform is a **planned goal**, not a current requirement

### Why OpenGL?
- **Mature**: Well-established, widely supported
- **Explicit**: More control over GPU operations compared to higher-level APIs
- **Compatibility**: Works across different GPU vendors
- **Version**: OpenGL 4.5/4.6 provides compute shaders and modern features needed for the terrain system

## Architecture Decisions

### Scene Graph: OOP with CRTP
**Decision**: Use OOP-based scene graph with CRTP (Curiously Recurring Template Pattern) for static binding.

**Rationale**:
- First version used a single node type with static binding for high performance
- CRTP enables polymorphic behavior without virtual function overhead
- Concerns about flexibility led to exploring ECS as an alternative

**Status**: **CURRENT**, used throughout the codebase

**Future**: Migration to ECS is planned (draft exists in `type/scene/ecs/` but is incomplete and unused)

### Terrain System: Hierarchical Clipmaps
**Decision**: Implement terrain LOD using a hierarchical Clipmap structure with GPU procedural generation.

**Rationale**:
- Enables **infinite terrain streaming** without loading screens
- **Dynamic tile management**: Only tiles entering view need loading, tiles exiting can be disposed
- **GPU acceleration**: Compute shaders generate height, gradient, color, and mixmap data in real-time
- **Smooth transitions**: Optimal LOD blending eliminates popping and gaps between tiles
- **Evolution**: Started from CPU-based OpenGL 1.3 fixed-function, evolved to OpenGL 4.6 compute shaders

**Tradeoffs**:
- **Pros**: No loading times, seamless infinite world, optimal performance
- **Cons**: Heightmap-only (cannot encode overhangs), limited to terrain representation

**Future Considerations**:
- Runtime tessellation for vertical detail (hangovers, cliffs)
- Enhanced geometric complexity while maintaining performance

### Shader Embedding: bin2c Tool
**Decision**: Embed shaders as C char arrays using the bin2c tool.

**Rationale**:
- **Zero-dependency deployment**: All assets bundled in the executable
- **Simplified distribution**: No need to ship separate shader files
- **Integration**: CPU and GPU code can be mixed together in the project
- **Workflow**: GLSL source → binary → C header → executable

**Status**: bin2c is a placeholder that could evolve into a more sophisticated shader compiler with advanced capabilities (future possibility, not guaranteed)

## Historical Context

### Origins
Gdevice began as a **proof of concept** for:
1. A terrain level-of-detail algorithm
2. A minimalist 3D engine

Over time, it evolved into an **open-ended journey** through intersections of real-time rendering and parallel computing, beyond typical enterprise constraints.

### Evolution Path
- **CPU Era**: Initial terrain generation on CPU with OpenGL 1.3 fixed-function shaders
- **GPU Era**: Migration to OpenGL 4.6 with compute shaders for procedural generation
- **Architecture Era**: Experimentation with different scene graph approaches (OOP → ECS exploration)
- **Future Era**: Planned improvements in realism, testing, portability

### Key Insights Gained
- Header-only libraries can scale to production use
- C++03 is sufficient for complex, modern and high-performance scenarios
- CRTP provides static polymorphism benefits
- GPU procedural generation enables infinite, seamless terrain
- Minimal dependencies simplify maintenance

## Technical Debt & Future Roadmap

### High Priority (Next Steps)
1. **ECS Scene Graph Migration**
   - Replace OOP mixin-based system with pure ECS
   - Benefits: Better cache locality, data-driven design, parallelizable systems
   - Scope: Remove `type/scene/oop/`, implement ECS framework, migrate Heightmap/Clipmap/Tile
   - Success: Walker renders identically, all F1-F12 debug modes functional

2. **Unit/Integration/Performance Tests**
   - Currently: Only manual testing via walker.exe
   - Goal: Comprehensive automated test coverage
   - Scope: Test terrain LOD calculations, matrix transforms, tile invalidation logic
   - Benefit: Prevent regressions, enable confident refactoring

### Medium Priority (Phase 2)
3. **CMake Build System + Cross-Platform**
   - Migrate from MSBuild .vcproj to CMakeLists.txt
   - Abstract Win32 APIs (window, threading, timer)
   - Enable Linux/macOS builds
   - Automate shader .glsl → .glsl.h generation

4. **Data-Driven Material System**
   - Replace hardcoded shader uniforms with flexible material definitions
   - Allow per-tile material assignment
   - Keep F1-F12 debug modes as shader variations

5. **Async Terrain Tile Streaming**
   - Move tile GPU generation to worker threads
   - Prevent frame stalls during LOD transitions
   - Goal: Maintain 60 FPS during camera movement

6. **Impostor LOD System**
   - Pre-render distant terrain to impostor textures
   - Reduce geometry at far LODs
   - Success metric: 30%+ frame time reduction on low-end GPUs

### Long-Term / Experimental
- **Editor tool**: Scene tree visualization, parameter tweaking
- **Profiler**: GPU/CPU profiling overlay
- **Shader graph editor**: Visual shader authoring
- **Voxel support**: Alongside heightmap
- **Networking**: Multiplayer terrain sync
- **Spatial audio**: LOD-based audio system

## Domain Context: Real-Time Rendering

### Key Challenges Addressed
- **Infinite terrain**: How can continuous, unsandboxed terrain LOD operate?
- **Procedural generation**: Can materials be generated dynamically at runtime?
- **Quasi PBR**: How realistic can physically-based rendering become without full PBR?
- **Global illumination**: Is GI achievable without pre-baking?
- **Atmospheric scattering**: How to implement efficiently?
- **Volumetric rendering**: Techniques for atmospheric and shadow effects

### Competitive Landscape
Unlike Unity, Unreal, or Godot:
- Gdevice is **not a complete engine** - it's a foundation
- **No editor** (planned future work)
- **No asset pipeline** (assets are procedural or simple)
- **No plugin system** (yet)
- **Header-only** approach vs. compiled libraries
- **C++03** vs. modern C++

### Research Value
Gdevice serves as a **testbed** for:
- Terrain LOD algorithm research
- GPU procedural generation techniques
- Static polymorphism patterns (CRTP)
- Minimal-dependency engine architecture
- OpenGL 4.x feature utilization

## Why This Matters

Gdevice represents an important exploration in engine design:

**Can simplicity scale?** Most engines grow complex over time. Gdevice asks: can we keep the core simple while still being production-ready?

**Can constraints drive creativity?** By limiting to C++03, Windows, OpenGL, and minimal dependencies, Gdevice demonstrates that sophisticated rendering is possible without modern complexities.

**Can foundations be flexible?** Gdevice is designed to be both a usable engine AND a building block for other engines - a dual purpose that requires careful architectural decisions.

**Can open development work for engines?** Unlike proprietary engines, Gdevice evolves in the open, allowing the community to learn from and contribute to its development.

---

*See ROADMAP.md for detailed milestone planning and BUILD.md for development setup instructions.*
