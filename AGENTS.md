# AGENTS.md - AI Agent Guidelines for gdevice

## Project Overview

Gdevice is a lightweight, open-source ecosystem of modular and reusable components for real-time C++ applications. It aims to occupy a middle ground: providing a flexible foundation for building engine stacks without starting from zero, while serving as a platform for experimental 3D rendering. The project demonstrates that a production-ready engine can also be a foundation for other engines built on top.

**Target audience**: Researchers, indie developers, AAA studios, and open-source contributors.

## Critical Constraints - NEVER VIOLATE

### Language & Version
- **C++03 ONLY** - This is a hard constraint for readability and backward compatibility
- **No modern C++ features** (C++11, C++14, C++17, C++20, etc.) unless absolutely necessary
- **Only STL from C++03** is permitted as a 3rd party library
- **No other 3rd party libraries** - write from scratch where possible
- **Minimal C++ syntax subset** - prefer simplest constructs that achieve the goal

### Platform & Build
- **Windows ONLY** (currently) - XP, Vista, 8.1, 10, 11 supported
- **Visual Studio 2008 + MSBuild** - the ONLY supported build system
- **bin2c project must be built first** - it's a dependency for shader embedding
- **OpenGL 4.3 required** at runtime

### Code Structure
- **Header-only library** - all implementation in headers, no .cpp files in include/
- **OOP scene graph is CURRENT** - CRTP-based, static binding
- **ECS scene graph is FUTURE** - draft exists in type/scene/ecs/ but incomplete and unused
- **AI MUST use OOP approach** until ECS migration is complete and merged

## Code Style - MATCH EXACTLY

### Formatting
- **Indentation**: Tabs only (no spaces)
- **Line length**: Maximum 100 characters
- **Braces**: Always aligned (Allman style)
  ```cpp
  if (condition)
  {
      // code
  }
  ```

### Naming Conventions
| Entity Type | Convention | Example |
|-------------|------------|---------|
| Types, Classes, Functions, Constants | PascalCase | `class Heightmap`, `void Render()`, `const int MaxSize` |
| Variables, Member Variables, Parameters | pascalCase | `int tileResolution`, `vec3 position` |
| Collections/Vectors | PascalCase or pascalCase + "s" | `tiles`, `clipmaps` |
| GLSL-like types | Match GLSL naming | `vec2`, `vec3`, `mat4` |

### Comments
- **Minimal comments** - code should "speak for itself"
- **Documentation comments** - YES, use them for public APIs
- **No `// TODO:` or `// FIXME:`** - use proper issue tracking instead

## Memory & Error Handling

### Memory Management
- **RAII preferred** - use constructors/destructors for resource management
- **Avoid raw pointers and new/delete** where possible
- **Ownership semantics** - be explicit about who owns what

### Error Handling
- **NO exceptions** - disabled for maximum performance (demoscene-style)
- **Assertions** - use as pre-conditions to validate parameters
- **Logging** - for critical and non-critical states (system still evolving)
- **Debugging aids** - use existing macros in `ut/` where appropriate

## Development Workflow

### For AI Contributions
- **Treat AI like a human contributor**
- **Create branches** for changes
- **Commit changes** with descriptive messages
- **Open Pull Requests** for review
- **Do NOT push directly to trunk/main**

### Build Process
1. Build `bin2c` project first
2. Build `gdevice` project
3. Build `walker` project (depends on bin2c output)

### Testing
- **Current**: Only `walker.exe` manual testing
- **No unit/integration/performance tests** exist yet (planned future work)
- **AI must verify via tests** - no automated CI/CD currently
- **Visual verification**: Running walker.exe to spot visual regressions is possible


## Module Organization

| Directory | Purpose |
|-----------|---------|
| `include/` | Exported headers (header-only library) |
| `include/gl/` | OpenGL abstraction layer |
| `include/type/` | Core data types and algorithms |
| `include/type/scene/oop/` | **CURRENT** OOP scene graph (CRTP-based) |
| `include/type/scene/ecs/` | **FUTURE** ECS scene graph (draft, unused) |
| `include/os/` | Platform abstraction (Win32) |
| `include/io/` | File I/O abstractions |
| `include/ut/` | Utility & diagnostics |
| `tests/walker/` | Test/demo application |
| `tests/walker/source/shaders/` | GLSL shader sources |
| `binaries/` | Build output (git-ignored) |

## Key Subsystems

### Terrain LOD System
- **Hierarchical Clipmap** structure
- **Tiles** are loaded/disposed dynamically (no loading screens)
- **GPU procedural generation** - compute shaders generate height, gradient, color, mixmap
- **OpenGL 4.6** features fully embraced
- **Optimal LOD blending** - smooth transitions, no popping, no gaps
- **Tessellation** via OpenGL shaders
- **Limitation**: Heightmap-only (no overhangs/vertical detail currently)

### Shader Embedding
- **bin2c tool** converts `.glsl` -> C header -> char array
- **Shaders embedded in executable** - zero-dependency deployment
- **Workflow**:
  1. Create/edit `.glsl` in `tests/walker/source/shaders/`
  2. Build walker (triggers bin2c automatically)
  3. Include generated `.glsl.h` in C++ source

### Scene Graph (Current: OOP)
- **CRTP-based** with static binding for performance
- **Mixin-based** node system: Transform, Parent, Child, Geometry, Renderable, Updatable
- **Hierarchy**: Scene -> Heightmap -> Clipmap -> Tile

## What AI Should Avoid

- ❌ Modern C++ features (auto, lambda, smart pointers, etc.)
- ❌ 3rd party libraries (except C++03 STL)
- ❌ Exceptions
- ❌ Dynamic memory allocation (new/delete) when RAII is possible
- ❌ Modifying `type/scene/ecs/` (unused, incomplete)
- ❌ Changing the OOP scene graph without discussion
- ❌ Adding dependencies or build system changes
- ❌ Cross-platform code (Windows-only for now)

## What AI Can Do

- ✅ Write header-only C++03 code
- ✅ Use existing patterns and conventions
- ✅ Add to OOP scene graph
- ✅ Extend terrain/rendering systems
- ✅ Create new shaders (with bin2c workflow)
- ✅ Improve documentation
- ✅ Open PRs for review

## Open Questions / Future Work

These are documented for awareness but should not be implemented without discussion:
- Unit/integration/performance tests
- CMake build system port
- Cross-platform support
- ECS scene graph migration
- Unified Thread class (CPU/GPU)
- Volumetric shadows
- Renderer realism improvements (PBR, etc.)

## How to Get Help

- **For code questions**: Ask for clarification on specific subsystems
- **For design decisions**: Reference existing ADRs or create new ones
- **For build issues**: Verify VS2008 + MSBuild environment
