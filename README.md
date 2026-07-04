# Gdevice

The 3D engine landscape is broadly divided into two categories:
- Proprietary in-house engines, which power successful games but are used exclusively within the studios that develop them.
- Proprietary general-purpose engines, which showcase impressive technology and are available to all developers, but often face challenges in helping the developer achieving full production maturity across diverse use cases.

GDevice aims to occupy a middle ground by providing a lightweight, open-source ecosystem of modular and reusable components for real-time C++ applications, designed to unify CPU and GPU execution, and become a flexible foundation for building an engine stack without starting from zero. The project includes an example 3D engine for large open worlds with seamless real-time streaming and procedurally generated content at run-time.

## Core Technical Pillars

### 1. Terrain LOD system
The engine implements a sophisticated terrain LOD (Level of Detail) system based on nested grids, conceptually similar to Geometry Clipmaps but independently designed.
It enables **Infinite terrain streaming** and **smooth level transitions** in its GPU-based implementation.

### 2. Procedural Generation
The terrain is dynamically generated on the GPU, computing heights, gradients and textures in real time.

### 3. GLSL-like Math Library for C++
The engine provides `vec2`, `vec3`, `mat4`, and related types with syntax closely matching GLSL, including operator overloading and common functions like `mix`, `clamp`, and `fract`. This design allows math logic to be easily shared or ported between CPU and GPU.

### 4. Scene Management
The current `walker` demo primarily uses an object-oriented approach to scene management, while a transition to an ECS (Entity Component System) architecture is planned for a more data-oriented design.

### 5. Rendering
A **quasi physically-based rendering (PBR)** pipeline is implemented, along with **global illumination** and **atmospheric scattering**. 

### 6. Platform & Performance Focus
- Highly optimized for Windows, leveraging SSE (Streaming SIMD Extensions) for math operations and Multi-threading.
- It uses specialized `VertexBuffer` and `IndexBuffer` abstractions, which also support textures-based data storage (e.g. storing vertex attributes such as gradients and blend maps in textures).
- A custom tool compiles GLSL shaders into C headers, embedding the shader source directly into the executable for zero-dependency deployment.

## Project Details

### Dependencies:
- C++
- OpenGL 4.5
- Windows

### Motivation
This project revisits and simplifies many best practices commonly used in professional environments, aiming to explore fundamental design questions such as:

- Can a 3D engine solution **build in seconds**—say, five?
- What is the **simplest design** that still works?
- Can C++ source code be written in a **DRY** style, similar to Java?
- How can a **continuous, unsandboxed terrain LOD** scheme possibly operate?
- How realistic can **quasi physically-based rendering** become?
- Is some form of **global illumination** achievable without pre-baking?
- How could **atmospheric scattering** and general **volumetric rendering** be implemented efficiently?
- What roles can **procedurally generated content** serve?
- Can **materials** (albedo, normals, etc.) be **generated dynamically at runtime**?

### Sponsor
What began as a proof of concept for a terrain level-of-detail algorithm and a minimalist 3D engine gradually evolved into an open-ended journey through the uncharted intersections of real-time rendering and parallel computing, beyond the boundaries of enterprise settings where risk-taking is tightly constrained.

If you find what I do useful, or if you simply care about empowering creative people to keep building and innovating, please consider supporting a few hours of development or just offering a cup of coffee.

Thank you!

![donate](https://github.com/user-attachments/assets/d4f11a78-c1b6-4c04-b032-c415947f2d0b)

### Screenshots
![Sun](https://github.com/user-attachments/assets/f061be3f-a293-44a2-b02d-c0e4bd9fed98)
![Materials1](https://github.com/user-attachments/assets/1cab3e01-b878-43fe-9e62-9348598797ef)
![Materials2](https://github.com/user-attachments/assets/57e0a62e-6da3-4bfa-9341-4e11d44f6276)
![Volumetric](https://github.com/user-attachments/assets/f97e25a5-8c3f-4df5-ba02-64bb927b4685)
![Terrain1](https://github.com/user-attachments/assets/56676b41-7cb0-4834-988a-5256f3787158)
![Terrain2](https://github.com/user-attachments/assets/743756f9-2791-404a-8d45-71192bc4f742)
![Terrain3](https://github.com/user-attachments/assets/6a7091ea-22c6-4d0d-9846-e164404d735c)
![Sunset](https://github.com/user-attachments/assets/4c3a4106-feac-41e2-befa-4f71630435ca)

## Documentation

- **[CONTEXT.md](CONTEXT.md)** - Project vision, technical philosophy, architecture decisions, and historical context
- **[AGENTS.md](AGENTS.md)** - Guidelines for AI agents contributing to this repository
- **[BUILD.md](BUILD.md)** - Build instructions and environment setup
- **[ROADMAP.md](ROADMAP.md)** - Future plans and milestones

## Getting Started
See [BUILD.md](BUILD.md) for build instructions.
