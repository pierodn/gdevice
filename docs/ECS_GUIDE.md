# ECS Design Guide — Better Approaches & C++ Examples

> Saved for resumption (2026-06-13). Companion to `ROADMAP.md` P1 #4 (ECS re-platforming).

## Limits of the current `entity.h` design

| Issue | Why it matters |
|--------|----------------|
| `map<Entity, vector<T>>` per type | Slow lookups, poor cache locality, heap churn |
| `GetComponents(entity)` auto-inserts | Accidental “ghost” entries via `map::operator[]` |
| `EntityDestroyerFn` registry | Workaround: no central registry of component types |
| Header-only static maps | Risk of one storage copy per translation unit |
| Entity = monotonic `uint64` | No generation counter → stale handles after destroy/reuse |

`EntityDestroyerFn` is fine for a prototype. It is not how most production ECS libraries scale.

---

## Better approaches (most → least common in games)

### 1. Sparse set ECS (EnTT-style) — recommended default

Each component type uses a **dense array** + **sparse index** by entity id. Destroy = swap-and-pop. Iteration = linear scan (cache-friendly).

```cpp
#include <vector>
#include <cstdint>
#include <cassert>

struct Entity {
    std::uint32_t index;
    std::uint32_t generation;
};

template<typename T>
class ComponentPool {
public:
    T& emplace(Entity e, T value = {}) {
        ensureSparse(e.index);
        if (sparse_[e.index] == null_) {
            sparse_[e.index] = static_cast<std::uint32_t>(dense_.size());
            dense_.push_back(std::move(value));
            entities_.push_back(e);
        } else {
            dense_[sparse_[e.index]] = std::move(value);
        }
        return dense_[sparse_[e.index]];
    }

    T* get(Entity e) {
        if (e.index >= sparse_.size() || sparse_[e.index] == null_) return nullptr;
        Entity stored = entities_[sparse_[e.index]];
        if (stored.generation != e.generation) return nullptr; // stale handle
        return &dense_[sparse_[e.index]];
    }

    void remove(Entity e) {
        if (!get(e)) return;
        std::uint32_t slot = sparse_[e.index];
        std::uint32_t last = static_cast<std::uint32_t>(dense_.size() - 1);

        dense_[slot] = std::move(dense_[last]);
        entities_[slot] = entities_[last];
        sparse_[entities_[slot].index] = slot;

        dense_.pop_back();
        entities_.pop_back();
        sparse_[e.index] = null_;
    }

    template<typename Fn>
    void each(Fn&& fn) {
        for (std::size_t i = 0; i < dense_.size(); ++i)
            fn(entities_[i], dense_[i]);
    }

private:
    static constexpr std::uint32_t null_ = ~0u;
    std::vector<std::uint32_t> sparse_;
    std::vector<T> dense_;
    std::vector<Entity> entities_;

    void ensureSparse(std::uint32_t index) {
        if (sparse_.size() <= index) sparse_.resize(index + 1, null_);
    }
};

class Registry {
public:
    Entity create() {
        if (!free_.empty()) {
            std::uint32_t idx = free_.back(); free_.pop_back();
            return { idx, generations_[idx] };
        }
        std::uint32_t idx = static_cast<std::uint32_t>(generations_.size());
        generations_.push_back(0);
        return { idx, 0 };
    }

    void destroy(Entity e) {
        transform_.remove(e);
        velocity_.remove(e);
        // ... or iterate registered pools
        ++generations_[e.index];
        free_.push_back(e.index);
    }

    ComponentPool<struct Transform>& transforms() { return transform_; }
    ComponentPool<struct Velocity>& velocities() { return velocity_; }

private:
    std::vector<std::uint32_t> generations_;
    std::vector<std::uint32_t> free_;
    ComponentPool<struct Transform> transform_;
    ComponentPool<struct Velocity> velocity_;
};

struct Transform { float x, y; };
struct Velocity  { float vx, vy; };

void movementSystem(Registry& world) {
    world.transforms().each([&](Entity e, Transform& t) {
        if (auto* v = world.velocities().get(e))
            t.x += v->vx, t.y += v->vy;
    });
}
```

**Why it’s better:** O(1) add/get/remove, fast iteration, no `EntityDestroyerFn`, no `std::map`.

References: [EnTT](https://github.com/skypjack/entt)

---

### 2. Archetype / chunk ECS — best for huge entity counts

Entities with the **same component signature** live in one memory chunk (SoA: arrays of `Transform`, `Velocity`, …). Systems query archetypes.

```cpp
struct Archetype {
    std::vector<Entity> entities;
    std::vector<Transform> transforms;
    std::vector<Velocity>  velocities;
};

class World {
    std::vector<Archetype> archetypes_;
public:
    void moveAll() {
        for (Archetype& a : archetypes_) {
            for (std::size_t i = 0; i < a.entities.size(); ++i) {
                a.transforms[i].x += a.velocities[i].vx;
                a.transforms[i].y += a.velocities[i].vy;
            }
        }
    }
};
```

Used by: flecs, Unity DOTS, Bevy (same idea in Rust).

---

### 3. Type-erased component registry — middle ground

One `Registry` owns `ComponentStorage*` keyed by type id. `destroy()` walks storages; no per-T destroyer function list in user code.

```cpp
struct IComponentStorage {
    virtual ~IComponentStorage() = default;
    virtual void remove(Entity e) = 0;
};

template<typename T>
struct ComponentStorage : IComponentStorage {
    ComponentPool<T> pool;
    void remove(Entity e) override { pool.remove(e); }
};

class Registry {
    std::unordered_map<std::type_index, std::unique_ptr<IComponentStorage>> storages_;

public:
    template<typename T>
    ComponentPool<T>& pool() {
        auto key = std::type_index(typeid(T));
        auto it = storages_.find(key);
        if (it == storages_.end()) {
            auto storage = std::make_unique<ComponentStorage<T>>();
            auto* ptr = &storage->pool;
            storages_.emplace(key, std::move(storage));
            return *ptr;
        }
        return static_cast<ComponentStorage<T>&>(*it->second).pool;
    }

    void destroy(Entity e) {
        for (auto& [_, storage] : storages_)
            storage->remove(e);
    }
};
```

---

## EnTT — production-ready library example

```cpp
#include <entt/entt.hpp>

struct Transform { float x, y, z; };
struct Velocity  { float vx, vy, vz; };

int main() {
    entt::registry registry;

    auto entity = registry.create();
    registry.emplace<Transform>(entity, 0.f, 0.f, 0.f);
    registry.emplace<Velocity>(entity, 1.f, 0.f, 0.f);

    registry.view<Transform, Velocity>().each([](Transform& t, Velocity& v) {
        t.x += v.vx;
        t.y += v.vy;
        t.z += v.vz;
    });

    registry.destroy(entity);
}
```

---

## Recommendation for gdevice

| Phase | Action |
|-------|--------|
| **Now** | Don’t extend `map` + `EntityDestroyerFn` much further |
| **P1 #4** | Move to **sparse-set pools** (hand-rolled or EnTT) |
| **Later** | Archetypes only if entity count / multithreading demands it |

Target API shape for terrain migration:

```cpp
Entity e = world.create();
world.emplace<Transform>(e, ...);
world.emplace<TerrainTile>(e, ...);

TerrainGenSystem(world);
TerrainRenderSystem(world);
world.destroy(e);
```

Suggested systems for clipmap migration:

- `LODSystem` — clipmap scroll, tile invalidation
- `TerrainGenSystem` — async compute dispatch
- `TerrainRenderSystem` — draw tiles with materials

---

## Approach comparison

| Approach | Best for |
|----------|----------|
| Map per type + destroyer fn (`entity.h` today) | Learning, tiny demos |
| Sparse set (EnTT) | Most C++ game/engine ECS needs |
| Archetype chunks | Huge entity counts, SIMD, multithreading |
| Type-erased registry | Tooling, plugins, dynamic component sets |

---

## Resume checklist (tomorrow)

- [ ] Read `include/type/scene/ecs/entity.h` and map current API to sparse-set equivalent
- [ ] Decide: **hand-rolled `gd::Registry`** vs. **vendoring EnTT**
- [ ] List components needed for `Heightmap → Clipmap → Tile` migration
- [ ] Prototype `ComponentPool<T>` in `include/type/scene/ecs/` (new header, keep old API until migration)
- [ ] Cross-link progress in `ROADMAP.md` P1 #4

**Related files:** `include/type/scene/ecs/entity.h`, `include/type/scene/ecs/ECS_test.h`, `include/type/scene/oop/node.h`
