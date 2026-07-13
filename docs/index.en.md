title: BEACH Documentation

Lang: [English](index.en.md) | [日本語](index.md)

# BEACH

BEACH simulates charge accumulation on insulating surfaces and charged-particle trajectories in the electric
field produced by that accumulated charge.

It evaluates surface charging, potential evolution, particle absorption, and escape on triangular
boundary-element meshes.

**[Start the 10-minute tutorial](Tutorial.en.html)** · [Installation](Installation.en.html)

## Run the shortest example

```bash
beachx config init beach.toml
beachx lint beach.toml
beach beach.toml
beachx inspect outputs/latest
```

Create and check a configuration, run the simulation, and inspect the results written to `outputs/latest`.

## From input to output

```text
Particles, surface mesh, and boundary conditions
                         ↓
                       BEACH
                         ↓
Surface charge, potential, particle statistics, and batch history
```

In each batch, BEACH computes the electric field, advances particles, and deposits the charge of absorbed
particles onto triangular surface elements. The accumulated charge contributes to the field in the next batch.

<div align="center">
  <img src="images/potential_history.gif" alt="Evolution of the potential distribution on an insulating mesh under electron-beam irradiation" width="80%">
  <p><i>Potential evolution on an insulating mesh under electron-beam irradiation</i></p>
  <sub>3D model: <a href="https://www.turbosquid.com/ja/3d-models/rubber-duck-pbr-game-ready-model-2001526">Rubber Duck PBR Game Ready</a> (TurboSquid)</sub>
</div>

Use the navigation sidebar for detailed usage and numerical-model documentation.
