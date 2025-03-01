# ShaderToy Ray Tracer for Random Spheres Scene
### This repository contains a Ray Tracing shader written for ShaderToy. The scene features multiple spheres of different materials, including diffuse, metal, and dielectric (glass). It implements a simple path-tracing algorithm to render the scene with realistic lighting and reflections.
## Link 
you can see the Shadertoy program from here :https://www.shadertoy.com/view/Wf23WV
## Features
- Ray Tracing: Computes the color of each pixel by simulating rays of light bouncing through the scene.
-  Material Types: Supports three material types:
    - Diffuse: Randomly scatters rays.
    - Metal: Reflects rays with a fuzz effect.
    - Dielectric (Glass): Refracts rays and reflects based on the Schlick approximation.
- Depth of Field: Implements a lens radius to simulate depth of field effects, blurring distant objects based on camera settings.
- Random Spheres: Includes a grid of small random spheres, with varied materials, placed around a few larger spheres.
- Anti-Aliasing: Uses jittering of pixel coordinates for smooth and high-quality rendering.
- Gamma Correction: Corrects the output for better visual representation on standard displays
