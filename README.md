# Racetrack generation (WIP)
Based on our paper "Generating Race Tracks With Repulsive Curves" https://ieeexplore.ieee.org/document/10645670

![overview_borders](https://github.com/user-attachments/assets/6e707ff1-99c5-46be-9a43-141380418ab5)

## Notes
1. The project is still in prototyping phase. Feel free to play around with it, but don't expect a polished user interface. Unexpected issues are expected.
2. There is no possibility yet to automatically generate a track from start to finish.
3. The lack of quality control in the intersection introduction can mess up tracks.
You can use the "Generate Bezier Spline" button anytime to revert to the initial spline with no intersectinos.
4. The generation process works in the Scene view, and only while the scene is redrawn.
To force Unity to constantly redraw the scene, switch to the "View Tool" (hand icon) and place your cursor inside the Scene view.

## Setup
1. The project has only one scene, "Main". Open it. "Main" contains three objects key to the generation:
"RoadObject" contains some wrapper MonoBehaviors to serialize mesh data even when you close the Unity project.
All other data in the tool is not stored here and will only be serialized within one session.
"CurveGenBackground" is a white background. Enable it during curve generation.
"3DModelBackground" is a green background. Enable it for the 3D model.
3. Open "Window > Map Gen Tool". This is the editor window which gives you control over the whole generation process.
4. Use the editor window to generate race tracks.

## Overview Map Gen Tool (editor window)
The Map Gen Tool currently features two sections:
1. "Advanced" exposes pretty much all internals of the generation process.
If you want to override the constrained space, add constraints or potentials, adjust the way crossings are formed (+ much more), this is where to look.
2. "Manual" houses a small subset of the features under "Advanced". In general, these are the only buttons you need for generating a curve and 3D model.
