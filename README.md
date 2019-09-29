# NatureWorks

![](https://img.shields.io/badge/minzipped_size-4.2MB-blue)
![](https://img.shields.io/badge/node-v10.14.1-yellow)
![](https://img.shields.io/badge/npm-6.4.1-yellow)
![](https://img.shields.io/badge/webpack-4.38.0-yellow)
![](https://img.shields.io/badge/three.js-r106-green)

**NatureWorks** is a web-based 3D Nature graphics engine for a realistic real-time rendering of Nature, for example, sky, clouds, snows, rain, fish, bubbles, smoke, fire, tree, flowers, grass, water, ocean, ground, terrain, and their combinations, based on the popular javascript library, [three.js](https://threejs.org/). It will be upgraded for more realistic rendering, and more features of nature will continue to be inserted for more powerful modeling. You can create any natural features or phenomena you imagine in a way you wish through NatureWorks.

# Main functions

This application provides various functions such as undo/redo, auto saving, scene graph, multiple selections, different controls (camera controls, transform controls, and drag controls), texture map libraries (matcap material library & normal map library), boolean operations (union, intersect, and subtract), and importing 3D model files (11 file formats supported). In addition, it includes user interface elements such as drop-down menus, sidebar, context menu, tooltip, modal popups, and drag & drop. It also has a structure for multiple languages and allows five themes (blue, dark, light, teal, transparent).

# Highlighted features

- Supported natures<br>
This application provides 23 different kinds of natures for non-experts. They are sunlight, sky (sky dome), clouds (cloud dome), snows, rain, boids(e.g., a school of fish), curl noise, curl particles, bubbles (sprite bubbles), smoke, volume fire, fire particles, tree, tufts(e.g., flower, grass), water, ocean, ground, terrain, gpu particles, and raymarch-based natures (e.g., terrain, tunnel, galaxy, mountains, canyon, caves, sea, river, forest, fish, etc)”. Note that the raymarch-based natures were inspired by [Shadertoy](https://www.shadertoy.com).

- Undo / redo<br>
The undo function is used to erase the last change done to the document reverting it to an older state. And redo is used to rerun the recent actions you undid. Note that there is no limit to the number of undo operations.

- Auto save<br>
Like Google documents, all changed data is automatically saved. Note that if the data size exceeds the capacity of your web browser used, the storage may fail. In this case, you need to use the pop-up menu (File > Save As) to save the changed data.

- Scene graph<br>
In a scene graph with a tree structure, a node marked with a check box is a parent node having a child node, and a node without a check box is a leaf node. Each node has a colored box preceded by its name, and if it is a mesh, there is a blue node pointing to the geometry and a green node pointing to the material next to it. The selected node is displayed in red, and you can edit the details of that node through the right-click context menu.

- MatCap material library<br>
You can apply Matcap material to an object using the Matcap material library of 50 textures. Note that a custom material map can be also applied by dropping it to the last element of material drop zone.

- Normal map library<br>
You can adjust the normal map of material through the 26 normal textures and normal scale adjustments.

- Boolean operations
    - union: it obtains the union of two meshes.
    - intersect: it gets the intersection of two meshes.
    - subtract: tt obtains the difference set between two meshes.

- Context menu<br>
The context menus appears when selecting objects by mouse right-clicking. You can edit various properties of the 3D model using this context menu.

- Importing files<br>
You can load 3D model files using the pop-up menu (FILE > Open Files...) or drag & drop them to your web browser. Supported file formats are json, obj(+mtl), stl, ply, dae(collada), gltf, glb, amf, 3mf, wrl(vrml), and fbx.

## User interface

- Undo / redo
    - ctrl + z : undo
    - ctrl + y : redo
- Selection
    - mouse left click : single selection
    - mouse left click again : deselection of the clicked object
    - ctrl + mouse left click : multiple selection
    - alt + mouse left click : ancestor selection of the clicked object
- Camera controls
    - mouse left click + drag : rotate
    - mouse scroll wheel : zoom (in/out)
    - mouse right click + drag : pan
- Transformation controls
    - w key : translate mode
    - e key : rotate mode
    - r key : scale mode
    - q key : toggle world / local space
    - x key : toggle x-space
    - y key : toggle y-space
    - z key : toggle z-space
    - space key: toggle visible
    - '+' key : increase axes
    - '-' key : decrease axes
    - hold ctrl down : translate or rotate snap
- Function keys
    - F1 key : toggle scene graph
    - F2 key : toggle all GUIs
    - F3 key : toggle drag controls
    - F4 key : toggle transform controls
    - F5 key : Toggle all editors
    - F6 key : screen capture
    - F8 key : place the selected in screen center
    - F9 key : clear the current setting
    - F10 key : clear all data in current screen
- Keyboard keys
    - delete key : delete the selected object
    - escape key : enable camera controls
    - c key : toggle FPS / orbit camera
    - l key : toggle all light helpers
    - m key : toggle sounds if available
    - o key : stop animation if available
    - p key : play next animation if available

## Application areas

NatureWorks can be used in such applications as video games, movie making, GIS simulation, digital art, web design, virtual reality, and other commercial works. For example, this app is able to make background scenes for video games played in natural sceneries. Also it can produce complex backgrounds that are often required to make a film about disasters.

## Question or suggestion

Please contact us at <info@nova-graphix.com> for any question or suggestion.

Thank you for reading the above description on **NatureWorks**, developed by [NovaGraphix, Co.](https://www.nova-graphix.com/) Note that we will continue to add new features and technologies.