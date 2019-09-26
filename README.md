# NatureWorks

![](https://img.shields.io/badge/minzipped_size-4.2MB-blue)
![](https://img.shields.io/badge/node-v10.14.1-yellow)
![](https://img.shields.io/badge/npm-6.4.1-yellow)
![](https://img.shields.io/badge/webpack-4.38.0-yellow)
![](https://img.shields.io/badge/three.js-r106-green)

**NatureWorks** is a web-based 3D Nature graphics engine for realistic real-time rendering of Nature, including sky, clouds, fire, flower, grass, ground, ocean, tree, smoke, snow, bubbles, sunlight, terrain, fire, water, a flock of birds, a school of fish and their combinations, based on the popular javascript library, [three.js](https://threejs.org/). It will be upgraded for more realistic rendering, and more features of nature will continue to be inserted for further analysis. You can create any natural features or phenomena you imagine in a way you wish through NatureWorks.

## Abstract
This application has a general purpose 3D rendering engine based on WebGL and JavaScript, which supports various functions such as undo/redo, auto saving, scene graph, multiple selections, three different controls (camera controls, transform controls, drag controls), two texture map libraries (matcap and normal map library), boolean operations (union, intersect, subtract), import/export files, and so on. In addition, it includes a general purpose CSS styling engine, which provides common user interface elements, such as drop-down menus, sidebar, context menu, tooltip, modal popups, drag & drop, and so on. It also has a structure for multiple languages and allows five themes (blue, dark, light, teal, transparent) using the selected language. Note that it was built on three.js and its associated js files.

# Highlighted Features

- Underlying Technologies<br>
WebGL (three.js), GLSL Shaders using ray-marching technique.

- Create from Template<br>
More than 15 different features of nature are provided for non-experts in JSON file format. For example, they are “sky, clouds, fire, flower, grass, ground, ocean, tree, smoke, snow, bubbles, sunlight, terrain, fire, water, a flock of birds, and a school of fish”.

- Drag & drop<br>
These files (listed above) can be uploaded to a web browser just in a drag & drop manner.

- Pop-up menu / context menu<br>
The position, size, number or properties of each feature can be edited through a simple GUI according to user's need. The main components of the GUI are the pop-up menus at the top of your web browser and the context menus that appear when selecting objects by right-clicking.

- Export or share<br>
All modified data can be stored in the JSON file format on a personal computer or shared with others via the Internet. 

- PBR & Shaders<br>
NatureWorks provides game developers with a framework for physically based rendering (PBR) and shader programming environment for procedural approach to render the huge amount of features in detail and in a realistic way.

- Shader Programming<br>
With developed shader programs, application developers are able to create a new scene just by calling the functions existing inside the shaders.


## User interface

- Undo / Redo
    - Ctrl + Z : Undo
    - Ctrl + Y : Redo
- Selections
    - Left Click : Single selection
    - Ctrl + Left Click : Multiple selection
    - Alt + Left Click : Ancestor selection of the clicked
    - Click Again : Cancel the selection
- Camera Controls
    - Left Click + Drag : Rotate
    - Scroll Wheel : Zoom
    - Right Click + Drag : Pan
- Transformation Controls
    - W key : Translate mode
    - E key : Rotate mode
    - R key : Scale mode
    - Q key : Toggle world / local space
    - X key : Toggle X-space
    - Y key : Toggle Y-space
    - Z key : Toggle Z-space
    - Space key: Toggle visible
    - '+' key : Increase axes
    - '-' key : Decrease axes
    - Hold Ctrl down : Translate or rotate snap
- Function Keys
    - F1 key : Toggle scene graph
    - F2 key : Toggle all GUIs
    - F3 key : Toggle drag controls
    - F4 key : Toggle transform controls
    - F5 key : Toggle all editors
    - F6 key : Screen capture
    - F8 key : Place the selected in screen center
    - F9 key : Clear the current setting
    - F10 key : Clear all data in current screen
- Keyboard Keys
    - Delete key : Delete the selected object
    - Escape key : Enable camera controls
    - C key : Toggle FPS / orbit camera
    - L key : Toggle all light helpers
    - M key : Toggle sounds if available
    - O key : Stop animation if available
    - P key : Play next animation if available

## Application Areas

NatureWorks can be used in various fields: scientific visualization, video games, movie making, GIS simulation, digital art, web design, virtual reality (VR), computer education, and other commercial works. For example, this engine will be used widely in making background for games played in natural sceneries. Additionally, it will provide the background scenes necessary for film that deal with disasters.

- Movie Making<br>
You can easily build background scenes that are expensive to produce.

- Digital Art<br>
You can easily design an amazing digital world that you imagine.

- Video Games<br>
You can create a variety of natural phenomena and render realistic images in a simple way.

- Virtual Reality (VR)<br>
You can design realistic natural phenomena in 3D space for various VR applications.


## Copyright Notice
This presentation is protected by U.S. and International copyright laws. Reproduction and distribution of the presentation without written permission of the sponsor is prohibited.<br>
Copyright (c) 2019 NovaGraphix, Co. All rights reserved.