# Fog Renderer

This project implements the fog renderer. To better understand this project, please visit the [website](https://58191554.github.io/fog_renderer.html) for it.

![pipelin1.png](https://github.com/58191554/FogRenderer/blob/main/images/pipeline1.png)

## Run Code

```
cd \out\build\x64-Release
.\pathtracer.exe -t 8 -s 16 -l 16 -m 6 -g -0.8  -m 8  ..\..\..\dae\sky\<dae file>
```

### Command line options

| Flag and parameters    | Description                                                  |
| ---------------------- | :----------------------------------------------------------- |
| `-s <INT>`             | Number of camera rays per pixel (default=1, should be a power of 2) |
| `-l <INT>`             | Number of samples per area light (default=1)                 |
| `-t <INT>`             | Number of render threads (default=1)                         |
| `-m <INT>`             | Maximum ray depth (default=1)                                |
| `-f <FILENAME>`        | Image (.png) file to save output to in windowless mode       |
| `-r <INT> <INT>`       | Width and height in pixels of output image (if windowless) or of GUI window |
| `-p <x> <y> <dx> <dy>` | Used with the -f flag (windowless mode) to render a cell with its upper left corner at [x,y] and spanning [dx, dy] pixels. |
| `-c <FILENAME>`        | Load camera settings file (mainly to set camera position when windowless) |
| `-a <INT> <FLOAT>`     | Samples per batch and tolerance for adaptive sampling        |
| `-H`                   | Enable hemisphere sampling for direct lighting               |
| `-h`                   | Print command line help message                              |

### Keyboard Commands

|   Key   | Action                                                    |
| :-----: | :-------------------------------------------------------- |
|    f    | Fog mode                                                  |
|    g    | Lens mode                                                 |
| `q`/`w` | Increase or decrease phase function parameter value in HG |
|    R    | Start rendering                                           |
|    S    | Save a screenshot                                         |
|  - / +  | Decrease/increase area light samples                      |
|  [ / ]  | Decrease/increase camera rays per pixel                   |
|  < / >  | Decrease/increase maximum ray depth                       |
|    H    | Toggle uniform hemisphere sampling                        |

### Rendered Image

![dif_phase.drawio (1)](https://github.com/58191554/FogRenderer/blob/main/images/dif_phase.drawio%20(1).png)
