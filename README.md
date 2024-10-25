# 2DIrregularPacking
---

# 排样程序说明

## 环境要求
程序运行需要安装 Python 环境，具体要求如下：
- **Python 版本**：`Python 3.10`
- **依赖包**：
  - `numpy`
  - `pandas`
  - `scipy`
  - `shapely`
  - `pyclipper`
  - `matplotlib`

依赖包可以通过以下命令安装：
```bash
pip install numpy pandas scipy shapely pyclipper matplotlib
```

## 一、整体说明

该排样程序用于优化多边形的布局，确保将多个零件在有限的容器中进行高效排样。程序支持基于 `ConvexHull` 或者 `MBR`（最小外接矩形）的排样方式，同时支持对多边形进行分组，在排样过程中基于BLF算法进行360°旋转的贪心排样。程序运行时，需提供输入文件路径以及相应的排样参数。

### 功能说明：
- **容器宽度设置**：用户可以设置排样容器的宽度。
- **零件与零件、零件与板材的距离**：可以设置零件间的最小间隔和与边界板材的最小间隔。
- **排样方式**：可以选择基于 `ConvexHull` 或 `MBR` 的排样方式。
- **分组功能**：用户可以选择是否对多边形进行分组处理。
- **结果输出**：排样结果将输出到用户指定的路径，支持输出为CSV格式。

## 二、数据准备说明

- 输入文件必须为 `.csv` 格式，且每个文件应包含多边形信息。
- 输入数据应包含每个多边形的顶点坐标，具体格式如下：
  - 每个多边形的数据列包含该类多边形的数量、以及以 JSON 格式存储的顶点坐标。
  - 示例：
    ```csv
    num,polygon
    4,"[[0,0],[10,0],[10,10],[0,10]]"
    3,"[[5,5],[10,5],[7.5,8]]"
    ```
- 根据江南方以往给定的数据格式，我们编写了dataprocess.py程序用以批量数据处理，直接将数据集中多组数据对应转化为`.csv`标准格式文件。文件位于./jndata/dataprocess.py


## 三、程序使用说明

程序的运行格式为通过命令行传递参数来执行排样操作，具体格式如下：

```bash
python pack.py --container_width <容器宽度> --part_gap <零件间距> --board_gap <板材间距> --packing_method <排样方式> --grouping <是否分组> --input_path <输入文件路径> --output_path <输出文件路径>
```

### 参数说明：
- `--container_width <容器宽度>`：设置排样容器的宽度，单位为毫米（mm），该参数为必填。
- `--part_gap <零件间距>`：设置零件与零件之间的最小间距，单位为毫米，默认值为 `15mm`。
- `--board_gap <板材间距>`：设置零件与容器板材之间的最小距离，单位为毫米，默认值为 `10mm`。
- `--packing_method <排样方式>`：设置排样方式，支持两种方式：
  - `ConvexHull`：基于凸包进行排样。
  - `MBR`：基于最小外接矩形（MBR）进行排样，默认选项。
- `--grouping`：是否对多边形进行分组处理，不传递该参数则默认为不分组，传递该参数则表示对多边形进行分组。
- `--input_path <输入文件路径>`：待排样数据的文件路径，必须是有效的 `.csv` 文件。
- `--output_path <输出文件路径>`：排样结果输出的文件路径，必须是有效的路径。

### 示例命令：

以下为一个完整的程序运行示例：
```bash
python pack.py --container_width 3000 --part_gap 15 --board_gap 10 --packing_method ConvexHull --grouping --input_path /home/user4/pack/jndata/Test-Data/41.csv --output_path ./res/
```

此命令将执行以下操作：
- 容器的宽度设置为 `3000mm`。
- 零件之间的间距为 `15mm`，零件与板材之间的间距为 `10mm`。
- 排样方式选择为 `ConvexHull`（基于凸包进行排样）。
- 使用了分组功能。
- 输入的待排样数据文件路径为 `/home/user4/pack/jndata/Test-Data/41.csv`。
- 排样结果将输出到当前目录的 `./res/` 文件夹。

### 输出结果：
- 程序将生成一个 `.csv` 文件，包含排样的详细结果，包括每个零件的位置坐标。还将生成一个`.png`文件，为排样结果样例图。

---