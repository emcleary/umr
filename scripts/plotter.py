import sys
import os
import json
import vtkmodules.all as vtk
import matplotlib.pyplot as plt
from pathlib import Path
from multiprocessing import Process


def vtk_parser(filename):
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(filename)
    reader.Update()
    output = reader.GetOutput()

    points = []
    for i in range(output.GetNumberOfPoints()):
        pt = output.GetPoint(i)
        points.append((pt[0], pt[1]))

    cells = []
    edges = set()
    for i in range(output.GetNumberOfCells()):
        cell = output.GetCell(i)
        cells.append(cell)
        for j in range(cell.GetNumberOfEdges()):
            edge = cell.GetEdge(j)
            assert(edge.GetNumberOfPoints() == 2)
            pt0 = edge.GetPoints().GetPoint(0)[:2]
            pt1 = edge.GetPoints().GetPoint(1)[:2]
            if (pt0, pt1) not in edges:
                if (pt1, pt0) not in edges:
                    edges.add((pt0, pt1))

    
    source_points = []
    source_segments = []
    fields = output.GetFieldData()
    if fields.GetNumberOfArrays() == 2:
        assert fields.GetArrayName(0) == 'SourcePoints'
        assert fields.GetArrayName(1) == 'SourceSegments'
        data_sp = fields.GetArray(0)
        data_ss = fields.GetArray(1)
        for i in range(data_sp.GetDataSize()):
            v = data_sp.GetValue(i)
            source_points.append(points[v])
        for i in range(0, data_ss.GetDataSize(), 2):
            v = data_ss.GetValue(i)
            w = data_ss.GetValue(i+1)
            source_segments.append((points[v], points[w]))

    return points, list(edges), source_points, source_segments


def json_parser(filename):
    data = json.load(open(filename, 'r'))
    points = data['points']
    edges = data['edges']
    sources = data['source_points']
    segments = data['source_segments']
    for i, p in enumerate(points):
        points[i] = tuple(p)
    for i, e in enumerate(edges):
        edges[i] = (points[ei] for ei in e)
    for i, p in enumerate(sources):
        sources[i] = points[p]
    for i, e in enumerate(segments):
        segments[i] = (points[ei] for ei in e)
    return points, edges, sources, segments


def plotter(points, edges, source_points, segments, outfile=None, lw=1.0, ps=4.0, include_sources=False):
    fig = plt.figure(num=1, clear=True)
    ax = fig.add_subplot()

    ax.scatter(*list(zip(*points)), c='black', s=ps)

    x = []
    y = []
    for org, dest in edges:
        x += (org[0], dest[0])
        x.append(None)
        y += (org[1], dest[1])
        y.append(None)
    plt.plot(x, y, c='black', lw=lw)

    if include_sources:
        if source_points:
            ax.scatter(*list(zip(*source_points)), c='red', s=ps, zorder=10)

        x.clear()
        y.clear()
        for org, dest in segments:
            x += (org[0], dest[0])
            x.append(None)
            y += (org[1], dest[1])
            y.append(None)
        if x:
            plt.plot(x, y, c='red', lw=lw, zorder=10)

    if outfile:
        fig.tight_layout()
        plt.gca().set_aspect('equal')
        print('Saving plot as', outfile) 
        plt.savefig(outfile)
    else:
        plt.show()


def get_cmd_line_arguments():
    import argparse
    description = '2D plotting tool for VTK and JSON grid files'
    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('filename', nargs='+', help='input file, must have vtk extension')
    parser.add_argument('-w', '--line-width', default=1.0, type=float, help='plot line width')
    parser.add_argument('-p', '--point-size', default=4.0, type=float, help='plot point size')
    parser.add_argument('-d', '--display', action='store_true', help='show image rather than write to file')
    parser.add_argument('-n', '--processors', default=8, type=int, help='multiprocessor parallelization')
    parser.add_argument('-s', '--sources', action='store_true', help='include sources')
    
    return parser.parse_args()


def plot_files(*files):
    for filename in files:
        _, ext = os.path.splitext(filename)
        if ext == '.vtk':
            points, edges, sources, segments = vtk_parser(filename)
        elif ext == '.json':
            points, edges, sources, segments = json_parser(filename)
        else:
            print(f"Not setup for files with the {ext} extension")
            continue

        if not points:
            print('No points in', filename)
            continue
        outfile = None if args.display else str(Path(filename).with_suffix('.png'))
        plotter(points, edges, sources, segments, outfile,
                args.line_width, args.point_size, args.sources)


if __name__ == '__main__':
    args = get_cmd_line_arguments()
    filenames = args.filename
    nproc = args.processors
    nfiles = len(filenames) // nproc

    processes = []
    for i in range(nproc):
        if i == nproc-1:
            p = Process(target=plot_files, args=filenames[i*nfiles:])
        else:
            p = Process(target=plot_files, args=filenames[i*nfiles:nfiles*(i+1)])
        p.start()
        processes.append(p)

    for p in processes:
        p.join()
